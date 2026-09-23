#!/usr/bin/env python3
"""Observed profile-Hessian/delta-method SEs for an existing PredLMM fit."""
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from predlmm_ace.model import ExactProfileML, LowRankACE

def central_hessian(function, z, step):
    z = np.asarray(z, float); p = z.size; eye = np.eye(p) * step
    h = np.empty((p, p), float); f0 = float(function(z))
    for i in range(p):
        h[i, i] = (function(z + eye[i]) - 2*f0 + function(z - eye[i])) / step**2
        for j in range(i):
            value = (function(z+eye[i]+eye[j]) - function(z+eye[i]-eye[j])
                     - function(z-eye[i]+eye[j]) + function(z-eye[i]-eye[j])) / (4*step**2)
            h[i, j] = h[j, i] = value
    return (h + h.T) / 2

def read_keep(path, ids):
    if path is None:
        return np.ones(len(ids), dtype=bool)
    keep = pd.read_csv(path, sep=r"\s+", header=None, dtype=str)
    if keep.shape[1] == 1:
        wanted = set(keep.iloc[:, 0]); mask = ids.IID.isin(wanted).to_numpy()
    else:
        wanted = set(map(tuple, keep.iloc[:, :2].to_numpy()))
        mask = np.fromiter(((a,b) in wanted for a,b in ids[["FID","IID"]].to_numpy()), bool, len(ids))
    if mask.sum() != len(keep):
        raise ValueError(f"Keep rows={len(keep)} but matched factor samples={mask.sum()}")
    return mask

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fit-json", required=True)
    p.add_argument("--factor-dir", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--phenotype", required=True)
    p.add_argument("--family-column", required=True)
    p.add_argument("--continuous", default="")
    p.add_argument("--categorical", default="")
    p.add_argument("--keep", default=None)
    p.add_argument("--hessian-step", type=float, default=0.002)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    old = json.loads(Path(a.fit_json).read_text())
    model = old["model"]
    ids = pd.read_csv(Path(a.factor_dir)/"sample_ids.tsv", sep="\t", dtype=str)
    mask = read_keep(a.keep, ids)
    meta = pd.read_csv(a.metadata, sep="\t", dtype=str)
    data = ids.merge(meta, on=["FID","IID"], how="left", validate="one_to_one", indicator=True)
    if (data._merge != "both").any(): raise ValueError("Metadata missing for factor samples")
    data = data.drop(columns="_merge").loc[mask].reset_index(drop=True)
    continuous = [x for x in a.continuous.split(",") if x]
    categorical = [x for x in a.categorical.split(",") if x]
    y = pd.to_numeric(data[a.phenotype], errors="raise").to_numpy(float)
    design = pd.DataFrame({"intercept":np.ones(len(data))})
    for column in continuous:
        v = pd.to_numeric(data[column], errors="raise").to_numpy(float)
        design[column] = (v-v.mean())/v.std(ddof=0)
    if categorical:
        design = pd.concat([design, pd.get_dummies(data[categorical].astype(str), drop_first=True, dtype=float)], axis=1)
    x = design.to_numpy(float); labels = data[a.family_column].astype(str).to_numpy()
    b = np.load(Path(a.factor_dir)/"nystrom_factor.npy", mmap_mode="r")[mask]
    delta = np.load(Path(a.factor_dir)/"diagonal_correction.npy")[mask]
    engine = ExactProfileML(LowRankACE(b, delta, labels))
    rhs = np.column_stack((y, x)); n = len(y)

    if model == "AE":
        h2 = float(old["h2"]); e2 = 1-h2; zhat = np.array([np.log(h2/e2)])
        def unpack(z):
            h = 1/(1+np.exp(-float(z[0]))); return h, 0.0, 1-h
    else:
        h2, c2 = float(old["h2"]), float(old["c2"]); e2 = 1-h2-c2
        zhat = np.log([h2/e2, c2/e2])
        def unpack(z):
            ez=np.exp(np.r_[z,0]-np.max(np.r_[z,0])); q=ez/ez.sum(); return tuple(q)

    def objective(z):
        h,c,_ = unpack(z); solved,logdet=engine._solve_logdet(h,c,rhs)
        vy,vx=solved[:,0],solved[:,1:]; beta=np.linalg.solve(x.T@vx,x.T@vy)
        resid=y-x@beta; q=float(resid@(vy-vx@beta))
        return 0.5*(logdet+n*np.log(q/n))

    objective_recomputed = float(objective(zhat))
    objective_recorded = old.get("profile_objective")
    objective_difference = (
        None if objective_recorded is None
        else objective_recomputed - float(objective_recorded)
    )
    hessian = central_hessian(objective, zhat, a.hessian_step)
    eigen = np.linalg.eigvalsh(hessian); condition=float(np.linalg.cond(hessian))
    if eigen.min() <= 0: raise RuntimeError(f"Hessian not positive definite: eigenvalues={eigen}")
    covz=np.linalg.inv(hessian)
    if model == "AE":
        d=h2*(1-h2); jac=np.array([[d],[0.0],[d]])
    else:
        jac=np.array([[h2*(1-h2),-h2*c2],[-h2*c2,c2*(1-c2)],
                      [h2*e2,c2*e2]])
    cov=jac@covz@jac.T; se=np.sqrt(np.maximum(np.diag(cov),0))
    old.update({
        "se_h2":float(se[0]), "se_c2":None if model=="AE" else float(se[1]),
        "se_h2_plus_c2":float(se[2]), "covariance_proportions":cov.tolist(),
        "hessian_logratio":hessian.tolist(), "hessian_step":a.hessian_step,
        "hessian_condition_number":condition, "hessian_eigenvalues":eigen.tolist(),
        "se_method":"observed profile-likelihood Hessian; delta method",
        "profile_objective_recomputed":objective_recomputed,
        "profile_objective_difference":objective_difference,
        "source_fit_json":str(Path(a.fit_json).resolve()),
    })
    Path(a.out).write_text(json.dumps(old,indent=2)+"\n")
    print(json.dumps({k:old[k] for k in ["model","n","rank","h2","se_h2","c2","se_c2",
          "se_h2_plus_c2","hessian_condition_number","hessian_eigenvalues"]},indent=2))
if __name__ == "__main__": main()
