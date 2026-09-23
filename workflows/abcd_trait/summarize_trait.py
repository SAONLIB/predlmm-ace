#!/usr/bin/env python3
import argparse, json, math, re
from pathlib import Path
import pandas as pd
from scipy.stats import chi2

p=argparse.ArgumentParser(); p.add_argument("--trait",required=True); p.add_argument("--out-dir",required=True)
a=p.parse_args(); root=Path(a.out_dir); res=root/"results"; rows=[]

def ci(est,se):
    if est is None or se is None: return (None,None)
    return (max(0,est-1.96*se),min(1,est+1.96*se))

for model,sample,rank in [("ACE","full",4000),("ACE","full",5000),("AE","full",5000),("AE","onefamily",5000)]:
    f=res/f"{a.trait}_PredLMM_{model}_{sample}_rank{rank}_se.json"
    if not f.exists(): continue
    x=json.loads(f.read_text()); lo,hi=ci(x.get("h2"),x.get("se_h2"))
    clo,chi=ci(x.get("c2"),x.get("se_c2"))
    rows.append(dict(method=f"PredLMM-{model}",sample=sample,rank=rank,n=x.get("n"),
      h2=x.get("h2"),se_h2=x.get("se_h2"),h2_ci_low=lo,h2_ci_high=hi,
      c2=x.get("c2"),se_c2=x.get("se_c2"),c2_ci_low=clo,c2_ci_high=chi,
      converged=x.get("converged"),objective=x.get("profile_objective")))

hsq=res/f"{a.trait}_GCTA_onefamily.hsq"
if hsq.exists():
    tab=pd.read_csv(hsq,sep="\t"); hit=tab[tab.Source.eq("V(G)/Vp")].iloc[0]
    rows.append(dict(method="GCTA-GREML",sample="onefamily",rank=None,n=None,h2=hit.Variance,
      se_h2=hit.SE,h2_ci_low=max(0,hit.Variance-1.96*hit.SE),h2_ci_high=min(1,hit.Variance+1.96*hit.SE),
      c2=None,se_c2=None,c2_ci_low=None,c2_ci_high=None,converged=True,objective=None))

for log in sorted(res.glob(f"{a.trait}_fastGWA_*.log")):
    text=log.read_text(errors="replace"); m=re.search(r"Heritability\s*=\s*([0-9.eE+-]+)",text)
    h=float(m.group(1)) if m else None
    # The log does not directly report the delta-method SE of h2. Do not use
    # SE(Vg)/Vp because it ignores uncertainty in Ve and their covariance.
    seh=None
    lo,hi=ci(h,seh); label=log.stem.rsplit("_",1)[-1]
    rows.append(dict(method=f"fastGWA-{label}",sample="full",rank=None,n=None,h2=h,se_h2=seh,
      h2_ci_low=lo,h2_ci_high=hi,c2=None,se_c2=None,c2_ci_low=None,c2_ci_high=None,
      converged="fastGWA-REML converged" in text,objective=None))

df=pd.DataFrame(rows)
df.to_csv(res/f"{a.trait}_final_method_comparison.csv",index=False)

# Same-rank, same-sample ACE-vs-AE boundary LRT. Objectives are negative profile log likelihoods.
try:
    ace=df[(df.method=="PredLMM-ACE")&(df["rank"]==5000)&(df["sample"]=="full")].iloc[0]
    ae=df[(df.method=="PredLMM-AE")&(df["rank"]==5000)&(df["sample"]=="full")].iloc[0]
    lrt=2*(ae.objective-ace.objective); pval=.5*chi2.sf(max(lrt,0),1)
    pd.DataFrame([{"comparison":"ACE versus AE, full sample, rank 5000","LRT":lrt,"boundary_p":pval}]).to_csv(
      res/f"{a.trait}_ACE_AE_LRT.csv",index=False)
except Exception as e: print("LRT unavailable:",e)

show=df.copy()
for c in ["h2","se_h2","h2_ci_low","h2_ci_high","c2","se_c2","c2_ci_low","c2_ci_high"]:
    if c in show: show[c]=show[c].map(lambda z:"" if pd.isna(z) else f"{z:.3f}")
cols=["method","sample","rank","n","h2","se_h2","h2_ci_low","h2_ci_high","c2","se_c2","c2_ci_low","c2_ci_high"]
(res/f"{a.trait}_final_method_comparison.tex").write_text(show[cols].to_latex(index=False,escape=True,na_rep=""))
print(show[cols].to_string(index=False))
