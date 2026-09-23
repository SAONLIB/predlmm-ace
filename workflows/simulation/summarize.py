#!/usr/bin/env python3
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
def main():
 p=argparse.ArgumentParser();p.add_argument('--root',required=True);a=p.parse_args();root=Path(a.root)
 reps=pd.read_csv(root/'manifests/replicates.tsv',sep='\t'); fac=pd.read_csv(root/'manifests/factors.tsv',sep='\t'); rows=[]
 for _,f in fac.iterrows():
  path=root/f"results/predlmm_{int(f.factor_id):06d}_se.json"; base=root/f"results/predlmm_{int(f.factor_id):06d}.json"
  if not path.exists(): path=base
  if path.exists():
   z=json.loads(path.read_text()); rows.append({**f.to_dict(),'method':'PredLMM-ACE','estimate_h2':z['h2'],'estimate_c2':z['c2'],'se_h2':z.get('se_h2'),'se_c2':z.get('se_c2'),'runtime_sec':z.get('runtime_sec'),'success':z.get('converged',False)})
 for _,r in reps.iterrows():
  path=root/f"external/task_{int(r.task_id):06d}/results.csv"
  if path.exists():
   for _,z in pd.read_csv(path).iterrows(): rows.append({**r.to_dict(),**z.to_dict()})
 d=pd.DataFrame(rows); d.to_csv(root/'all_results.csv',index=False)
 keys=['simulation','N','h2','c2','method','knot_method','knot_fraction','sparse_cutoff']; [d.__setitem__(k,np.nan) for k in keys if k not in d]
 out=[]
 for key,x in d.groupby(keys,dropna=False):
  est=pd.to_numeric(x.estimate_h2,errors='coerce'); truth=float(x.h2.iloc[0]); se=pd.to_numeric(x.get('se_h2'),errors='coerce')
  out.append({**dict(zip(keys,key)),'n':len(x),'mean_h2':est.mean(),'bias_h2':(est-truth).mean(),'rmse_h2':np.sqrt(((est-truth)**2).mean()),'empirical_sd_h2':est.std(),'mean_se_h2':se.mean(),'coverage_h2':(((est-1.96*se)<=truth)&(truth<=(est+1.96*se))).mean()})
 pd.DataFrame(out).to_csv(root/'summary.csv',index=False); print(f'COMPLETE rows={len(d)}',flush=True)
if __name__=='__main__':main()
