#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np, pandas as pd
def main():
 p=argparse.ArgumentParser(); p.add_argument('--metadata',required=True); p.add_argument('--rank',type=int,required=True); p.add_argument('--method',choices=['random','structure','family'],required=True); p.add_argument('--seed',type=int,required=True); p.add_argument('--out',required=True); a=p.parse_args()
 d=pd.read_csv(a.metadata,sep='\t',dtype={'FID':str}); rng=np.random.default_rng(a.seed); n=len(d); chosen=[]
 if a.method=='random': chosen=rng.choice(n,a.rank,replace=False)
 else:
  sizes=d.groupby('family_id').size(); fams=sizes[sizes>1].index.to_numpy();
  if a.method=='structure':
   strata=d.population.astype(str)
  else: strata=pd.Series('1',index=d.index)
  candidates=[]
  for fam in fams:
   idx=d.index[d.family_id==fam].to_numpy(); candidates.append(int(rng.choice(idx)))
  if len(candidates)>a.rank:
   by=pd.DataFrame({'idx':candidates,'stratum':strata.iloc[candidates].to_numpy()}); chosen=[]
   for _,grp in by.groupby('stratum'): chosen.extend(rng.choice(grp.idx,max(1,round(a.rank*len(grp)/len(by))),replace=False).tolist())
   chosen=np.array(chosen[:a.rank])
  else:
   chosen=list(candidates); remaining=np.setdiff1d(np.arange(n),chosen); need=a.rank-len(chosen)
   if need: chosen.extend(rng.choice(remaining,need,replace=False).tolist())
 chosen=np.sort(np.asarray(chosen,dtype=np.int64)); d.iloc[chosen][['FID','IID']].to_csv(a.out,sep='\t',index=False)
 print(f'COMPLETE method={a.method} rank={len(chosen)}',flush=True)
if __name__=='__main__': main()
