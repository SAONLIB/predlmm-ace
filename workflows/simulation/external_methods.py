#!/usr/bin/env python3
import argparse, json, re, subprocess, time
from pathlib import Path
import numpy as np
def run(cmd,log):
 t=time.time()
 with open(log,'w') as f: rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
 return rc,time.time()-t
def hsq(path):
 h=se=np.nan
 if Path(path).exists():
  for line in Path(path).read_text(errors='replace').splitlines():
   x=line.split()
   if len(x)>=3 and x[0] in {'V(G)/Vp','V(1)/Vp'}: h,se=map(float,x[1:3])
 return h,se
def fast(path):
 if not Path(path).exists(): return np.nan
 text=Path(path).read_text(errors='replace'); num=r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)'
 z=re.findall(rf'Heritability\s*=\s*{num}',text,re.I)
 return float(z[-1]) if z else np.nan
def main():
 p=argparse.ArgumentParser(); p.add_argument('--data',required=True); p.add_argument('--grm-prefix',required=True); p.add_argument('--simulation',type=int,required=True); p.add_argument('--gcta',required=True); p.add_argument('--threads',type=int,default=8); p.add_argument('--cutoffs',required=True); p.add_argument('--out',required=True); a=p.parse_args()
 d=Path(a.data); out=Path(a.out); out.mkdir(parents=True,exist_ok=True); grm=Path(a.grm_prefix); b=d/'geno'; t0=0.0
 gout=out/'gcta'; cmd=[a.gcta,'--reml','--grm',str(grm),'--pheno',str(d/'phenotype.txt'),'--qcovar',str(d/'pcs.qcovar'),'--thread-num',str(a.threads),'--out',str(gout)]
 if a.simulation==1: cmd += ['--keep',str(d/'unrelated.keep')]
 rcg,tg=run(cmd,out/'gcta.driver.log'); gh,gs=hsq(str(gout)+'.hsq'); rows=[dict(method='GCTA-GREML',estimate_h2=gh,se_h2=gs,runtime_sec=t0+tg,success=rcg==0)]
 for c in [float(x) for x in a.cutoffs.split(',')]:
  sp=out/f'sparse_{c:g}'; rc1,t1=run([a.gcta,'--grm',str(grm),'--make-bK-sparse',str(c),'--thread-num',str(a.threads),'--out',str(sp)],out/f'sparse_{c:g}.log')
  fh=np.nan; rc2=-1; t2=0
  if rc1==0:
   fo=out/f'fastgwa_{c:g}'; rc2,t2=run([a.gcta,'--bfile',str(b),'--grm-sparse',str(sp),'--fastGWA-mlm','--model-only','--est-vg','REML','--pheno',str(d/'phenotype.txt'),'--qcovar',str(d/'pcs.qcovar'),'--thread-num',str(a.threads),'--out',str(fo)],out/f'fastgwa_{c:g}.driver.log'); fh=fast(str(fo)+'.log')
  rows.append(dict(method='fastGWA-REML',sparse_cutoff=c,estimate_h2=fh,runtime_sec=t1+t2,success=rc1==0 and rc2==0))
 import pandas as pd; pd.DataFrame(rows).to_csv(out/'results.csv',index=False); print('COMPLETE external',flush=True)
if __name__=='__main__': main()
