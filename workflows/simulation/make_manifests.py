#!/usr/bin/env python3
import argparse, csv
from pathlib import Path

def split(x, cast): return [cast(v) for v in x.split(',') if v]
def write(path, rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(rows)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--reps',type=int,required=True)
    p.add_argument('--fractions',required=True); p.add_argument('--methods',required=True)
    p.add_argument('--out',required=True); p.add_argument('--include-n4000',action='store_true'); a=p.parse_args()
    ns=[4000,8000,12000] if a.include_n4000 else [8000,12000]
    reps=[]; rid=0
    for sim in (1,2):
      for n in ns:
       for h in (.25,.70):
        for c in ((.10,.25) if sim==1 else (0.0,)):
         for rep in range(1,a.reps+1):
          reps.append(dict(task_id=rid,simulation=sim,N=n,h2=h,c2=c,replicate=rep)); rid+=1
    factors=[]; fid=0
    for r in reps:
      for frac in split(a.fractions,float):
       for method in split(a.methods,str):
        factors.append(dict(factor_id=fid,replicate_task_id=r['task_id'],simulation=r['simulation'],N=r['N'],
                            h2=r['h2'],c2=r['c2'],replicate=r['replicate'],knot_fraction=frac,
                            rank=max(2,int(r['N']*frac)),knot_method=method)); fid+=1
    root=Path(a.out); write(root/'replicates.tsv',reps); write(root/'factors.tsv',factors)
    print(f'replicates={len(reps)} factors={len(factors)}')
if __name__=='__main__': main()

