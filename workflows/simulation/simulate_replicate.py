#!/usr/bin/env python3
"""Simulate one replicate without ever forming an N x N GRM."""
import argparse, json, math
from pathlib import Path
import numpy as np, pandas as pd

def transmit(parent,rng): return rng.binomial(1,parent/2).astype(np.int8)

def randomized_pcs(w,k,seed,oversample=8):
    rng=np.random.default_rng(seed); n,m=w.shape; q=k+oversample
    omega=rng.normal(size=(m,q)); y=np.zeros((n,q))
    for s in range(0,n,2048): y[s:s+2048]=np.asarray(w[s:s+2048],float)@omega
    qmat,_=np.linalg.qr(y,mode='reduced'); b=np.zeros((q,m))
    for v in range(0,m,2048): b[:,v:v+2048]=qmat.T@np.asarray(w[:,v:v+2048],float)
    uh,_,_=np.linalg.svd(b,full_matrices=False); return qmat@uh[:,:k]

def write_plink(prefix,g,ids):
    n,m=g.shape; ids[['FID','IID']].assign(PAT=0,MAT=0,SEX=0,PHENO=-9).to_csv(str(prefix)+'.fam',sep='\t',header=False,index=False)
    with open(str(prefix)+'.bim','w') as f:
        for j in range(m): f.write(f'1\trs{j+1}\t0\t{j+1}\tA\tG\n')
    code=np.array([0,2,3],np.uint8); nb=(n+3)//4; padded=4*nb
    with open(str(prefix)+'.bed','wb') as f:
      f.write(bytes([0x6c,0x1b,0x01]))
      for j in range(m):
       x=code[np.asarray(g[:,j],int)]
       if padded>n:
        z=np.zeros(padded,np.uint8); z[:n]=x; x=z
       x=x.reshape(-1,4); f.write((x[:,0]|x[:,1]<<2|x[:,2]<<4|x[:,3]<<6).astype(np.uint8).tobytes())

def main():
 p=argparse.ArgumentParser(); p.add_argument('--simulation',type=int,required=True); p.add_argument('--N',type=int,required=True)
 p.add_argument('--M',type=int,required=True); p.add_argument('--h2',type=float,required=True); p.add_argument('--c2',type=float,required=True)
 p.add_argument('--replicate',type=int,required=True); p.add_argument('--seed',type=int,required=True); p.add_argument('--n-pcs',type=int,default=5); p.add_argument('--out',required=True); a=p.parse_args()
 out=Path(a.out); out.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(a.seed); n,m=a.N,a.M
 g=np.lib.format.open_memmap(out/'genotypes.npy',mode='w+',dtype='int8',shape=(n,m)); fids=[]; pops=[]
 if a.simulation==1:
  p0=rng.uniform(.1,.5,m); nf=n//10
  for f in range(nf):
   mo=rng.binomial(2,p0).astype(np.int8); fa=rng.binomial(2,p0).astype(np.int8)
   g[2*f]=transmit(mo,rng)+transmit(fa,rng); g[2*f+1]=transmit(mo,rng)+transmit(fa,rng); fids += [f'F{f:07d}']*2; pops += [1,1]
  for i in range(2*nf,n): g[i]=rng.binomial(2,p0); fids.append(f'U{i:07d}'); pops.append(1)
 else:
  if n%4: raise ValueError('Simulation 2 N must be divisible by 4')
  pbase=rng.uniform(.1,.5,m); theta=.15; alpha=pbase*(1-theta)/theta; beta=(1-pbase)*(1-theta)/theta
  pf=np.vstack([rng.beta(alpha,beta) for _ in range(4)]); cursor=0; per=n//4; fam_people=round(.30*per)
  if fam_people%3: raise ValueError('0.30*N/4 must be divisible by 3')
  for pop in range(4):
   for f in range(fam_people//3):
    mo=rng.binomial(2,pf[pop]).astype(np.int8); fa=rng.binomial(2,pf[pop]).astype(np.int8)
    fam=f'P{pop+1}F{f:06d}'; g[cursor]=mo; g[cursor+1]=transmit(mo,rng)+transmit(fa,rng); g[cursor+2]=transmit(mo,rng)+transmit(fa,rng)
    fids += [fam]*3; pops += [pop+1]*3; cursor+=3
   for u in range(per-fam_people): g[cursor]=rng.binomial(2,pf[pop]); fids.append(f'P{pop+1}U{u:07d}'); pops.append(pop+1); cursor+=1
 g.flush(); ids=pd.DataFrame({'FID':fids,'IID':[f'I{i+1:08d}' for i in range(n)],'family_id':fids,'population':pops})
 mean=np.zeros(m); var=np.zeros(m)
 for v in range(0,m,2048):
  x=np.asarray(g[:,v:v+2048],float); mean[v:v+2048]=x.mean(0); var[v:v+2048]=x.var(0)
 sd=np.sqrt(np.maximum(var,1e-8)); w=np.lib.format.open_memmap(out/'standardized_genotypes.npy',mode='w+',dtype='float32',shape=(n,m))
 for s in range(0,n,1024): w[s:s+1024]=((np.asarray(g[s:s+1024],float)-mean)/sd).astype(np.float32)
 w.flush(); pcs=randomized_pcs(w,a.n_pcs,a.seed+1)
 causal=rng.choice(m,max(1,round(.1*m)),replace=False); raw=np.asarray(w[:,causal],float)@rng.normal(size=len(causal)); X=np.column_stack([np.ones(n),pcs]); raw-=X@np.linalg.lstsq(X,raw,rcond=None)[0]
 gv=raw*math.sqrt(a.h2/np.var(raw)); fam_eff={f:rng.normal(0,math.sqrt(a.c2)) for f in pd.unique(np.asarray(fids, dtype=object))} if a.c2 else {}
 cv=np.array([fam_eff.get(f,0.) for f in fids]); y=gv+cv+rng.normal(0,math.sqrt(1-a.h2-a.c2),n)
 ids['phenotype']=y
 for j in range(a.n_pcs): ids[f'PC{j+1}']=pcs[:,j]
 ids.to_csv(out/'metadata.tsv',sep='\t',index=False); ids[['FID','IID','phenotype']].to_csv(out/'phenotype.txt',sep='\t',header=False,index=False)
 q=ids[['FID','IID']+[f'PC{i+1}' for i in range(a.n_pcs)]]; q.to_csv(out/'pcs.qcovar',sep='\t',header=False,index=False)
 sizes=ids.groupby('FID').size(); ids.loc[ids.FID.map(sizes)==1,['FID','IID']].to_csv(out/'unrelated.keep',sep='\t',header=False,index=False)
 write_plink(out/'geno',g,ids); (out/'complete.json').write_text(json.dumps(vars(a),indent=2)+'\n'); print(f'COMPLETE simulation={a.simulation} N={n} M={m}',flush=True)
if __name__=='__main__': main()
