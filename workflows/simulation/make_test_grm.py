#!/usr/bin/env python3
import argparse,numpy as np,pandas as pd
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--data',required=True);a=p.parse_args();d=Path(a.data)
w=np.load(d/'standardized_genotypes.npy',mmap_mode='r'); A=(np.asarray(w,float)@np.asarray(w,float).T)/w.shape[1]; pref=d/'grm/full';pref.parent.mkdir(exist_ok=True)
A[np.tril_indices(len(A))].astype(np.float32).tofile(str(pref)+'.grm.bin');np.full(len(A)*(len(A)+1)//2,w.shape[1],np.float32).tofile(str(pref)+'.grm.N.bin')
pd.read_csv(d/'metadata.tsv',sep='\t')[['FID','IID']].to_csv(str(pref)+'.grm.id',sep='\t',header=False,index=False)

