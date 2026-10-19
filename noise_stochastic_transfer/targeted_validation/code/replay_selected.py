"""Replay two fixed 50-record chunks and compare to saved records."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,numpy as np
from pathlib import Path
from fresh_study import job
R=Path(__file__).resolve().parent.parent
meta=json.load(open(R/'results/fresh_meta.json'));saved=np.load(R/'results/fresh_records.npz');rows=[]
for i in [0,71]:
 _,_,a=job((meta['cells'][i],i,0,50));b=saved['cell'+str(i)][:50]
 ok=np.allclose(a,b,rtol=1e-11,atol=1e-11,equal_nan=True)
 rows.append(dict(cell=i,draws=50,allclose=bool(ok),max_abs_difference=float(np.nanmax(np.abs(a-b)))))
 assert ok,(i,rows[-1])
json.dump(rows,open(R/'results/replay_selected.json','w'),indent=2)
print(rows)
