import json,csv,numpy as np
from pathlib import Path
from backend import *
ROOT=Path(__file__).resolve().parent.parent;meta=json.load(open(ROOT/'results/fresh_meta.json'));rows=[]
for i,c in enumerate(meta['cells']):
 name,z,lam=c['cluster'],c['z'],c['lambda'];w,o,r,hh=CL[name];h=mean(name,z,lam);x=(h+h[::-1])/2+1j*(h-h[::-1])/2
 original=(np.array(w)[None,:]*np.exp(1j*z*np.outer(S,o))).sum(1)-lam*hh*z**(2*r)*np.exp(1j*B*S)
 model_error=float(np.max(abs(x-original)));anti=float(np.linalg.norm((x-np.conj(x[::-1]))/(2j)))
 reg=regression(h,r);a,b=VP(h,r).wells(reg['nodes']);d=a['cost']-b['cost'];label='A' if d<0 else 'B'
 rows.append(dict(cell=i,cluster=name,z=z,nu=c['nu'],position=c['position'],truth=c['truth'],cost_label=label,cost_difference=d,mean_max_error=model_error,anti_mean_norm=anti))
with open(ROOT/'results/fresh_mean_checks.csv','w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
print('Noiseless label mismatches',sum(a['truth']!=a['cost_label'] for a in rows),'max model error',max(a['mean_max_error'] for a in rows),'max anti mean',max(a['anti_mean_norm'] for a in rows),flush=True)
