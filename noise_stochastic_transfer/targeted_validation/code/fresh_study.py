"""The 72 fixed fresh cells; 50 draws per independently seeded chunk."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import sys,json,time
from pathlib import Path
import numpy as np
from multiprocessing import Pool
from backend import *
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;CFG=json.load(open(ROOT/'config.json'))
def job(task):
 cell,index,chunk,nd=task;r=CL[cell['cluster']][2];h0=mean(cell['cluster'],cell['z'],cell['lambda']);x=(h0+h0[::-1])/2+1j*(h0-h0[::-1])/2;sigma=cell['sigma'];rng=np.random.default_rng(np.random.SeedSequence([CFG['master_seed'],index,chunk]));rows=[]
 assert np.max(np.abs(x-np.conj(x[::-1])))<1e-13
 for k in range(nd):
  xi=(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2);y=x+sigma*xi
  yh=(y+np.conj(y[::-1]))/2;ya=(y-np.conj(y[::-1]))/(2j);h=yh.real+yh.imag;aa=ya.real+ya.imag;sh=np.sqrt(2*(aa@aa)/N)
  try:
   s=statistics(h,sh,r)
   if s['status']<=0:raise RuntimeError('solver status nonpositive')
   row=[s['tL'],s['TC'],s['A_hat'],s['m0_hat'],sh,s['Hcond'],s['status'],s['nfev'],s['optimality']]
  except (RuntimeError,ValueError,np.linalg.LinAlgError,FloatingPointError) as e:
   row=[0.,0.,-1.,1.,sh,np.nan,-1.,0.,np.nan]
  rows.append(row)
 return index,chunk,np.array(rows)
def main():
 roots=json.load(open(ROOT/'results/roots.json'));rd={(a['cluster'],a['z']):float(a['lambda_c']) for a in roots};cells=[]
 for name in CFG['clusters']:
  r=CL[name][2];hh=CL[name][3];c=CONST[r]
  for z in CFG['spacings']:
   for nu in CFG['noise_levels']:
    for position in CFG['positions']:
     rho=nu*c['unit'];A=hh*z**(2*r);width=rho*c['sX']/c['Phi1'];lc=rd[name,z]
     cells.append(dict(cluster=name,z=z,nu=nu,position=position,lambda_c=lc,lambda_=lc+position*width,truth='A' if position<0 else 'B',A=A,rho=rho,sigma=rho*A,unit=c['unit'],width=width))
 for c in cells:c['lambda']=c.pop('lambda_')
 json.dump(dict(cells=cells,columns=['tL','TC','A_hat','m0_hat','sigma_hat','Hcond','status','nfev','optimality']),open(ROOT/'results/fresh_meta.json','w'),indent=2)
 tasks=[]
 for i,c in enumerate(cells):
  for ch in range(CFG['draws_per_cell']//50):tasks.append((c,i,ch,50))
 if '--smoke' in sys.argv:tasks=[(cells[0],0,0,2),(cells[-1],len(cells)-1,0,2)]
 t=time.time();save={};completed=0
 with Pool(CFG['workers']) as pool:
  for i,ch,a in pool.imap_unordered(job,tasks,chunksize=1):
   save[i,ch]=a;completed+=1
   if completed%12==0 or '--smoke' in sys.argv:print('chunks',completed,'/',len(tasks),'elapsed',time.time()-t,'failed',sum(int((v[:,6]<0).sum()) for v in save.values()),flush=True)
 arrays={}
 for i in range(len(cells)):
  aa=[save[i,ch] for ch in sorted(ch for ii,ch in save if ii==i)]
  if aa:arrays['cell'+str(i)]=np.vstack(aa)
 fn='fresh_smoke_DISCARD.npz' if '--smoke' in sys.argv else 'fresh_records.npz'
 np.savez_compressed(ROOT/'results'/fn,**arrays)
 print('complete',time.time()-t,flush=True)
if __name__=='__main__':main()
