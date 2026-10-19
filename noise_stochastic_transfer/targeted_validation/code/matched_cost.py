"""Fixed-grid operating curves with paired, within-cell bootstrap."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import json,csv,numpy as np
from scipy.stats import t as student
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;CFG=json.load(open(ROOT/'config.json'))
AL=np.array(CFG['alphas']);Q=student.ppf(1-AL,21);B=CFG['bootstrap_replicates']
def putcsv(name,rows):
 if not rows:return
 with open(ROOT/'results'/name,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def wilson(k,n):
 if n==0:return None,None
 z=1.95996398454;p=k/n;den=1+z*z/n;c=(p+z*z/(2*n))/den;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return max(0,c-h),min(1,c+h)
def score(tL,TC,sel,truth):
 chosen=np.where(sel,TC,tL);agree=np.where(tL*TC>0,np.sign(tL)*np.minimum(abs(tL),abs(TC)),0.)
 out={}
 for rule,st in [('S',chosen),('I',agree),('C',TC),('L',tL)]:
  d=np.sign(st)*(1 if truth=='B' else -1);d=np.where(np.isfinite(st),d,0)
  passed=np.abs(np.nan_to_num(st,nan=0.))[:,None]>Q[None,:]
  out[rule]=dict(correct=(d[:,None]>0)&passed,wrong=(d[:,None]<0)&passed,d=d,bins=passed.sum(1))
 return out

def interp(w,c,target):
 x=[];y=[]
 for xx,yy in zip(w,c):
  if x and xx==x[-1]:y[-1]=max(y[-1],yy)
  else:x.append(xx);y.append(yy)
 if len(x)<2 or target<x[0] or target>x[-1]:return np.nan
 return float(np.interp(target,x,y))

def analyze(cells,name):
 rng=np.random.default_rng(np.random.SeedSequence([2026101910,0 if name=='old' else 1]))
 sums={};boots={};sizes={};individual=[]
 for ci,c in enumerate(cells):
  out=score(c['tL'],c['TC'],c['sel'],c['truth']);n=len(c['tL'])
  groups=['all','nu='+str(float(c['nu']))]
  if name=='fresh':groups+=['z='+str(float(c['z'])),'case='+c['cluster']+'@'+str(float(c['z']))]
  for g in groups:
   sums.setdefault(g,{r:np.zeros((len(AL),2)) for r in out});boots.setdefault(g,np.zeros((B,2,len(AL),2)));sizes[g]=sizes.get(g,0)+n
   for r in out:sums[g][r]+=np.column_stack([out[r]['correct'].sum(0),out[r]['wrong'].sum(0)])
  joint=np.column_stack([out['S']['d'],out['S']['bins'],out['C']['d'],out['C']['bins']]).astype(int)
  cat,cnt=np.unique(joint,axis=0,return_counts=True);weights=rng.multinomial(n,cnt/n,size=B)
  vals=np.zeros((len(cat),2,len(AL),2))
  for j,row in enumerate(cat):
   for k in range(2):
    sign,bin_=row[2*k:2*k+2];passed=np.arange(len(AL))>=len(AL)-bin_
    vals[j,k,:,0]=(sign>0)*passed;vals[j,k,:,1]=(sign<0)*passed
  bs=(weights@vals.reshape(len(cat),-1)).reshape(B,2,len(AL),2)
  for g in groups:boots[g]+=bs
  ia=int(np.where(AL==CFG['alpha'])[0][0])
  for r in out:
   co=int(out[r]['correct'][:,ia].sum());wr=int(out[r]['wrong'][:,ia].sum());lo,hi=wilson(wr,n)
   individual.append(dict(cluster=c['cluster'],z=c['z'],nu=c['nu'],position=c['position'],rule=r,n=n,correct=co,wrong=wr,abstain=n-co-wr,wrong_rate=wr/n,wrong_low=lo,wrong_high=hi,conditional_wrong=wr/(co+wr) if co+wr else None,selected_C=float(np.mean(c['sel'])),failures=int(c.get('failures',0))))
 curves=[];matched=[]
 for g in sums:
  M=sizes[g];boot=boots[g]/M
  for rule in sums[g]:
   rates=sums[g][rule]/M
   for i,a in enumerate(AL):curves.append(dict(group=g,rule=rule,alpha=float(a),n=M,correct=float(rates[i,0]),wrong=float(rates[i,1]),abstain=float(1-rates[i].sum()),conditional_wrong=float(rates[i,1]/rates[i].sum()) if rates[i].sum() else None))
  for target in CFG['matched_errors']:
   base={r:interp(sums[g][r][:,1]/M,sums[g][r][:,0]/M,target) for r in ['S','C']}
   differences=[]
   for b in range(B):
    a=interp(boot[b,0,:,1],boot[b,0,:,0],target);i=interp(boot[b,1,:,1],boot[b,1,:,0],target)
    if np.isfinite(a+i):differences.append(a-i)
   available=len(differences)>=.95*B
   lo,hi=np.quantile(differences,[.025,.975]) if available else (np.nan,np.nan)
   matched.append(dict(group=g,target_wrong=target,S_correct=base['S'],C_correct=base['C'],difference=base['S']-base['C'],low=float(lo),high=float(hi),bootstrap_attainable=len(differences),bootstrap_total=B))
 putcsv(name+'_matched_cost.csv',matched)
 summary=[r for r in matched if r['group']=='all'];print(name,'matched pooled',summary,flush=True)
 return individual

def loadold():
 D=json.load(open(ROOT.parent/'adaptive_selection/results/raw.json'))
 cells=[]
 for c in D['configs']:
  rr=c['records'];cells.append(dict(cluster=c['cluster'],z=c['z'],nu=c['nu'],position=c['pos'],truth=c['truth'],tL=np.array([a['tL'] for a in rr]),TC=np.array([a['TC'] for a in rr]),sel=np.array([a['sel']=='C' for a in rr])))
 return cells

def loadfresh():
 meta=json.load(open(ROOT/'results/fresh_meta.json'));D=np.load(ROOT/'results/fresh_records.npz');cells=[]
 for i,c in enumerate(meta['cells']):
  arr=D['cell'+str(i)];r=1 if c['cluster']=='pair' else 2;A=arr[:,2];m0=arr[:,3];sig=arr[:,4];unit=c['unit'];ok=(A>0)&(m0>0)
  sel=np.zeros(len(arr),bool);sel[ok]=sig[ok]/(unit*A[ok])<=(A[ok]/m0[ok])**(.25/r)
  cells.append(dict(cluster=c['cluster'],z=c['z'],nu=c['nu'],position=c['position'],truth=c['truth'],tL=arr[:,0],TC=arr[:,1],sel=sel,failures=int((arr[:,6]<=0).sum())))
 return cells

def sensitivity():
 meta=json.load(open(ROOT/'results/fresh_meta.json'));D=np.load(ROOT/'results/fresh_records.npz');rows=[];q=student.ppf(.95,21)
 for i,c in enumerate(meta['cells']):
  a=D['cell'+str(i)];r=1 if c['cluster']=='pair' else 2;ok=(a[:,2]>0)&(a[:,3]>0)
  for fraction in CFG['gamma_fractions']:
   for mult in CFG['selector_multipliers']:
    sel=np.zeros(len(a),bool);sel[ok]=a[ok,4]/(c['unit']*a[ok,2])<=mult*(a[ok,2]/a[ok,3])**(fraction/r)
    t=np.where(sel,a[:,1],a[:,0]);d=np.where(abs(t)>q,np.sign(t),0)*(1 if c['truth']=='B' else -1)
    co=int((d>0).sum());wr=int((d<0).sum());lo,hi=wilson(wr,len(a))
    rows.append(dict(cluster=c['cluster'],z=c['z'],nu=c['nu'],position=c['position'],gamma_fraction=fraction,multiplier=mult,n=len(a),correct=co,wrong=wr,abstain=len(a)-co-wr,wrong_rate=wr/len(a),wrong_low=lo,wrong_high=hi,selected_C=float(sel.mean())))
 putcsv('sensitivity.csv',rows)

if __name__=='__main__':
 import sys
 if sys.argv[1]=='old':analyze(loadold(),'old')
 else:analyze(loadfresh(),'fresh');sensitivity()
