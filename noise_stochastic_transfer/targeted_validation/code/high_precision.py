"""Independent direct 21-sample regression, high-precision stationary wells."""
import sys,os,json,time
from pathlib import Path
import mpmath as mp
import numpy as np
from backend import *
HERE=Path(__file__).resolve().parent;OUT=HERE.parent/'results'
def direct(name,z,lam,r,y=None,dps=70):
 mp.mp.dps=dps
 zz=mp.mpf(str(z));w,o,_,hh=CL[name];ss=[mp.mpf(j)/21 for j in range(-10,11)]
 mu=[sum(mp.mpf(str(a))*mp.mpf(str(v))**n for a,v in zip(w,o)) for n in range(2*r+1)]
 H0=mp.matrix([[mu[i+j] for j in range(r)] for i in range(r)]);bv0=mp.matrix(mu[r:2*r]);hr=mu[-1]-(bv0.T*mp.lu_solve(H0,bv0))[0];A=hr*zz**(2*r)
 def atom(v):return mp.matrix([mp.cos(v*t)+mp.sin(v*t) for t in ss])
 def datum(ll):
  if y is not None:return mp.matrix([mp.mpf(str(t)) for t in y])
  return sum((mp.mpf(str(a))*atom(zz*mp.mpf(str(v))) for a,v in zip(w,o)),mp.zeros(21,1))-ll*A*atom(mp.mpf(10))
 def evaluate(u,ll):
  E=mp.matrix([[mp.cos(v*t)+mp.sin(v*t) for v in u] for t in ss]);yy=datum(ll)
  beta=mp.lu_solve(E.T*E,E.T*yy);res=yy-E*beta
  grad=mp.matrix([-2*beta[j]*sum(t*(mp.cos(u[j]*t)-mp.sin(u[j]*t))*res[k] for k,t in enumerate(ss))/A**2 for j in range(len(u))])
  cost=mp.fdot(res,res)/A**2
  return cost,grad,res,beta
 return A,datum,evaluate

def root_case(name,z,dps=70,start=None):
 r=CL[name][2];ls=LAM[r]
 previous=json.load(open(HERE.parent.parent/'reliability_followup/results/validation_raw.json'))['lambda_c'];guess=previous.get(f'{name}@{z}',ls)
 if start is None:
  h=mean(name,z,guess);reg=regression(h,r);ab=VP(h,r).wells(reg['nodes']);u0=[*ab[0]['u'],*ab[1]['u'],guess]
 else:u0=[mp.mpf(t) for t in start]
 A,datum,ev=direct(name,z,ls,r,dps=dps);K=r+1
 def equations(*x):
  ca,ga,_,_=ev(x[:K],x[-1]);cb,gb,_,_=ev(x[K:2*K],x[-1]);return [*ga,*gb,ca-cb]
 root=mp.findroot(equations,tuple(u0),tol=mp.mpf(10)**(-dps+18),maxsteps=100,solver='mdnewton')
 lam=root[-1];ua=list(root[:K]);ub=list(root[K:2*K]);ca,ga,ra,ba=ev(ua,lam);cb,gb,rb,bb=ev(ub,lam)
 # Noiseless linear-contrast root, computed by a separate direct polynomial regression.
 ss=[mp.mpf(j)/21 for j in range(-10,11)];X=mp.matrix([[t**n for n in range(2*r+1)]+[mp.cos(10*t)+mp.sin(10*t)] for t in ss])
 def contrast(ll):
  th=mp.lu_solve(X.T*X,X.T*datum(ll));m=[mp.factorial(j)*sig(j)*th[j] for j in range(2*r+1)]
  H=mp.matrix([[m[i+j] for j in range(r)] for i in range(r)]);bv=mp.matrix(m[r:2*r]);aa=m[-1]-(bv.T*mp.lu_solve(H,bv))[0]
  return -th[-1]-mp.mpf(str(ls))*aa
 ll=mp.findroot(contrast,(lam,lam+mp.mpf(str(ls))*mp.mpf('.001')),tol=mp.mpf(10)**(-dps+18))
 hess=[]
 for u in [ua,ub]:
  mat=mp.matrix([[mp.diff(lambda t:ev([t if k==j else u[k] for k in range(K)],lam)[1][i],u[j]) for j in range(K)] for i in range(K)])
  hess.append(mp.nstr(min(mp.eigsy(mat,eigvals_only=True)),24))
 def side(delta):
  la=lam+mp.mpf(str(ls))*delta
  aa=mp.findroot(lambda *u:tuple(ev(u,la)[1]),tuple(ua),tol=mp.mpf(10)**(-dps+18),maxsteps=25)
  bb=mp.findroot(lambda *u:tuple(ev(u,la)[1]),tuple(ub),tol=mp.mpf(10)**(-dps+18),maxsteps=25)
  return mp.nstr(ev(aa,la)[0]-ev(bb,la)[0],25)
 kappa=2
 return dict(cluster=name,z=z,dps=dps,lambda_c=mp.nstr(lam,dps-8),lambda_L=mp.nstr(ll,dps-8),offset_scaled=mp.nstr((lam-ll)/mp.mpf(str(z))**(2 if name!='asym3' else 1),30),uA=[mp.nstr(t,dps-8) for t in ua],uB=[mp.nstr(t,dps-8) for t in ub],gradient_max=mp.nstr(max([abs(t) for t in [*ga,*gb]]),20),equal_cost_error=mp.nstr(abs(ca-cb),20),min_hessian=hess,side_negative=side(mp.mpf('-1e-6')),side_positive=side(mp.mpf('1e-6')))

def main():
 cfg=json.load(open(HERE.parent/'config.json'));cases=[(n,z) for n in cfg['clusters'] for z in cfg['spacings']]+[('asym3',.3),('asym3',.15)]
 records=[]
 for name,z in cases:
  t=time.time();a=root_case(name,z,70);b=root_case(name,z,100,start=[*a['uA'],*a['uB'],a['lambda_c']])
  b['lambda_c_delta_70_100']=str(abs(mp.mpf(a['lambda_c'])-mp.mpf(b['lambda_c'])))
  records.append(b);json.dump(records,open(OUT/'roots.json','w'),indent=2)
  print(name,z,b['lambda_c'],b['offset_scaled'],'seconds',time.time()-t,flush=True)
if __name__=='__main__':main()
