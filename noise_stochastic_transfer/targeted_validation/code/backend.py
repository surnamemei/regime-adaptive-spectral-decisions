"""Scaled real variable projection and analytic moment-contrast derivatives."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1'); os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import solve_triangular
N=21; B=10.; S=np.arange(-10,11)/21
LAM={1:.032990205563495685,2:2.5936704891874683e-4}
ARCS={1:((-8.,-1.),(10.,15.)),2:((-15.,-5.),(12.,20.))}
CONST={1:dict(Phi1=1.114344120,sX=.249435265,unit=.0370563141),2:dict(Phi1=2.112466994e-3,sX=6.645863625e-4,unit=1.916883326e-4)}
CL={'pair':([1.,1.],[-1.,1.],1,2.),'sym3':([1.,1.,1.],[-1.,0.,1.],2,2/3),'asym3':([1.,2.,.5],[-1.,.3,1.5],2,1.9144115796098178)}
def sig(n):return 1 if n%4 in (0,1) else -1
def psi(u):
 a=np.outer(S,np.atleast_1d(u));return np.cos(a)+np.sin(a)
def dpsi(u):
 a=np.outer(S,np.atleast_1d(u));return S[:,None]*(np.cos(a)-np.sin(a))
def mean(name,z,lam):
 w,o,r,h=CL[name];A=h*z**(2*r)
 return psi(np.array(o)*z)@np.array(w)-lam*A*psi(B)[:,0]
def regression(h,r):
 import math
 X=np.column_stack([S**j for j in range(2*r+1)]+[psi(B)[:,0]])
 Q,R=np.linalg.qr(X);Ri=solve_triangular(R,np.eye(len(R)));th=solve_triangular(R,Q.T@h)
 m=np.array([math.factorial(j)*sig(j)*th[j] for j in range(2*r+1)])
 H=np.array([[m[i+j] for j in range(r)] for i in range(r)]);bv=m[r:2*r]
 q=np.linalg.solve(H,bv);A=float(m[2*r]-bv@q)
 da=np.zeros(2*r+1);da[-1]=1
 for i in range(r):
  da[r+i]-=2*q[i]
  for j in range(r):da[i+j]+=q[i]*q[j]
 g=np.r_[-LAM[r]*da*np.array([math.factorial(j)*sig(j) for j in range(2*r+1)]),-1.]
 se1=np.sqrt(g@(Ri@Ri.T)@g/2)
 f=float(-th[-1]-LAM[r]*A)
 coeff=np.r_[-q,1.];nodes=np.roots(coeff[::-1])
 if np.max(np.abs(nodes.imag))>1e-7:raise ValueError('nonreal moment initialization')
 return dict(A=A,m0=float(m[0]),f=f,se1=float(se1),nodes=np.sort(nodes.real),Hcond=float(np.linalg.cond(H)))
class VP:
 def __init__(self,h,r,sigma=0.):
  self.h=np.asarray(h,float);self.r=r
  Q,_=np.linalg.qr(np.column_stack([S**j for j in range(2*r)]))
  self.scale=max(np.linalg.norm(h-Q@(Q.T@h)),sigma,1e-14*np.linalg.norm(h))
 def values(self,u):
  E=psi(u);Q,R=np.linalg.qr(E);a=solve_triangular(R,Q.T@self.h)
  res=self.h-Q@(Q.T@self.h)
  return Q,R,a,res
 def residual(self,u):return self.values(u)[3]/self.scale
 def jac(self,u):
  Q,R,a,res=self.values(u);D=dpsi(u);out=np.empty((N,len(u)))
  for j in range(len(u)):
   t=D[:,j]*a[j];ek=np.zeros(len(u));ek[j]=D[:,j]@res
   out[:,j]=-(t-Q@(Q.T@t)+Q@solve_triangular(R.T,ek,lower=True))/self.scale
  return out
 def fit(self,u0,free=True,arc=None):
  u0=np.asarray(u0,float);idx=np.arange(len(u0) if free else len(u0)-1)
  def full(v):
   u=u0.copy();u[idx]=v;return u
  lo=np.full(len(idx),-np.inf);hi=np.full(len(idx),np.inf)
  if free and arc is not None:lo[-1],hi[-1]=arc
  x=np.clip(u0[idx],lo+1e-10,hi-1e-10)
  out=least_squares(lambda x:self.residual(full(x)),x,jac=lambda x:self.jac(full(x))[:,idx],bounds=(lo,hi),method='trf' if arc else 'lm',x_scale='jac',ftol=1e-11,xtol=1e-11,gtol=1e-11,max_nfev=400)
  u=full(out.x);res=self.values(u)[3]
  return dict(u=u,cost=float(res@res),res=res,status=int(out.status),nfev=int(out.nfev),optimality=float(out.optimality))
 def wells(self,nodes,step=.25):
  answer=[]
  for arc in ARCS[self.r]:
   loc=nodes.copy();profile=[]
   for v in np.arange(arc[0],arc[1]+1e-8,step):
    try:
     out=self.fit(np.r_[loc,v],free=False);loc=out['u'][:-1];profile.append(out)
    except (ValueError,np.linalg.LinAlgError,FloatingPointError):loc=nodes.copy()
   candidates=[]
   for out in sorted(profile,key=lambda x:x['cost'])[:3]:
    try:candidates.append(self.fit(out['u'],arc=arc))
    except (ValueError,np.linalg.LinAlgError,FloatingPointError):pass
   if not candidates:raise RuntimeError('all restricted starts failed')
   good=[a for a in candidates if a['status']>0]
   answer.append(min(good or candidates,key=lambda x:x['cost']))
  return answer
 def quick_wells(self,nodes,starts):
  ans=[]
  for arc,start in zip(ARCS[self.r],starts):
   ans.append(self.fit(np.r_[nodes,start],arc=arc))
  return ans

def statistics(h,sigma_hat,r,step=.25):
 reg=regression(h,r);vp=VP(h,r,sigma_hat);a,b=vp.wells(reg['nodes'],step)
 den=np.sqrt(2)*sigma_hat*np.linalg.norm(a['res']-b['res'])
 return dict(tL=reg['f']/(sigma_hat*reg['se1']),TC=(a['cost']-b['cost'])/den,A_hat=reg['A'],m0_hat=reg['m0'],sigma_hat=sigma_hat,Hcond=reg['Hcond'],status=min(a['status'],b['status']),nfev=a['nfev']+b['nfev'],optimality=max(a['optimality'],b['optimality']),uA=a['u'],uB=b['u'],JA=a['cost'],JB=b['cost'])
