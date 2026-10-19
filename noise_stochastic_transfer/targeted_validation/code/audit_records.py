"""Preselected old records, refitted independently at 70/100 digits."""
import sys,json,time,math,csv
from pathlib import Path
import numpy as np
from scipy.stats import t as tdist
from backend import *
from high_precision import direct
import mpmath as mp
ROOT=Path(__file__).resolve().parent.parent
D=json.load(open(ROOT.parent/'adaptive_selection/results/raw.json'));CFG=json.load(open(ROOT/'config.json'))
roots={(a['cluster'],a['z']):a for a in json.load(open(ROOT/'results/roots.json'))}
def mp_linear(h,r,sigma):
 ss=[mp.mpf(j)/21 for j in range(-10,11)];X=mp.matrix([[t**n for n in range(2*r+1)]+[mp.cos(10*t)+mp.sin(10*t)] for t in ss]);yy=mp.matrix([mp.mpf(str(t)) for t in h]);inv=(X.T*X)**-1;th=inv*X.T*yy
 mm=[mp.factorial(j)*sig(j)*th[j] for j in range(2*r+1)];H=mp.matrix([[mm[i+j] for j in range(r)] for i in range(r)]);bv=mp.matrix(mm[r:2*r]);q=mp.lu_solve(H,bv);aa=mm[-1]-(bv.T*q)[0]
 da=[mp.mpf(0)]*(2*r+1);da[-1]=1
 for i in range(r):
  da[r+i]-=2*q[i]
  for j in range(r):da[i+j]+=q[i]*q[j]
 gg=mp.matrix([-mp.mpf(str(LAM[r]))*da[j]*mp.factorial(j)*sig(j) for j in range(2*r+1)]+[-1]);se=sigma*mp.sqrt((gg.T*inv*gg)[0]/2)
 return (-th[-1]-mp.mpf(str(LAM[r]))*aa)/se,aa,mm[0]
def solve_h(name,z,h,sigma,initial,dps):
 r=CL[name][2];A,_,ev=direct(name,z,0,r,y=h,dps=dps);ans=[]
 for u in initial:
  rr=mp.findroot(lambda *x:tuple(ev(x,0)[1]),tuple(mp.mpf(str(v)) for v in u),tol=mp.mpf(10)**(-dps+20),maxsteps=70)
  ans.append((list(rr),ev(list(rr),0)))
 ca,ga,ra,_=ans[0][1];cb,gb,rb,_=ans[1][1]
 tc=(ca-cb)*A*A/(mp.sqrt(2)*mp.mpf(str(sigma))*mp.norm(ra-rb))
 tl,aa,m0=mp_linear(h,r,mp.mpf(str(sigma)))
 nuhat=mp.mpf(str(sigma))/(mp.mpf(str(CONST[r]['unit']))*aa) if aa>0 else mp.inf
 sel=bool(aa>0 and m0>0 and nuhat<=(aa/m0)**(mp.mpf(1)/(4*r)))
 return dict(tL=mp.nstr(tl,50),TC=mp.nstr(tc,50),sel=sel,gradient=mp.nstr(max(abs(t) for t in [*ga,*gb]),15),nodes=[[mp.nstr(v,60) for v in u[0]] for u in ans])
def dec(t):return 'A' if float(t)<-float(tdist.ppf(.95,21)) else ('B' if float(t)>float(tdist.ppf(.95,21)) else 'inc')
rows=[];labels=[]
for ci,c in enumerate(D['configs']):
 if [c['cluster'],c['z']] not in CFG['hard_old_cases']:continue
 if not ((c['nu']==.005 and c['pos'] in ['m0.1','p0.1']) or (c['nu']==.5 and c['pos']=='between')):continue
 name,z=c['cluster'],c['z'];w,o,r,hr=CL[name];lc=mp.mpf(roots[name,z]['lambda_c']);gap=float(c['lam'])-float(lc)
 # Independently check noiseless cost signs at each audited cell.
 init=[roots[name,z]['uA'],roots[name,z]['uB']];A,_,ev=direct(name,z,c['lam'],r,dps=100)
 cc=[]
 for u in init:
  u=mp.findroot(lambda *x:tuple(ev(x,mp.mpf(str(c['lam'])))[1]),tuple(mp.mpf(v) for v in u),tol=mp.mpf('1e-75'),maxsteps=70)
  cc.append(ev(u,mp.mpf(str(c['lam'])))[0])
 labels.append(dict(cluster=name,z=z,nu=c['nu'],position=c['pos'],stored_truth=c['truth'],mp_truth='A' if cc[0]<cc[1] else 'B',normalized_cost_difference=mp.nstr(cc[0]-cc[1],30),lambda_distance=gap))
 rng=np.random.default_rng(np.random.SeedSequence([2026101902,ci,0]));xx=(np.array(w)[None,:]*np.exp(1j*z*np.outer(S,o))).sum(1)-c['lam']*c['A']*np.exp(1j*10*S)
 for di in [0,1]:
  noise=(rng.normal(size=21)+1j*rng.normal(size=21))/np.sqrt(2);yy=xx+c['sigma']*noise;yh=(yy+np.conj(yy[::-1]))/2;h=yh.real+yh.imag;sh=c['records'][di]['sigma_hat']
  start=time.time();rg=regression(h,r);ab=VP(h,r,sh).wells(rg['nodes'],step=.125)
  a=solve_h(name,z,h,sh,[v['u'] for v in ab],70);b=solve_h(name,z,h,sh,a['nodes'],100);old=c['records'][di]
  out=dict(cluster=name,z=z,nu=c['nu'],position=c['pos'],draw=di,old_tL=old['tL'],old_TC=old['TC'],mp_tL=float(b['tL']),mp_TC=float(b['TC']),delta_tL=float(b['tL'])-old['tL'],delta_TC=float(b['TC'])-old['TC'],agreement_70_100=max(abs(float(mp.mpf(a[k])-mp.mpf(b[k]))) for k in ['tL','TC']),gradient=b['gradient'],old_dS=old['dS'],mp_dS=dec(b['TC'] if b['sel'] else b['tL']),old_dI=old['dI'],mp_dI=dec(b['tL']) if dec(b['tL'])==dec(b['TC']) else 'inc')
  rows.append(out);json.dump(dict(labels=labels,records=rows),open(ROOT/'results/hard_audit.json','w'),indent=2)
  print(name,z,c['nu'],c['pos'],di,'dTC',out['delta_TC'],'decisions',out['old_dS'],out['mp_dS'],'seconds',time.time()-start,flush=True)
print('complete',len(labels),len(rows),flush=True)
