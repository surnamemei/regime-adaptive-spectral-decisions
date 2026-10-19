"""Figures and compact result tables from the fixed analyses."""
from pathlib import Path
import csv,json,math,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
def read(name):return list(csv.DictReader(open(ROOT/'results'/name)))
def fval(x):return float(x)
plt.rcParams.update({'font.family':'serif','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'pdf.fonttype':42})
COL={'S':'#286ec0','I':'#d65f23','C':'#28927b','L':'#8b58a6'}
D=read('old_curves.csv');M=read('old_matched.csv');DF=read('fresh_curves.csv');MF=read('fresh_matched.csv')
fig,axs=plt.subplots(2,3,figsize=(6.8,4.6),sharex=True,sharey=True)
for row,(data,match,prefix) in enumerate([(D,M,'Original'),(DF,MF,'Fresh')]):
 for j,(g,title) in enumerate([('nu=0.005',r'$\nu=0.005$'),('nu=0.5',r'$\nu=0.5$'),('nu=1.0',r'$\nu=1$')]):
  ax=axs[row,j]
  for rule in ['S','I','C','L']:
   rr=[r for r in data if r['group']==g and r['rule']==rule]
   x=[fval(r['wrong']) for r in rr];y=[fval(r['correct']) for r in rr];ax.plot(x,y,'o-',color=COL[rule],ms=2.3,lw=1,label=rule)
   i=next(i for i,r in enumerate(rr) if fval(r['alpha'])==.05);ax.plot(x[i],y[i],marker='s',mfc='white',mec=COL[rule],ms=4.5)
  ax.set_title(prefix+': '+title,fontsize=8.5);ax.grid(alpha=.2);ax.set_xlim(0,.08);ax.set_ylim(0,.8)
  for rule in ['S','I','C']:
   r=next(r for r in match if r['group']==g and fval(r['target_wrong'])==.02)
   value = fval(r[rule+'_correct']) if rule != 'C' else fval(next(t for t in read(('old' if row==0 else 'fresh')+'_matched_cost.csv') if t['group']==g and fval(t['target_wrong'])==.02)['C_correct'])
   if math.isfinite(value):ax.plot(.02,value,'x',color=COL[rule],ms=6)
for ax in axs[-1,:]:ax.set_xlabel('wrong probability')
for ax in axs[:,0]:ax.set_ylabel('correct probability')
axs[0,0].legend(ncol=4,loc='lower right',fontsize=6.5);fig.tight_layout();fig.savefig(OUT/'fig_tradeoff.pdf');plt.close(fig)

roots=json.load(open(ROOT/'results/roots.json'));F=read('fresh_cells.csv')
fig,axs=plt.subplots(2,2,figsize=(6.8,4.5),sharex='col')
for j,name in enumerate(['pair','sym3']):
 rr=sorted([r for r in roots if r['cluster']==name],key=lambda r:r['z']);axs[0,j].plot([r['z'] for r in rr],[float(r['offset_scaled']) for r in rr],'o-',color=COL['S'])
 axs[0,j].set_title('Equal pair' if name=='pair' else 'Symmetric triple');axs[0,j].set_ylabel(r'$(\lambda_c-\lambda_L)/z^2$');axs[0,j].ticklabel_format(axis='y',style='sci',scilimits=(0,0));axs[0,j].grid(alpha=.2)
 for nu,col in zip([.005,.5,1.],['#d65f23','#28927b','#286ec0']):
  xx=[];yy=[];err=[]
  for z in [.075,.15,.3]:
   rr=[r for r in F if r['cluster']==name and fval(r['z'])==z and fval(r['nu'])==nu and r['rule']=='S' and abs(fval(r['position']))==.1]
   pp=[int(r['wrong'])/int(r['n']) for r in rr];var=sum(p*(1-p)/(int(r['n'])-1) for p,r in zip(pp,rr))/4
   xx.append(z);yy.append(np.mean(pp));err.append(1.959964*np.sqrt(var))
  axs[1,j].errorbar(xx,yy,yerr=err,marker='o',ms=3,color=col,lw=1,label=r'$\nu=%g$'%nu,capsize=2)
 axs[1,j].axhline(.05,color='black',lw=.7);axs[1,j].set_ylim(0,.085);axs[1,j].set_xlabel('spacing z');axs[1,j].set_ylabel('near-tie wrong rate');axs[1,j].grid(alpha=.2)
axs[1,0].legend(ncol=3,fontsize=7,loc='lower right');fig.tight_layout();fig.savefig(OUT/'fig_convergence.pdf');plt.close(fig)

# Main matched-error table. Each row is averaged over its own fixed design.
lines=[r'\begin{tabular}{@{}llrrrr@{}}',r'\toprule',r' study & noise & S correct & C correct & I correct & S $-$ C [95\%]\\',r'\midrule']
for tag,match in [('Original',M),('Fresh',MF)]:
 cost=read(('old' if tag=='Original' else 'fresh')+'_matched_cost.csv')
 for g,label in [('all','All'),('nu=0.005',r'$0.005$'),('nu=0.5',r'$0.5$'),('nu=1.0',r'$1$')]:
  r=next(r for r in match if r['group']==g and fval(r['target_wrong'])==.02)
  c=next(r for r in cost if r['group']==g and fval(r['target_wrong'])==.02)
  lines.append(tag+' & '+label+' & '+f"{100*fval(r['S_correct']):.2f}\\% & {100*fval(c['C_correct']):.2f}\\% & {100*fval(r['I_correct']):.2f}\\% & {100*fval(c['difference']):+.2f} [{100*fval(c['low']):+.2f},{100*fval(c['high']):+.2f}]\\\\")
lines += [r'\bottomrule',r'\end{tabular}'];(OUT/'table_matched.tex').write_text('\n'.join(lines)+'\n')
# Compact sensitivity, all 36,000 fresh draws.
T=read('sensitivity.csv');lines=[r'\begin{tabular}{@{}rrrrrr@{}}',r'\toprule',r'$r\gamma$ & multiplier & correct & wrong & abstain & max. cell wrong\\',r'\midrule']
for frac in [.125,.25,.375]:
 for mul in [.5,1.,2.]:
  rr=[r for r in T if fval(r['gamma_fraction'])==frac and fval(r['multiplier'])==mul]
  co=sum(int(r['correct']) for r in rr);wr=sum(int(r['wrong']) for r in rr);n=sum(int(r['n']) for r in rr);ma=max(fval(r['wrong_rate']) for r in rr)
  lines.append(f'{frac:g} & {mul:g} & {co} & {wr} & {n-co-wr} & {100*ma:.1f}\\%\\\\')
lines += [r'\bottomrule',r'\end{tabular}'];(OUT/'table_sensitivity.tex').write_text('\n'.join(lines)+'\n')
# One row per new cell; retain all methods in the supplement.
lines=[r'\begin{longtable}{@{}lrrrlrrrr@{}}',r'\caption{Fresh-cell correct/wrong/abstention counts, 500 draws per cell. The last column is the selection rule\textquotesingle s wrong count divided by its declaration count.}\label{tab:newcells}\\',r'\toprule',r'case & $z$ & $\nu$ & position & truth & S & I & C & S wrong/decl.\\',r'\midrule\endfirsthead',r'\toprule',r'case & $z$ & $\nu$ & position & truth & S & I & C & S wrong/decl.\\',r'\midrule\endhead']
meta=json.load(open(ROOT/'results/fresh_meta.json'));lookup={(r['cluster'],fval(r['z']),fval(r['nu']),fval(r['position']),r['rule']):r for r in F}
for c in meta['cells']:
 key=(c['cluster'],c['z'],c['nu'],c['position']);rr=[lookup[(*key,r)] for r in ['S','I','C']];rs=rr[0]
 cells=[f"{r['correct']}/{r['wrong']}/{r['abstain']}" for r in rr];wrong=int(rs['wrong']);decl=int(rs['correct'])+wrong
 lines.append(('pair' if c['cluster']=='pair' else 'triple')+f" & {c['z']:g} & {c['nu']:g} & {c['position']:+g} & {c['truth']} & "+' & '.join(cells)+f' & {wrong}/{decl}\\\\')
lines += [r'\bottomrule',r'\end{longtable}'];(OUT/'table_fresh_cells.tex').write_text('\n'.join(lines)+'\n')
lines=[r'\begin{tabular}{@{}lrrrr@{}}',r'\toprule',r'case & $z$ & $\lambda_c$ & $(\lambda_c-\lambda_L)/z^\kappa$ & min. Hessian\\',r'\midrule']
for a in roots:lines.append(a['cluster']+f" & {a['z']:g} & {float(a['lambda_c']):.12g} & {float(a['offset_scaled']):+.5g} & {min(float(v) for v in a['min_hessian']):.3g}\\\\")
lines += [r'\bottomrule',r'\end{tabular}'];(OUT/'table_roots.tex').write_text('\n'.join(lines)+'\n')
# Full curve tables and data live in CSV; the manuscript takes the compact figures.
print('Figures and TeX tables generated.')
