"""Replot archived budget evidence; never runs decoding or overwrites paper figures."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis_outputs' / 'evidence_budget'
OUT.mkdir(parents=True, exist_ok=True)
sources = [ROOT/'results/paper1_extra/budget_ler.csv',
           ROOT/'results_soe_fast/timing_v2_p0.10/work_counts.csv',
           ROOT/'results/paper1_extra/budget_paired_matched.csv']
B, W, P = [pd.read_csv(p) for p in sources]
assert not B.duplicated(['d','p','config','inference']).any()
assert not W.duplicated(['d','config']).any()
D = B[B.inference == 'r2'].merge(W, on=['d','config'], validate='many_to_one', indicator=True)
assert D['_merge'].eq('both').all()
assert D.shots.eq(10000).all()
assert np.allclose(W.dE_evaluations, [8*int(c.split('S' if c.startswith('sb') else 'W')[1])*int(c.split('R' if c.startswith('sb') else 'C')[1].split('S' if c.startswith('sb') else 'W')[0])*m for c,m in zip(W.config,W.m)])
plt.rcParams.update({'font.size':9,'axes.titlesize':10,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none'})
groups = [(['sbR8S100','sbR8S200','sbR8S400'],'SB: 8 replicas, 100/200/400 steps','#2474bc','o','-'),
          (['sbR8S100','sbR16S100','sbR32S100','sbR64S100'],'SB: 100 steps, 8/16/32/64 replicas','#2474bc','s',':'),
          (['sbR32S400'],'SB: 32 replicas, 400 steps','#174478','*','none'),
          (['saC8W25','saC8W50','saC8W100','saC8W200'],'SA: 8 chains, 25/50/100/200 sweeps','#d86927','^','-'),
          (['saC16W100','saC32W200'],'SA: 16x100 / 32x200','#9d4716','P','none')]
def save(fig,name):
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/f'{name}.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)

for p in [.08,.10,.12,.14,.16,.18]:
    fig,axs=plt.subplots(2,3,figsize=(12,7),sharey=True)
    for col,d in enumerate([9,11,13]):
        cell=D[(D.d==d)&np.isclose(D.p,p)].set_index('config')
        for row,resource in enumerate(['dE_evaluations','serial_rounds']):
            ax=axs[row,col]
            for cfgs,label,color,mk,ls in groups:
                q=cell.reindex([c for c in cfgs if c in cell.index])
                if len(q):
                    # Replica changes share a serial-round coordinate; do not join them.
                    style='none' if row==1 and 'replicas' in label and '100 steps,' in label else ls
                    ax.plot(q[resource],q.LER_over_REF,color=color,marker=mk,ls=style,label=label,ms=5)
            ax.axhline(1,color='black',ls=':',lw=1)
            ax.set_xscale('log');ax.grid(alpha=.2)
            ax.set_title(f'd = {d}')
            ax.set_xlabel('Single-flip energy evaluations / decode' if row==0 else 'Sequential update rounds / decode')
            if col==0: ax.set_ylabel('Logical error rate / reference LER')
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1.00),ncol=3)
    fig.suptitle(f'Budget study: p = {p:.2f}, SOE radius 2, first 10,000 paired shots per cell',y=1.055)
    fig.subplots_adjust(top=.85,hspace=.35,wspace=.15)
    save(fig,f'resources_p{p:.2f}')

comparisons=[('sbR8S400_r2','sbR32S100_r2','More steps vs more replicas\nSB 8x400 minus SB 32x100'),
             ('sbR8S400_r2','saC8W100_r2','Moderate budgets\nSB 8x400 minus SA 8x100'),
             ('sbR32S400_r2','saC32W200_r2','Largest budgets\nSB 32x400 minus SA 32x200')]
fig,axs=plt.subplots(1,3,figsize=(13,4.5),sharey=True)
used=[]; report=[]
for ax,(a,b,title) in zip(axs,comparisons):
    q=P[(P.a==a)&(P.b==b)].copy()
    assert not q.duplicated(['d','p']).any()
    assert np.allclose(q.LER_a_minus_LER_b,(q.a_wrong_b_right-q.a_right_b_wrong)/q.shots)
    for d,color,mk,shift in [(9,'#2474bc','o',-.0015),(11,'#d86927','s',0),(13,'#68469e','^',.0015)]:
        z=q[q.d==d].sort_values('p')
        y=100*z.LER_a_minus_LER_b
        ax.errorbar(z.p+shift,y,yerr=[y-100*z.ci95_lo,100*z.ci95_hi-y],fmt='none',ecolor=color,capsize=3,lw=1)
        for _,r in z.iterrows():
            ax.plot(r.p+shift,100*r.LER_a_minus_LER_b,marker=mk,color=color,mfc=color if r.holm_p<.05 else 'white',ls='none',ms=6)
        ax.plot([],[],marker=mk,color=color,ls='none',label=f'd={d}')
    ax.axhline(0,color='black',ls='--',lw=1);ax.grid(axis='y',alpha=.2)
    ax.set_title(title);ax.set_xlabel('Physical error rate p (markers offset for readability)')
    ax.set_xticks(sorted(q.p.unique()));ax.legend()
    used.append(q)
    report.append({'a':a,'b':b,'cells':len(q),'a_lower_LER_cells':int((q.LER_a_minus_LER_b<0).sum()),'holm_significant_cells':int((q.holm_p<.05).sum()),'wins':int(q.a_right_b_wrong.sum()),'losses':int(q.a_wrong_b_right.sum())})
axs[0].set_ylabel('LER(A) - LER(B), percentage points\nNegative values favor A')
fig.suptitle('SOE radius 2: paired differences with archived conservative 95% intervals',y=1.04)
fig.text(.5,-.045,'Filled markers: archived Holm p < 0.05. Open markers: not significant. Intervals are not Holm-adjusted.',ha='center',fontsize=9)
fig.tight_layout()
save(fig,'paired_improvements_r2')
D.drop(columns='_merge').to_csv(OUT/'resource_plot_data.csv',index=False)
pd.concat(used).to_csv(OUT/'paired_plot_data.csv',index=False)
(OUT/'manifest.json').write_text(json.dumps({'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'inference':'r2','shots_per_cell':10000,'paired_summary':report,'resource_y':'LER_over_REF (point estimates; no ratio CI inferred)','paired_interval':'archived conservative Bonferroni Clopper-Pearson difference; see qec_sb/p5/statistics.py','holm':'archived families unchanged'},indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
print('Saved figures and provenance to',OUT)
