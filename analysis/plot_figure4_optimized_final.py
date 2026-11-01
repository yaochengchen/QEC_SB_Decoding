"""Candidate only: archived accuracy joined to optimized per-syndrome timing."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"analysis_outputs/figure4_optimized_final"
OUT.mkdir(parents=True,exist_ok=True)
paths={"fast":"results_soe_fast/timing_v2_p0.10/timing_budget.csv", "baseline_time":"results/paper1_extra_timing_p0.10/timing_budget.csv", "accuracy":"results/paper1_extra/budget_ler.csv", "mps":"results/paper1_grid/grid_ler.csv", "bp":"results/paper1_extra/bposd_corr_ler.csv"}
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
protected=[p for p in (ROOT/"paper_draft_v1").rglob("*") if p.is_file() and p.suffix in (".md",".tex",".pdf",".png")]
before={str(p.relative_to(ROOT)):digest(p) for p in protected}
data={k:pd.read_csv(ROOT/v) for k,v in paths.items()}
def one(df,**kw):
    for k,v in kw.items(): df=df[np.isclose(df[k],v) if isinstance(v,float) else df[k].eq(v)]
    assert len(df)==1,(kw,len(df))
    return df.iloc[0]
rows=[]
for p in (.10,.12,.16):
    sub=data["accuracy"].query("d==13 and inference in ['r1','r2']")
    sub=sub[np.isclose(sub.p,p)]
    assert len(sub)==26
    for _,r in sub.iterrows():
        t=one(data["fast"],d=13,p=.1,config=r.config,stage="end_to_end_"+r.inference)
        rows.append(dict(p=p,config=r.config,method=r.solver,inference=r.inference,shots=int(r.shots),LER=r.LER,lo=r.LER_ci95_lo,hi=r.LER_ci95_hi,time_s=t.median_s,p90_s=t.p90_s,timing_n=int(t.syndromes),timing_p=.1,accuracy_source=paths["accuracy"],timing_source=paths["fast"]))
    for key,method,cfg in [("mps","MPS_chi16","mps16"),("bp","BP_OSD_corr","bposd3")]:
        r=one(data[key],d=13,p=p,method=method)
        t=one(data["baseline_time"],d=13,p=.1,config=cfg,stage="total")
        rows.append(dict(p=p,config=cfg,method=method,inference="baseline",shots=int(r.shots),LER=r.LER,lo=r.LER_ci95_lo,hi=r.LER_ci95_hi,time_s=t.median_s,p90_s=t.p90_s,timing_n=int(t.syndromes),timing_p=.1,accuracy_source=paths[key],timing_source=paths["baseline_time"]))
df=pd.DataFrame(rows)
assert len(df)==84 and (df.lo<=df.LER).all() and (df.LER<=df.hi).all()
# Independently verify the archived intervals are Wilson 95% intervals.
n=df.shots.to_numpy(); q=df.LER.to_numpy(); z=1.959963984540054
mid=(q+z*z/(2*n))/(1+z*z/n); half=z*np.sqrt(q*(1-q)/n+z*z/(4*n*n))/(1+z*z/n)
assert np.allclose(mid-half,df.lo,atol=1e-7) and np.allclose(mid+half,df.hi,atol=1e-7)
df.to_csv(OUT/"plot_data.csv",index=False)
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
colors={"SB":"#0072B2","SA":"#D55E00","MPS_chi16":"#222222","BP_OSD_corr":"#009E73"}
markers={"SB":"o","SA":"s","MPS_chi16":"*","BP_OSD_corr":"D"}
selected=["sbR8S100","sbR8S400","sbR32S400","saC8W100","saC32W200"]
def label(c):
    return {"sbR8S100":"SB 8×100","sbR8S400":"SB 8×400","sbR32S400":"SB 32×400","saC8W100":"SA 8×100","saC32W200":"SA 32×200","mps16":"MPS χ=16","bposd3":"BP+OSD (corr.)"}.get(c,c)
def point(ax,r,hollow=False):
    ax.errorbar(r.time_s,100*r.LER,yerr=[[100*(r.LER-r.lo)],[100*(r.hi-r.LER)]],fmt=markers[r.method],color=colors[r.method],mfc="white" if hollow else colors[r.method],ms=10 if r.config=="mps16" else 5,capsize=2,lw=.8,zorder=4 if r.inference=="baseline" else 3)
def setup(ax):
    ax.set_xscale("log"); ax.set_xlabel("Median end-to-end time per syndrome (s)")
    ax.set_ylabel("Logical error rate (%) ↓"); ax.grid(alpha=.18)
def overview(ax,p):
    sub=df[np.isclose(df.p,p)&(df.inference!="r1")]
    for _,r in sub.iterrows(): point(ax,r)
    for c in ["bposd3","sbR8S100","sbR8S400","sbR32S400","saC32W200","mps16"]:
        r=one(sub,config=c)
        offsets={"mps16":(-20,-13),"sbR8S400":(-48,18),"sbR32S400":(-5,45),"saC32W200":(-65,22),"sbR8S100":(6,8),"bposd3":(8,8)}
        ax.annotate(label(c),(r.time_s,100*r.LER),xytext=offsets[c],textcoords="offset points",fontsize=8,arrowprops=dict(arrowstyle="-",color=".65",lw=.5))
    setup(ax); ax.set_xlim(.0035,15)

fig,axes=plt.subplots(1,2,figsize=(12,4.5),gridspec_kw={"width_ratios":[1.3,1]})
overview(axes[0],.1);axes[0].set_title("(a) Accuracy–time tradeoff (radius two)")
ax=axes[1]
for i,c in enumerate(selected):
    a=one(df,p=.1,config=c,inference="r1");b=one(df,p=.1,config=c,inference="r2")
    ax.plot([100*a.LER,100*b.LER],[i,i],color=colors[a.method],lw=2)
    ax.plot(100*a.LER,i,marker=markers[a.method],mfc="white",mec=colors[a.method],ms=6)
    ax.plot(100*b.LER,i,marker=markers[b.method],color=colors[b.method],ms=6)
    ax.text(1.02,i,f"+{1e6*(b.time_s-a.time_s):.0f} μs",transform=ax.get_yaxis_transform(),va="center",fontsize=9)
ax.set_yticks(range(len(selected)),[label(c) for c in selected]);ax.invert_yaxis();ax.set_ylim(5.2,-.7)
ax.set_xlabel("Logical error rate (%) ←");ax.set_title("(b) Radius one → radius two")
ax.grid(axis="x",alpha=.18);ax.set_xlim(1.05,2.15)
ax.text(1.02,1.02,"Δ median time",transform=ax.transAxes,fontsize=8)
ax.text(.02,.02,"Open: radius one; filled: radius two\nPoint estimates; no paired inference in this panel",transform=ax.transAxes,fontsize=8)
handles=[Line2D([],[],color=colors[m],marker=markers[m],ls="",label=l) for m,l in [("SB","CP-SB"),("SA","CP-SA"),("MPS_chi16","MPS χ=16"),("BP_OSD_corr","Correlated BP+OSD")]]
fig.legend(handles=handles,loc="upper center",ncol=4,bbox_to_anchor=(.5,.94),frameon=False)
fig.suptitle("d = 13, p = 0.10 | Optimized CP implementation and archived baselines",fontsize=12)
fig.tight_layout(rect=(0,0,.97,.84));fig.subplots_adjust(wspace=.42)
for ext in ("png","pdf","svg"):fig.savefig(OUT/("figure4_optimized."+ext),dpi=240,bbox_inches="tight")
plt.close(fig)
print(OUT)
