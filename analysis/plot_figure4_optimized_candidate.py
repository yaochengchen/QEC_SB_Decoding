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
OUT=ROOT/"analysis_outputs/figure4_optimized_candidate"
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
fig,axes=plt.subplots(1,2,figsize=(12.4,4.8))
overview(axes[0],.1); axes[0].set_title("(a) All search budgets, radius-two aggregation")
ax=axes[1]
for c in selected:
    sub=df[np.isclose(df.p,.1)&df.config.eq(c)]
    a=one(sub,inference="r1"); b=one(sub,inference="r2")
    ax.plot([a.time_s,b.time_s],100*np.array([a.LER,b.LER]),color=colors[a.method],lw=1.5)
    point(ax,a,True); point(ax,b)
    ax.annotate(label(c),(b.time_s,100*b.LER),xytext=({"sbR8S100":(5,12),"sbR8S400":(5,12),"sbR32S400":(-68,45),"saC8W100":(8,20),"saC32W200":(-62,16)}[c]),textcoords="offset points",fontsize=8)
point(ax,one(df,p=.1,config="mps16"))
setup(ax); ax.set_xlim(.085,14); ax.set_title("(b) Selected budgets: radius one → radius two")
handles=[Line2D([],[],color=colors[m],marker=markers[m],ls="",label=l) for m,l in [("SB","CP-SB"),("SA","CP-SA"),("MPS_chi16","MPS χ=16"),("BP_OSD_corr","Correlated BP+OSD")]]
fig.legend(handles=handles,loc="upper center",ncol=4,bbox_to_anchor=(.5,.94),frameon=False)
ax.text(.98,.98,"Open: radius one; filled: radius two\nVertical bars: pointwise Wilson 95% CI",transform=ax.transAxes,ha="right",va="top",fontsize=8)
fig.suptitle("Candidate Fig. 4 | d = 13, accuracy and timing at p = 0.10",y=1.01,fontsize=13)
fig.text(.5,.01,"Accuracy: 10,000 shots per CP configuration; 20,000 per baseline. Timing: 100 unique nontrivial syndromes; preprocessing excluded.",ha="center",fontsize=8)
fig.tight_layout(rect=(0,.04,1,.86))
for ext in ("png","pdf","svg"): fig.savefig(OUT/("figure4_candidate."+ext),dpi=220,bbox_inches="tight")
plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12.4,4.8))
for ax,p in zip(axes,(.12,.16)):
    overview(ax,p); ax.set_title(f"Accuracy p = {p:.2f}; timing p = 0.10")
fig.suptitle("Supporting sensitivity views | d = 13 | different-noise timing proxy",fontsize=13)
fig.legend(handles=handles,loc="upper center",ncol=4,bbox_to_anchor=(.5,.94),frameon=False)
fig.text(.5,.01,"Not matched-noise latency measurements. Bars: Wilson 95% CI. CP: N=10,000; baselines: N=20,000 (p=.12), 10,000 (p=.16).",ha="center",fontsize=8)
fig.tight_layout(rect=(0,.04,1,.86))
for ext in ("png","pdf","svg"): fig.savefig(OUT/("figure4_noise_slices."+ext),dpi=220,bbox_inches="tight")
plt.close(fig)
assert all(digest(ROOT/p)==h for p,h in before.items()),"Protected manuscript artifact changed"
manifest={"inputs":{v:digest(ROOT/v) for v in paths.values()},"protected_files_checked":len(before),"protected_unchanged":True,"points":len(df),"wilson_intervals_verified":True,"timing_definition":"median of per-syndrome seed+search+SOE sums; 100 unique nontrivial syndromes; per-code preprocessing excluded", "limitations":["Archived accuracy linked to optimized implementation timing; no new decoding experiment", "Baseline implementation timings retained from original campaign", "No timing confidence intervals; p90 is a percentile, not a CI", "No paired significance or equivalence test against MPS"]}
(OUT/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
print(df[np.isclose(df.p,.1)&df.config.isin(selected+["mps16","bposd3"])][["config","inference","LER","time_s","p90_s"]].to_string(index=False))
print("Verified",len(before),"protected artifacts unchanged")
