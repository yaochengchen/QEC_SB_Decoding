from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
import json, time
import numpy as np
import pandas as pd

from .rotated_code import build_rotated_surface_code, sample_depolarizing
from .p5.statistics import exact_paired_statistics, wilson_interval as wilson
from .scalable import CSSMWPMDecoder, ScalableCPConfig, collect_cp_class_pools_no_oracle, cp_map_prediction, correction_class
from .adaptive_soe import fixed_soe_prediction
from .importance_mass import ImportanceMassConfig, adaptive_p41_prediction
from .bridge_ratio import BennettConfig, MultiStageBennettConfig, hybrid_bennett_refine


@dataclass
class P42AConfig:
    distance: int = 5
    p: float = 0.08
    shots: int = 500
    seed: int = 20261021
    solver_seed: int = 20261022
    cp: ScalableCPConfig | None = None
    importance: ImportanceMassConfig | None = None
    bennett: MultiStageBennettConfig | None = None

    def __post_init__(self):
        if self.cp is None:
            self.cp = ScalableCPConfig(
                replicas=8, steps=100, coupling_multipliers=(1.0,),
                random_restarts=1, warm_restarts=1,
            )
        if self.importance is None:
            self.importance = ImportanceMassConfig(
                stage0_top_m=1, stage0_radius=1, stage0_ratio_threshold=1.05,
                top_centers_per_class=2, core_radius=2,
                proposal_rhos=(0.08,0.20,0.40,0.50), uniform_mix=0.05,
                batch_size=32, min_samples=64, max_samples=160, z_score=2.58,
                prohibit_full_enumeration=True,
            )
        if self.bennett is None:
            self.bennett = MultiStageBennettConfig(
                lambdas=(0.0,0.25,0.5,0.75,1.0), n_samples=48, n_chains=4, burnin=48, thin=1,
                pair_flip_prob=0.20, bootstrap_reps=20, refine_ratio_threshold=1.20,
                decision_log_margin=0.03, sign_confidence=0.85,
            )


def _syn_key(sx, sz):
    out = 0
    for bit in np.concatenate([np.asarray(sx,dtype=np.uint8), np.asarray(sz,dtype=np.uint8)]):
        out = (out << 1) | int(bit)
    return out


def _paired_diff_ci(a, b):
    stats = exact_paired_statistics(a, b)
    return (stats['paired_difference'], stats['paired_ci95_lo'], stats['paired_ci95_hi'],
            stats['discordant_wins'], stats['discordant_losses'])


def run_p42a(outdir, cfg=P42AConfig()):
    outdir=Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    code=build_rotated_surface_code(cfg.distance)
    samples=sample_depolarizing(code,cfg.p,cfg.shots,cfg.seed)
    mwpm=CSSMWPMDecoder(code)
    models={}
    syndrome_rows=[]
    t_all=time.time()

    # Decode each distinct syndrome once; no exact enumeration or oracle is used.
    for i in range(cfg.shots):
        sx=samples['sx'][i]; sz=samples['sz'][i]; key=_syn_key(sx,sz)
        if key in models: continue
        t0=time.time(); ex0,ez0=mwpm.decode(sx,sz); t_mw=time.time()-t0
        mw_cls=correction_class(code,ex0,ez0)
        t0=time.time(); pools,hams,cpstats=collect_cp_class_pools_no_oracle(
            code,ex0,ez0,cfg.solver_seed+key*1009,cfg.cp); t_cp=time.time()-t0
        cp_cls,cp_w=cp_map_prediction(pools,hams)
        t0=time.time(); local_cls,local_scores,local_counts,_=fixed_soe_prediction(pools,hams,cfg.p,1,1); t_local=time.time()-t0
        t0=time.time(); p41_cls,p41_diag=adaptive_p41_prediction(
            pools,hams,cfg.p,rng_seed=cfg.solver_seed+7000003+key*7919,cfg=cfg.importance); t_p41=time.time()-t0
        t0=time.time(); bar_cls,bar_diag=hybrid_bennett_refine(
            p41_cls,local_scores,pools,hams,cfg.p,seed=cfg.solver_seed+9000007+key*6151,cfg=cfg.bennett); t_bar=time.time()-t0
        models[key]={
            'mwpm':int(mw_cls),'cp_map':int(cp_cls),'local':int(local_cls),
            'p41':int(p41_cls),'bar':int(bar_cls),
            'p41_diag':p41_diag,'bar_diag':bar_diag,
        }
        syndrome_rows.append({
            'syndrome_key':int(key),
            'mwpm_class':int(mw_cls),'cp_map_class':int(cp_cls),'local_soe_class':int(local_cls),
            'p41_importance_class':int(p41_cls),'p42_bennett_class':int(bar_cls),
            'cp_best_weight_0':cpstats[0]['best_weight'],'cp_best_weight_1':cpstats[1]['best_weight'],
            'cp_best_weight_2':cpstats[2]['best_weight'],'cp_best_weight_3':cpstats[3]['best_weight'],
            'cp_pool_total':int(sum(cpstats[c]['unique_cp_candidates'] for c in range(4))),
            'local_candidates':int(local_counts.sum()),
            'p41_stage':int(p41_diag['stage']),
            'p41_samples':int(p41_diag['samples_drawn']),
            'p41_candidate_evaluations':int(p41_diag['candidate_evaluations']),
            'bar_refined':int(bar_diag.get('refined',0)),
            'bar_confirmed':int(bar_diag.get('confirmed',0)),
            'bar_confirmation_agrees':int(bar_diag.get('confirmation_agrees',0)),
            'bar_local_top2_ratio':float(bar_diag.get('local_top2_ratio',np.nan)),
            'bar_log_ratio':float(bar_diag.get('log_ratio_b_over_a',np.nan)),
            'bar_sign_confidence':float(bar_diag.get('sign_confidence',np.nan)),
            'bar_acceptance_mean':float(bar_diag.get('mean_acceptance',np.nan)),
            'bar_samples_per_state':int(bar_diag.get('samples_per_bridge_state',0)),
            'bar_fallback':bar_diag.get('fallback',''),
            'time_mwpm_s':t_mw,'time_cp_s':t_cp,'time_local_s':t_local,'time_p41_s':t_p41,'time_bar_s':t_bar,
        })

    methods=['MWPM','CP_SB_MAP','CP_local_SOE','P4_1_importance','P4_2_Bennett']
    correct={m:np.zeros(cfg.shots,dtype=np.uint8) for m in methods}
    shot_rows=[]
    for i in range(cfg.shots):
        key=_syn_key(samples['sx'][i],samples['sz'][i]); mdl=models[key]
        truth=int(samples['logical'][i,0])*2+int(samples['logical'][i,1])
        pred={'MWPM':mdl['mwpm'],'CP_SB_MAP':mdl['cp_map'],'CP_local_SOE':mdl['local'],
              'P4_1_importance':mdl['p41'],'P4_2_Bennett':mdl['bar']}
        row={'shot':i,'syndrome_key':int(key),'truth_class':truth}
        for m in methods:
            correct[m][i]=int(pred[m]==truth); row[f'pred_{m}']=int(pred[m]); row[f'correct_{m}']=int(correct[m][i])
        shot_rows.append(row)

    summary=[]
    for m in methods:
        k=int(correct[m].sum()); lo,hi=wilson(k,cfg.shots)
        summary.append({'distance':cfg.distance,'p':cfg.p,'shots':cfg.shots,'method':m,
                        'logical_success':k/cfg.shots,'logical_error':1-k/cfg.shots,'ci95_lo':lo,'ci95_hi':hi})
    sdf=pd.DataFrame(summary)
    comps=[]
    for a,b in [('P4_2_Bennett','MWPM'),('P4_2_Bennett','CP_SB_MAP'),('P4_2_Bennett','CP_local_SOE'),('P4_2_Bennett','P4_1_importance'),('P4_1_importance','MWPM')]:
        d,lo,hi,ao,bo=_paired_diff_ci(correct[a],correct[b]); comps.append({'method_a':a,'method_b':b,'accuracy_diff_a_minus_b':d,'paired_ci95_lo':lo,'paired_ci95_hi':hi,'a_correct_b_wrong':ao,'a_wrong_b_correct':bo})
    syn_df=pd.DataFrame(syndrome_rows)
    pd.DataFrame(shot_rows).to_csv(outdir/'shot_level.csv',index=False)
    syn_df.to_csv(outdir/'syndrome_level.csv',index=False); sdf.to_csv(outdir/'summary.csv',index=False); pd.DataFrame(comps).to_csv(outdir/'paired_comparisons.csv',index=False)
    meta={
        'config':asdict(cfg),'elapsed_s':time.time()-t_all,'unique_syndromes':len(models),
        'exact_enumeration_used':False,'oracle_used':False,
        'mean_unique_decode_time_s':float((syn_df[['time_mwpm_s','time_cp_s','time_local_s','time_p41_s','time_bar_s']].sum(axis=1)).mean()),
        'fraction_bar_refined':float(syn_df['bar_refined'].mean()),
        'mean_bar_acceptance':float(syn_df.loc[syn_df['bar_refined']==1,'bar_acceptance_mean'].mean()) if np.any(syn_df['bar_refined']==1) else float('nan'),
        'mean_p41_candidate_evaluations':float(syn_df['p41_candidate_evaluations'].mean()),
    }
    (outdir/'metadata.json').write_text(json.dumps(meta,indent=2,default=str),encoding='utf-8')
    return sdf,pd.DataFrame(comps),syn_df,meta
