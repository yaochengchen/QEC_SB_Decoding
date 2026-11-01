from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit



def _weights(ham, u):
    phys = ham.physical_from_bits(np.asarray(u, dtype=np.uint8))
    n = ham.n_data
    return np.count_nonzero(phys[:, :n] | phys[:, n:], axis=1).astype(np.int16)


class PosteriorChainSampler:
    """Metropolis sampler of one logical class in stabilizer coordinates.

    Target density is proportional to r**PauliWeight(u), where
    r=(p/3)/(1-p). Local 1/2-generator flips are supplemented by symmetric
    XOR jumps learned from diverse CP-SB candidates. Because proposals are
    symmetric, the MH acceptance ratio depends only on the Pauli-weight change.
    """
    def __init__(self, ham, p, cp_pool, seed, n_chains=8, burnin=96, thin=2,
                 pair_flip_prob=0.20, cp_jump_prob=0.15):
        self.ham = ham
        self.p = float(p)
        self.ratio = (self.p / 3.0) / (1.0 - self.p)
        self.logr = float(np.log(self.ratio))
        self.rng = np.random.default_rng(int(seed))
        self.n_chains = int(n_chains)
        self.thin = int(thin)
        self.pair_flip_prob = float(pair_flip_prob)
        self.cp_jump_prob = float(cp_jump_prob)
        pool = np.asarray(cp_pool, dtype=np.uint8)
        if pool.ndim != 2 or pool.shape[1] != ham.n_spins:
            raise ValueError("cp_pool shape mismatch")
        # Start chains from the best/diverse CP candidates, cyclically.
        w = _weights(ham, pool)
        order = np.argsort(w, kind="stable")
        pool = pool[order]
        ids = np.arange(self.n_chains) % len(pool)
        self.u = pool[ids].copy()
        self.phys = ham.physical_from_bits(self.u)
        n = ham.n_data
        self.w = np.count_nonzero(self.phys[:, :n] | self.phys[:, n:], axis=1).astype(np.int16)
        self.G = ham.generators.astype(np.uint8)
        # Symmetric CP-derived XOR jump masks.
        base = pool[0]
        masks = np.unique(pool ^ base[None, :], axis=0)
        masks = masks[np.any(masks, axis=1)]
        self.jump_masks = masks
        if len(masks):
            self.jump_phys = ((masks @ self.G) & 1).astype(np.uint8)
        else:
            self.jump_phys = np.zeros((0, self.G.shape[1]), dtype=np.uint8)
        self.samples = []
        self.proposals = 0
        self.accepted = 0
        for _ in range(int(burnin)):
            self._step()

    def _step(self):
        C, m = self.u.shape
        prop_u = self.u.copy()
        prop_phys = self.phys.copy()
        modes = self.rng.random(C)
        for i in range(C):
            if len(self.jump_masks) and modes[i] < self.cp_jump_prob:
                k = int(self.rng.integers(len(self.jump_masks)))
                prop_u[i] ^= self.jump_masks[k]
                prop_phys[i] ^= self.jump_phys[k]
            elif modes[i] < self.cp_jump_prob + self.pair_flip_prob and m >= 2:
                a, b = self.rng.choice(m, size=2, replace=False)
                prop_u[i, a] ^= 1; prop_u[i, b] ^= 1
                prop_phys[i] ^= self.G[a] ^ self.G[b]
            else:
                a = int(self.rng.integers(m))
                prop_u[i, a] ^= 1
                prop_phys[i] ^= self.G[a]
        n = self.ham.n_data
        prop_w = np.count_nonzero(prop_phys[:, :n] | prop_phys[:, n:], axis=1).astype(np.int16)
        delta = prop_w.astype(float) - self.w.astype(float)
        log_accept = delta * self.logr
        accept = np.log(self.rng.random(C)) < np.minimum(0.0, log_accept)
        self.u[accept] = prop_u[accept]
        self.phys[accept] = prop_phys[accept]
        self.w[accept] = prop_w[accept]
        self.proposals += C
        self.accepted += int(accept.sum())

    def ensure(self, n_samples: int):
        target = int(n_samples)
        while sum(len(x) for x in self.samples) < target:
            for _ in range(self.thin):
                self._step()
            self.samples.append(self.u.copy())
        arr = np.vstack(self.samples)
        return arr[:target].copy()

    @property
    def acceptance_rate(self):
        return self.accepted / max(self.proposals, 1)


def bar_log_ratio_from_d(d_a, d_b):
    """Bennett/bridge estimate of log(Z_b/Z_a) from equal-size samples.

    d(u)=log q_b(u)-log q_a(u). Samples d_a are evaluated on draws from p_a,
    and d_b on draws from p_b. With equal sample sizes the BAR equation is

      E_a[sigmoid(d-Delta)] = E_b[sigmoid(Delta-d)],

    where Delta=log(Z_b/Z_a).
    """
    d_a = np.asarray(d_a, dtype=float)
    d_b = np.asarray(d_b, dtype=float)
    if len(d_a) == 0 or len(d_b) == 0:
        raise ValueError("BAR needs samples from both states")
    def f(delta):
        return float(expit(d_a - delta).mean() - expit(delta - d_b).mean())
    lo, hi = -16.0, 16.0
    flo, fhi = f(lo), f(hi)
    while flo * fhi > 0 and max(abs(lo), abs(hi)) < 1024:
        lo *= 2; hi *= 2
        flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        # Pathological zero-overlap fallback: midpoint of one-sided FEP estimates.
        a = np.log(np.mean(np.exp(np.clip(d_a, -700, 700))))
        b = -np.log(np.mean(np.exp(np.clip(-d_b, -700, 700))))
        return float(0.5 * (a + b))
    return float(brentq(f, lo, hi, xtol=1e-10, rtol=1e-10, maxiter=200))


def _logq(ham, u, logr):
    return _weights(ham, u) * float(logr)


def bennett_log_masses(hams, samples_by_class, p, centers_by_class=None):
    """Estimate four log partition functions up to an additive constant.

    A QEC-aware transport map is used before BAR: for pair (a,b), stabilizer
    coordinates are XOR-shifted by the difference between CP-SB low-energy
    centers. This bijection leaves Z_b invariant but aligns dominant basins,
    greatly improving bridge overlap. Six pairwise transported BAR ratios are
    reconciled by least-squares projection onto globally consistent log masses.
    """
    logr = float(np.log((p / 3.0) / (1.0 - p)))
    if centers_by_class is None:
        centers_by_class = {c: np.zeros(hams[c].n_spins, dtype=np.uint8) for c in range(4)}
    rows = []
    vals = []
    pairs = []
    for a in range(4):
        for b in range(a + 1, 4):
            shift = np.asarray(centers_by_class[a], dtype=np.uint8) ^ np.asarray(centers_by_class[b], dtype=np.uint8)
            ua = np.asarray(samples_by_class[a], dtype=np.uint8)
            vb = np.asarray(samples_by_class[b], dtype=np.uint8)
            # Evaluate transported state b~(u)=b(u xor shift).
            d_a = _logq(hams[b], ua ^ shift[None, :], logr) - _logq(hams[a], ua, logr)
            # Samples v~p_b map to u=v xor shift, hence are samples from transported p_b~.
            ub = vb ^ shift[None, :]
            d_b = _logq(hams[b], vb, logr) - _logq(hams[a], ub, logr)
            delta = bar_log_ratio_from_d(d_a, d_b)  # log Z_b - log Z_a
            row = np.zeros(3, dtype=float)
            if b > 0: row[b - 1] += 1.0
            if a > 0: row[a - 1] -= 1.0
            rows.append(row); vals.append(delta); pairs.append((a,b,delta))
    A = np.asarray(rows, dtype=float)
    y = np.asarray(vals, dtype=float)
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    z = np.concatenate([[0.0], sol])
    z -= np.max(z)
    residual = float(np.sqrt(np.mean((A @ sol - y) ** 2)))
    return z, residual, pairs


def _bootstrap_gap(hams, samples, p, centers, rng, reps):
    n = min(len(v) for v in samples.values())
    gaps = []
    preds = []
    for _ in range(int(reps)):
        bs = {}
        for c in range(4):
            ids = rng.integers(0, n, size=n)
            bs[c] = samples[c][ids]
        z, _, _ = bennett_log_masses(hams, bs, p, centers)
        order = np.argsort(z)[::-1]
        gaps.append(float(z[order[0]] - z[order[1]]))
        preds.append(int(order[0]))
    return np.asarray(gaps), np.asarray(preds, dtype=int)


@dataclass
class BennettConfig:
    sample_schedule: tuple = (64, 128, 256)
    n_chains: int = 8
    burnin: int = 96
    thin: int = 2
    pair_flip_prob: float = 0.20
    cp_jump_prob: float = 0.15
    bootstrap_reps: int = 40
    ci_level: float = 0.90
    min_log_gap: float = 0.05



def adaptive_bennett_prediction(pools, hams, p, seed, cfg=BennettConfig()):
    samplers = {
        c: PosteriorChainSampler(
            hams[c], p, pools[c], seed=int(seed + c * 100003),
            n_chains=cfg.n_chains, burnin=cfg.burnin, thin=cfg.thin,
            pair_flip_prob=cfg.pair_flip_prob, cp_jump_prob=cfg.cp_jump_prob,
        )
        for c in range(4)
    }
    rng = np.random.default_rng(int(seed + 7000003))
    history = []
    final = None
    for n in cfg.sample_schedule:
        samples = {c: samplers[c].ensure(int(n)) for c in range(4)}
        centers = {c: np.asarray(pools[c][0], dtype=np.uint8) for c in range(4)}
        z, residual, pairs = bennett_log_masses(hams, samples, p, centers)
        order = np.argsort(z)[::-1]
        top, second = int(order[0]), int(order[1])
        gap = float(z[top] - z[second])
        gaps, preds = _bootstrap_gap(hams, samples, p, centers, rng, cfg.bootstrap_reps)
        alpha = (1.0 - float(cfg.ci_level)) / 2.0
        lo, hi = np.quantile(gaps, [alpha, 1.0 - alpha])
        top_freq = float(np.mean(preds == top))
        resolved = bool((lo > cfg.min_log_gap) and (top_freq >= cfg.ci_level))
        history.append({
            "samples_per_class": int(n),
            "prediction": top,
            "log_gap": gap,
            "bootstrap_gap_lo": float(lo),
            "bootstrap_gap_hi": float(hi),
            "bootstrap_top_frequency": top_freq,
            "bar_cycle_residual": residual,
            "resolved": int(resolved),
        })
        final = (top, z, residual, pairs, samples, gap, lo, hi, top_freq)
        if resolved:
            break
    top, z, residual, pairs, samples, gap, lo, hi, top_freq = final
    return int(top), {
        "prediction": int(top),
        "relative_log_masses": z,
        "samples_per_class": int(len(samples[0])),
        "log_gap": float(gap),
        "bootstrap_gap_lo": float(lo),
        "bootstrap_gap_hi": float(hi),
        "bootstrap_top_frequency": float(top_freq),
        "bar_cycle_residual": float(residual),
        "acceptance_rates": np.asarray([samplers[c].acceptance_rate for c in range(4)]),
        "history": history,
        "mcmc_proposals_total": int(sum(samplers[c].proposals for c in range(4))),
        "mcmc_accepted_total": int(sum(samplers[c].accepted for c in range(4))),
    }

@dataclass
class MultiStageBennettConfig:
    lambdas: tuple = (0.0, 0.25, 0.5, 0.75, 1.0)
    n_samples: int = 64
    n_chains: int = 4
    burnin: int = 64
    thin: int = 1
    pair_flip_prob: float = 0.20
    bootstrap_reps: int = 40
    refine_ratio_threshold: float = 1.20
    decision_log_margin: float = 0.03
    sign_confidence: float = 0.85
    confirm_on_disagreement: bool = True
    confirm_lambdas: tuple = (0.0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0)
    confirm_n_samples: int = 256
    confirm_n_chains: int = 8
    confirm_burnin: int = 128
    confirm_thin: int = 2
    confirm_bootstrap_reps: int = 30
    confirm_sign_confidence: float = 0.95
    confirm_log_margin: float = 0.03


def _interp_logq(hama, hamb, u, shift, logr, lam):
    u=np.asarray(u,dtype=np.uint8)
    la=_weights(hama,u)*logr
    lb=_weights(hamb,u ^ shift[None,:])*logr
    return (1.0-float(lam))*la + float(lam)*lb


def _sample_interpolated(hama,hamb,p,shift,center,lam,seed,n_samples,n_chains,burnin,thin,pair_flip_prob):
    rng=np.random.default_rng(int(seed)); m=hama.n_spins; C=int(n_chains)
    u=np.repeat(np.asarray(center,dtype=np.uint8)[None,:],C,axis=0)
    logr=float(np.log((p/3.0)/(1.0-p)))
    cur=_interp_logq(hama,hamb,u,shift,logr,lam)
    out=[]; accepted=0; proposals=0
    total_steps=int(burnin)+int(np.ceil(n_samples/C))*int(thin)
    for step in range(total_steps):
        prop=u.copy()
        pair=rng.random(C)<float(pair_flip_prob)
        for i in range(C):
            if pair[i] and m>=2:
                a,b=rng.choice(m,size=2,replace=False); prop[i,a]^=1; prop[i,b]^=1
            else:
                a=int(rng.integers(m)); prop[i,a]^=1
        val=_interp_logq(hama,hamb,prop,shift,logr,lam)
        acc=np.log(rng.random(C)) < np.minimum(0.0,val-cur)
        u[acc]=prop[acc]; cur[acc]=val[acc]
        accepted += int(acc.sum()); proposals += C
        if step>=burnin and ((step-burnin)%int(thin)==0):
            out.append(u.copy())
    arr=np.vstack(out)[:int(n_samples)]
    return arr, accepted/max(proposals,1)


def multistage_bennett_log_ratio(hama,hamb,p,center_a,center_b,seed,cfg=MultiStageBennettConfig()):
    """Transport-aligned multi-stage BAR estimate log(Z_b/Z_a).

    Intermediate geometric bridge states q_lambda ∝ q_a^(1-lambda) q_b~^lambda
    are sampled to maintain overlap. q_b~ is an XOR-transport of class b that
    aligns its CP-SB minimum with class a; the bijection preserves Z_b.
    """
    center_a=np.asarray(center_a,dtype=np.uint8); center_b=np.asarray(center_b,dtype=np.uint8)
    shift=center_a ^ center_b; logr=float(np.log((p/3.0)/(1.0-p)))
    lambdas=tuple(float(x) for x in cfg.lambdas)
    samples={}; acc={}
    for i,lam in enumerate(lambdas):
        samples[lam],acc[lam]=_sample_interpolated(
            hama,hamb,p,shift,center_a,lam,seed+10007*i,cfg.n_samples,cfg.n_chains,cfg.burnin,cfg.thin,cfg.pair_flip_prob)
    deltas=[]
    for l0,l1 in zip(lambdas[:-1],lambdas[1:]):
        ua=samples[l0]; ub=samples[l1]
        # log q_l1 - log q_l0 = (l1-l0)*(log q_b~ - log q_a)
        base_a=_weights(hamb,ua ^ shift[None,:])*logr - _weights(hama,ua)*logr
        base_b=_weights(hamb,ub ^ shift[None,:])*logr - _weights(hama,ub)*logr
        d0=(l1-l0)*base_a; d1=(l1-l0)*base_b
        deltas.append(bar_log_ratio_from_d(d0,d1))
    delta=float(np.sum(deltas))
    # Bootstrap samples within each bridge state for sign confidence.
    rng=np.random.default_rng(int(seed+9000001)); boots=[]
    for _ in range(int(cfg.bootstrap_reps)):
        bd=[]
        res={lam:samples[lam][rng.integers(0,len(samples[lam]),size=len(samples[lam]))] for lam in lambdas}
        for l0,l1 in zip(lambdas[:-1],lambdas[1:]):
            ua=res[l0]; ub=res[l1]
            ba=_weights(hamb,ua ^ shift[None,:])*logr - _weights(hama,ua)*logr
            bb=_weights(hamb,ub ^ shift[None,:])*logr - _weights(hama,ub)*logr
            bd.append(bar_log_ratio_from_d((l1-l0)*ba,(l1-l0)*bb))
        boots.append(float(np.sum(bd)))
    boots=np.asarray(boots)
    sign_conf=float(max(np.mean(boots>0),np.mean(boots<0)))
    lo,hi=np.quantile(boots,[0.05,0.95])
    return delta,{
        'log_ratio_b_over_a':delta,'bootstrap_lo90':float(lo),'bootstrap_hi90':float(hi),
        'sign_confidence':sign_conf,'mean_acceptance':float(np.mean(list(acc.values()))),
        'stage_deltas':tuple(float(v) for v in deltas),'n_bridge_states':len(lambdas),
        'samples_per_bridge_state':int(cfg.n_samples),
    }


def hybrid_bennett_refine(p41_class, local_scores, pools, hams, p, seed, cfg=MultiStageBennettConfig()):
    order=np.argsort(np.asarray(local_scores,dtype=float))[::-1]
    a,b=int(order[0]),int(order[1])
    ratio=float(local_scores[a]/max(float(local_scores[b]),1e-300))
    if ratio>=float(cfg.refine_ratio_threshold):
        return int(p41_class),{'refined':0,'confirmed':0,'local_top2_ratio':ratio,'pair':(a,b),'fallback':'local_confident'}
    delta,diag=multistage_bennett_log_ratio(hams[a],hams[b],p,pools[a][0],pools[b][0],seed,cfg)
    if diag['sign_confidence']>=cfg.sign_confidence and abs(delta)>=cfg.decision_log_margin:
        low_pred=b if delta>0 else a
    else:
        low_pred=int(p41_class)
    diag.update({'refined':1,'confirmed':0,'local_top2_ratio':ratio,'pair':(a,b)})

    # A low-budget BAR disagreement is exactly where Monte-Carlo overlap error is
    # most dangerous. Escalate only these rare cases; require the high-budget
    # transported BAR to confirm the same direction before overriding P4.1.
    if int(low_pred) != int(p41_class) and cfg.confirm_on_disagreement:
        ccfg=MultiStageBennettConfig(
            lambdas=cfg.confirm_lambdas,
            n_samples=cfg.confirm_n_samples,
            n_chains=cfg.confirm_n_chains,
            burnin=cfg.confirm_burnin,
            thin=cfg.confirm_thin,
            pair_flip_prob=cfg.pair_flip_prob,
            bootstrap_reps=cfg.confirm_bootstrap_reps,
            refine_ratio_threshold=cfg.refine_ratio_threshold,
            decision_log_margin=cfg.confirm_log_margin,
            sign_confidence=cfg.confirm_sign_confidence,
            confirm_on_disagreement=False,
        )
        delta2,d2=multistage_bennett_log_ratio(
            hams[a],hams[b],p,pools[a][0],pools[b][0],seed+17000011,ccfg
        )
        if d2['sign_confidence']>=cfg.confirm_sign_confidence and abs(delta2)>=cfg.confirm_log_margin:
            high_pred=b if delta2>0 else a
        else:
            high_pred=int(p41_class)
        confirmed=(int(high_pred)==int(low_pred))
        pred=int(low_pred) if confirmed else int(p41_class)
        diag.update({
            'confirmed':1,'confirmation_agrees':int(confirmed),
            'confirm_log_ratio_b_over_a':float(delta2),
            'confirm_bootstrap_lo90':float(d2['bootstrap_lo90']),
            'confirm_bootstrap_hi90':float(d2['bootstrap_hi90']),
            'confirm_sign_confidence':float(d2['sign_confidence']),
            'confirm_mean_acceptance':float(d2['mean_acceptance']),
            'fallback':'bar_confirmed' if confirmed else 'p41_after_disagreement_check',
            'prediction':pred,
        })
        return pred,diag

    pred=int(low_pred)
    diag.update({'fallback':'bar_decisive' if pred!=int(p41_class) else 'p41_fallback','prediction':pred})
    return pred,diag

