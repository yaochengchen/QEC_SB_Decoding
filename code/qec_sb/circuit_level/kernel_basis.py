"""Test 1: a sparse, local generating set of the DEM stabilizer group  ker([H; L]).

Why: CP-SB needs generator variables u with  e(u) = e0 XOR u.G  where every row of G is a
low-weight set of error mechanisms. The old Core prototype used a plain GF(2) null-space
basis (row weight ~ hundreds) and lost to PyMatching. Here each generator is found by
*self-decoding*: for mechanism i, decode its own detector pattern D_i with i excluded;
the min-weight completion T_i gives the kernel element g_i = {i} U T_i. BP+OSD on [H; L]
(observable rows included) keeps the logical class fixed, so g_i is a stabilizer element.

Steps
  1. candidates: 3 decodes per mechanism (uniform weights, true LLR weights, jittered)
  2. greedy GF(2) basis by ascending weight (incremental elimination on int bitsets)
  3. completion: if rank < n_mech - rank([H;L]), fill from the global null space and REPORT
     how many heavy completions were needed -> this number is the locality verdict.
  4. logical generators: null(H) modulo the stabilizer basis (k of them).
"""
from __future__ import annotations

import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp

from . import gf2
from .dem import DEMMatrices

_WORKER = {}


def _import_bposd():
    """ldpc's package __init__ imports sinter, which is not needed here and may be missing (--no-deps install)."""
    try:
        from ldpc import BpOsdDecoder
    except ModuleNotFoundError as exc:
        if exc.name != "sinter":
            raise
        import sys
        import types
        stub = types.ModuleType("sinter")
        stub.Decoder = object
        stub.CompiledDecoder = object
        sys.modules["sinter"] = stub
        from ldpc import BpOsdDecoder
    return BpOsdDecoder



# ---------------------------------------------------------------- candidates (per worker)
def _init_worker(HL_csr, weights, osd_order, max_iter):
    BpOsdDecoder = _import_bposd()
    n = HL_csr.shape[1]
    base_p = 1.0 / (1.0 + np.exp(weights))                # p_i from w_i
    _WORKER["HL"] = HL_csr
    _WORKER["n"] = n
    _WORKER["profiles"] = {
        "uniform": np.full(n, 0.01),
        "llr": np.clip(base_p, 1e-9, 0.49),
    }
    _WORKER["dec"] = {
        name: BpOsdDecoder(HL_csr, channel_probs=prob, max_iter=max_iter, bp_method="minimum_sum",
                           ms_scaling_factor=0.9, osd_method="osd_cs", osd_order=osd_order)
        for name, prob in _WORKER["profiles"].items()
    }
    _WORKER["seed"] = 0


def _self_decode_chunk(idx_chunk, jitter_trials):
    """For each mechanism i in chunk: kernel elements {i} U T_i from several weight profiles."""
    HL = _WORKER["HL"]
    n = _WORKER["n"]
    HLc = HL.tocsc()
    out = []          # (mechanism, weight, support tuple)
    rng = np.random.default_rng(int(idx_chunk[0]) + 7919)
    profiles = [(name, dec, _WORKER["profiles"][name]) for name, dec in _WORKER["dec"].items()]
    for trial in range(jitter_trials):                      # one jittered decoder per chunk
        base = _WORKER["profiles"]["llr"] * np.exp(rng.normal(0, 0.5, n))
        dec, jit_probs = _make_jit_decoder(base)
        profiles.append((f"jit{trial}", dec, jit_probs))
    for i in idx_chunk:
        col = HLc.indices[HLc.indptr[i]:HLc.indptr[i + 1]]
        s = np.zeros(HL.shape[0], dtype=np.uint8)
        s[col] = 1
        for name, dec, prob in profiles:
            p_excl = prob.copy()
            p_excl[i] = 1e-12                               # exclude i itself
            _set_probs(dec, p_excl)
            t = np.asarray(dec.decode(s), dtype=np.uint8)
            if t[i]:
                continue                                    # decoder fell back on i: no completion
            sup = np.nonzero(t)[0]
            # verify exactly: HL (t + e_i) = 0
            chk = (np.asarray(HL[:, sup].sum(axis=1)).ravel() + s) % 2
            if chk.any():
                continue
            g = tuple(sorted(sup.tolist() + [int(i)]))
            out.append((int(i), len(g), g))
    return out


def _set_probs(dec, p):
    try:
        dec.update_channel_probs(p)
    except AttributeError:
        dec.channel_probs = p


def _make_jit_decoder(prob):
    BpOsdDecoder = _import_bposd()
    probs = np.clip(prob, 1e-9, 0.49)
    dec = BpOsdDecoder(_WORKER["HL"], channel_probs=probs, max_iter=30,
                       bp_method="minimum_sum", ms_scaling_factor=0.9, osd_method="osd_cs", osd_order=5)
    return dec, probs        # Cython decoder objects cannot carry extra attributes


# ---------------------------------------------------------------- main routine
@dataclass
class GeneratorSet:
    G: sp.csr_matrix                 # (m, n_mech) stabilizer generators
    Lg: sp.csr_matrix                # (k, n_mech) logical generators (null(H) mod stabilizers)
    weights: np.ndarray              # row weights of G
    n_completion: int                # generators that had to come from the global null space
    completion_weights: list = field(default_factory=list)
    self_found_rate: float = 0.0     # fraction of mechanisms whose own decode yielded a kernel element
    stats: dict = field(default_factory=dict)


def build_generator_set(dm: DEMMatrices, workers: int = 32, osd_order: int = 7, max_iter: int = 30,
                        jitter_trials: int = 1, chunk: int = 200, log=print, completion: str = "full") -> GeneratorSet:
    """completion = "full": lighten each null-space completion by a randomised local information-set search
    (that search took ~2 h of the 2.2 h T1 wall time on bb72 in Phase 0); "fast": greedy pairwise lightening only.
    The completions are reported separately either way (they are the 'global elements', see t1_inspect.py)."""
    t0 = time.time()
    HL = dm.HL()
    n = dm.n_mech
    rows_HL = gf2.rows_from_csr(HL)
    rank_HL = gf2.rank(rows_HL)
    target = n - rank_HL
    log(f"[kb] n_mech={n} n_det={dm.n_det} n_obs={dm.n_obs} rank(HL)={rank_HL} -> stabilizer dim {target}")

    # 1. candidates
    chunks = [np.arange(a, min(a + chunk, n)) for a in range(0, n, chunk)]
    cands: set[tuple] = set()                       # supports as tuples (small); ints built lazily
    found_for = np.zeros(n, dtype=bool)
    with ProcessPoolExecutor(workers, initializer=_init_worker,
                             initargs=(HL, dm.weights, osd_order, max_iter)) as ex:
        for res in ex.map(_self_decode_chunk, chunks, [jitter_trials] * len(chunks), chunksize=1):
            for i, w, g in res:
                found_for[i] = True
                cands.add(g)
    log(f"[kb] candidates: {len(cands)} unique, self-found rate {found_for.mean():.4f}, {time.time()-t0:.0f}s")

    # 2. greedy basis by weight
    basis = gf2.Basis()
    chosen: list[int] = []
    for g in sorted(cands, key=len):
        v = 0
        for j in g:
            v |= 1 << j
        if basis.insert(v):
            chosen.append(v)
        if basis.rank == target:
            break
    log(f"[kb] greedy basis rank {basis.rank}/{target}, {time.time()-t0:.0f}s")

    # 3. completion from the global null space (lazy, only if needed)
    completion_w: list[int] = []
    if basis.rank < target:
        piv = gf2.rref(rows_HL)
        pcols = sorted(piv)
        pset = set(pcols)
        R = gf2.dense_from_rows([piv[c] for c in pcols], n)
        pc = np.asarray(pcols)
        for f in range(n):
            if f in pset:
                continue
            v = 1 << f
            for c in pc[R[:, f].astype(bool)]:
                v |= 1 << int(c)
            if basis.reduce(v) == 0:                 # already in the span of the sparse basis
                continue
            v = _lighten(v, chosen, passes=6)         # same coset (chosen subset of the span), lighter
            if completion == "full":
                v = _lighten_search(v, chosen, basis, n)  # randomised information-set search in the local span
            ok = basis.insert(v)                      # proper insertion: reduces, never overwrites a pivot
            assert ok, "lightened completion vector became dependent"
            chosen.append(v)
            completion_w.append(gf2.popcount(v))
            if len(completion_w) % 20 == 0:
                log(f"[kb]   {len(completion_w)} completions, rank {basis.rank}/{target}, {time.time()-t0:.0f}s")
            if basis.rank == target:
                break
        log(f"[kb] completed with {len(completion_w)} null-space vectors, weights "
            f"{sorted(completion_w)[:10]}... max {max(completion_w) if completion_w else 0}")
    assert basis.rank == target, (basis.rank, target)

    # 4. logical generators: solve [H; L] v = [0; e_a] for each observable a
    lg = []
    for a in range(dm.n_obs):
        rhs = np.zeros(len(rows_HL), dtype=np.uint8)
        rhs[dm.n_det + a] = 1
        v = gf2.solve(rows_HL, rhs, n)
        if v is None:
            log(f"[kb] WARNING: observable {a} has no pure-logical representative")
            continue
        lg.append(_lighten(v, chosen))

    G = _csr_from_ints(chosen, n)
    Lg = _csr_from_ints(lg, n)
    w = np.asarray(G.sum(axis=1)).ravel()
    colw = np.asarray(G.sum(axis=0)).ravel()
    n_self = len(chosen) - len(completion_w)                      # self-decoded generators come first
    w_self = w[:n_self] if n_self else w
    stats = {
        "n_mech": n, "n_det": dm.n_det, "n_obs": dm.n_obs, "rank_HL": rank_HL, "dim_stabilizer": target,
        "gen_weight_mean": float(w.mean()), "gen_weight_p50": float(np.percentile(w, 50)),
        "gen_weight_p99": float(np.percentile(w, 99)), "gen_weight_max": int(w.max()),
        "selfdecoded_weight_max": int(w_self.max()), "selfdecoded_frac_weight_le8": float((w_self <= 8).mean()),
        "frac_weight_le8": float((w <= 8).mean()), "frac_weight_le12": float((w <= 12).mean()),
        "col_weight_mean": float(colw.mean()), "col_weight_max": int(colw.max()),
        "mech_uncovered": int((colw == 0).sum()),
        "n_completion": len(completion_w), "completion_weight_max": int(max(completion_w)) if completion_w else 0,
        "self_found_rate": float(found_for.mean()), "n_candidates": len(cands),
        "logical_weights": [int(gf2.popcount(v)) for v in lg], "seconds": time.time() - t0,
        "completion_mode": completion,
    }
    bad = int(((HL.astype(np.int64) @ G.T.astype(np.int64)).toarray() % 2).any(axis=0).sum())
    assert bad == 0, f"{bad} generators are not in ker([H;L])"
    return GeneratorSet(G, Lg, w, len(completion_w), completion_w, float(found_for.mean()), stats)


def _lighten(v: int, sparse_vectors: list[int], passes: int = 2) -> int:
    """Greedy pairwise reduction of v against sparse vectors sharing support (keeps v in the same coset)."""
    for _ in range(passes):
        improved = False
        wv = gf2.popcount(v)
        for g in sparse_vectors:
            if v & g:
                w2 = gf2.popcount(v ^ g)
                if w2 < wv:
                    v ^= g
                    wv = w2
                    improved = True
        if not improved:
            break
    return v


def _lighten_search(v: int, chosen: list[int], basis: "gf2.Basis", n: int, trials: int = 40, hops: int = 2) -> int:
    """Lightest vector found in span({v} U {sparse generators within `hops` of v}) that is still
    independent of the current basis (i.e. a valid replacement for v). Falls back to v."""
    rng = np.random.default_rng(gf2.popcount(v) + 12345)
    sup = v
    local = []
    for _ in range(hops):
        local = [g for g in chosen if g & sup]
        for g in local:
            sup |= g
        if len(local) > 400:
            break
    found = gf2.low_weight_vectors([v] + local, n, trials, rng, max_weight=gf2.popcount(v))
    best, bw = v, gf2.popcount(v)
    for x, w in sorted(found.items(), key=lambda kv: kv[1]):
        if w >= bw:
            break
        if basis.reduce(x):
            best, bw = x, w
            break
    return best


def _csr_from_ints(vecs: list[int], n: int) -> sp.csr_matrix:
    indptr = [0]
    indices = []
    for v in vecs:
        indices.extend(gf2.support(v))
        indptr.append(len(indices))
    return sp.csr_matrix((np.ones(len(indices), dtype=np.uint8), indices, indptr), shape=(len(vecs), n))


def verdict(stats: dict, max_ok=12, frac_le8_ok=0.99, completion_ok=0.01) -> tuple[bool, str]:
    """Pre-declared Test-1 gate."""
    reasons = []
    if stats["gen_weight_max"] > max_ok:
        reasons.append(f"max generator weight {stats['gen_weight_max']} > {max_ok}")
    if stats["frac_weight_le8"] < frac_le8_ok:
        reasons.append(f"only {stats['frac_weight_le8']:.3f} of generators have weight <= 8")
    if stats["n_completion"] > completion_ok * stats["dim_stabilizer"]:
        reasons.append(f"{stats['n_completion']} heavy null-space completions (> {completion_ok:.0%})")
    note = ""
    if reasons and stats.get("selfdecoded_weight_max", 99) <= max_ok and stats.get("selfdecoded_frac_weight_le8", 0) >= frac_le8_ok \
            and stats["n_completion"] <= completion_ok * stats["dim_stabilizer"]:
        note = " (would PASS if the few null-space completions were excluded from the weight gate)"
    return (not reasons), ("PASS" if not reasons else "FAIL: " + "; ".join(reasons) + note)
