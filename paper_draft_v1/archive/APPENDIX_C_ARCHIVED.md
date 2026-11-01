# Archived Appendix C: structural circuit-level checks

Project record only; excluded from the current submission and its Supplemental Material. The source text below is the complete pre-relocation English Appendix C. Its relative data links retain their original manuscript-level meaning.

## Appendix C. Structural circuit-level checks

### C.1. Fault-space coordinates and construction

These checks examine whether fault-space coordinates can preserve detector and hypothesized logical-observable records; they do not evaluate a circuit-level decoder. Let H and L denote the corresponding binary maps, and let \(A=[H;L]\) denote their vertical concatenation. All operations are over \(\mathbb F_2\). For a fault-mechanism vector f, the records are \(Hf^T\) and \(Lf^T\). Generators satisfying \(AG^T=0\) preserve both records under \(f\mapsto f\oplus uG\). A complete independent basis spans \(\ker A\), with dimension \(N_{\rm mech}-\operatorname{rank}_{\mathbb F_2}A\). Here \(N_{\rm mech}\) is the number of fault mechanisms, \(H\in\mathbb F_2^{N_{\rm det}\times N_{\rm mech}}\), and \(L\in\mathbb F_2^{N_{\rm obs}\times N_{\rm mech}}\), where \(N_{\rm det}\) and \(N_{\rm obs}\) count detectors and tracked logical observables. If the kernel dimension is \(M\), then \(G\in\mathbb F_2^{M\times N_{\rm mech}}\). Here L is a response matrix, not the logical-class label used in the main text. These are fault-space coordinates, rather than the data-qubit stabilizer coordinates of the main benchmark.

The archived instances are memory-Z experiments: rotated surface codes with d rounds at d=3,5,7,9,11,13, and a \([[72,12,6]]\) bivariate-bicycle (BB) circuit with six rounds of the depth-eight syndrome-extraction cycle associated with Ref. [34]. The circuit noise parameter is p=0.003. The surface-code circuits use Stim's generated reset, measurement, post-Clifford depolarization, and pre-round data-depolarization settings [35]. The implemented BB circuit includes post-CNOT two-qubit depolarization, single-qubit depolarization on idle data during gate layers, post-reset and pre-measurement flips, and end-of-cycle data depolarization, with initialization and final readout handled explicitly. The extraction-cycle depth does not include all experiment boundaries. The circuit parameter is not a uniform probability for the merged detector-error mechanisms.

Low-weight candidates are generated heuristically by applying BP+OSD to individual mechanism responses under the stacked map, suppressing the mechanism being reconstructed and rejecting candidates that fail the kernel constraint. Independent candidates are selected over \(\mathbb F_2\), and remaining directions are supplied by global null-space completion. The candidate search imposes no fixed spatial or temporal window. The code checks the target rank and kernel membership, but the construction is not a minimum-weight-basis algorithm. Here “global” identifies the completion procedure, not a guarantee that every completion element spans the circuit. Generator weight counts detector-error mechanisms, not physical data-qubit Pauli weight.

### C.2. Observed basis structure

Table C1 summarizes the archived reports. In the selected surface-code bases, the generators retained before completion have median weight three and 99th-percentile weight at most four. The construction adds respectively 0, 4, 6, 6, 13, and 14 completion elements. These counts describe the generated candidate sets, rather than a lower bound for every possible basis. For BB, it adds 261 elements to a basis of dimension 15,726.

**Table C1. Archived circuit-level basis statistics.** Weight percentiles refer only to selected generators before completion; basis dimension and maximum mechanism incidence refer to the completed basis. Incidence counts generator rows containing a given mechanism. The dash denotes no completion elements. The noninteger 99th percentile at d=3 results from percentile interpolation; individual weights are integers.

| Instance | Mechanisms | Basis dimension | Pre-completion weight median / 99th percentile | Completion count | Maximum completion weight | Maximum mechanism incidence |
|---|---:|---:|---:|---:|---:|---:|
| Surface, d=3 | 219 | 194 | 3 / 3.07 | 0 | — | 21 |
| Surface, d=5 | 1,677 | 1,556 | 3 / 4 | 4 | 14 | 24 |
| Surface, d=7 | 5,471 | 5,134 | 3 / 4 | 6 | 26 | 26 |
| Surface, d=9 | 12,705 | 11,984 | 3 / 4 | 6 | 30 | 27 |
| Surface, d=11 | 24,483 | 23,162 | 3 / 4 | 13 | 36 | 29 |
| Surface, d=13 | 41,909 | 39,724 | 3 / 4 | 14 | 57 | 26 |
| BB \([[72,12,6]]\) | 16,164 | 15,726 | 3 / 4 | 261 | 104 | 95 |

At distance 13, ten completion elements touch all fourteen detector-time layers and one touches thirteen layers; the maximum completion weight is 57 within a basis of 39,724 elements. These counts refer to the union of detector coordinates touched by constituent mechanisms before their responses cancel modulo two, not to a nonzero net syndrome or a fault duration. The BB maximum completion weight is 104. Its basis uses greedy lightening, whereas the surface-code runs additionally use heuristic search; neither procedure certifies minimum weights.

Generator row weights and mechanism incidence are distinct. The maximum incidence is 21–29 in the surface-code instances and 95 in BB. Under an independent-mechanism likelihood, incidence gives the spin order of each mechanism's energy term before combining like terms. Small row weights alone therefore do not establish low-order interactions or the bounded support and overlap conditions used in Appendix A. These finite-instance measurements establish neither asymptotic locality nor its impossibility for alternative bases.

### C.3. Scope and reproducibility limits

No circuit-level decoding LER or online decoder latency was measured in these structural checks. The results motivate a fault-space construction but do not validate likelihood aggregation, search performance, or end-to-end complexity. In particular, nonuniform mechanism weights and global basis elements require treatment before the main search-and-aggregation approach can be assessed at circuit level.

The [summary table](../../results/circuit_locality/summary/locality_summary.csv), [distance-13 inspection report](../../results/circuit_locality/sc13_p0.0030_r13/T1_INSPECT.md), and [BB construction report](../../results/circuit_locality/bb72B_p0.0030_r6/T1_REPORT.md) retain the reported structural evidence. The [construction script](../../code/scripts/run_circuit_locality.py) and [kernel-basis implementation](../../code/qec_sb/circuit_level/kernel_basis.py) specify the procedure and its internal checks. The present project snapshot does not include the underlying circuit, detector-error model, or matrix artifacts. Consistency of the retained reports can therefore be checked, but independent recomputation of the archived ranks, kernel membership, and full supports requires recovering or regenerating those artifacts. A new heuristic run need not reproduce the same basis element by element.

