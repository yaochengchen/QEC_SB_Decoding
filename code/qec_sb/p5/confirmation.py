"""Fixed, held-out Paper 1 confirmation with complete-family inference."""

from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from ..p42a import run_p42a
from .config import build_frozen_p42a_config, frozen_config_sha256
from .provenance import source_manifest
from .statistics import adjust_p_values, exact_paired_statistics, wilson_interval


DEFAULT_BASE_SEED = 2027101300
METHODS = ("MWPM", "CP_SB_MAP", "CP_local_SOE", "P4_1_importance", "P4_2_Bennett")
REQUIRED_ARTIFACTS = ("metadata.json", "shot_level.csv", "summary.csv",
                      "syndrome_level.csv", "paired_comparisons.csv")


@dataclass(frozen=True)
class ConfirmationBatch:
    batch_id: str
    distance: int
    p: float
    replicate: int
    seed: int
    solver_seed: int
    role: str
    family_index: int
    expected_shots: int


def build_confirmation_plan(base_seed: int = DEFAULT_BASE_SEED) -> tuple[ConfirmationBatch, ...]:
    # Reserve a separate seed namespace, above all retained exploratory grid seeds.
    if isinstance(base_seed, bool) or not isinstance(base_seed, int) or not 1_000_000_000 <= base_seed <= 2**32 - 1_000_012:
        raise ValueError("base seed must be an integer in [1000000000, 4293967284]")
    batches = []
    for distance, shots in ((5, 1500), (7, 1000), (9, 500)):
        for probability in (.10, .12):
            family_index = len(batches) // 2
            for replicate in range(2):
                offset = 2 * family_index + replicate
                batches.append(ConfirmationBatch(
                    batch_id=f"d{distance}_p{probability:.3f}_r{replicate}",
                    distance=distance, p=probability, replicate=replicate,
                    seed=base_seed + offset, solver_seed=base_seed + 1_000_000 + offset,
                    role="primary" if distance == 5 else "secondary",
                    family_index=family_index, expected_shots=shots,
                ))
    return tuple(batches)


def validate_confirmation_registry(records, plan, *, allow_incomplete=False):
    if not isinstance(records, list) or any(not isinstance(row, dict) for row in records):
        raise ValueError("registry must be a list of batch records")
    counts = Counter(row.get("batch_id") for row in records)
    expected = {batch.batch_id: batch for batch in plan}
    if any(count > 1 for count in counts.values()):
        raise ValueError("duplicate registered batches")
    if set(counts) - set(expected):
        raise ValueError("unexpected registered batches")
    if set(expected) - set(counts):
        raise ValueError("missing registered batches")
    for row in records:
        for key, value in asdict(expected[row["batch_id"]]).items():
            if row.get(key) != value:
                raise ValueError(f"registry plan mismatch: {row['batch_id']} {key}")
        allowed = {"pending", "running", "failed", "completed"} if allow_incomplete else {"completed"}
        if row.get("status") not in allowed:
            raise ValueError(f"batch {row['batch_id']} has status {row.get('status')}")
    for key, message in (("frozen_config_sha256", "mixed-config"),
                         ("experiment_source_manifest_sha256", "mixed-generating-source")):
        values = {row.get(key) for row in records}
        if len(values) != 1 or any(not isinstance(value, str) or len(value) != 64 for value in values):
            raise ValueError(f"{message} hashes in registry")
    return records


def _json(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"missing or invalid JSON artifact: {path}") from exc


def _write_json(path, payload):
    _write_text(path, json.dumps(payload, indent=2))


def _write_text(path, content):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def _write_bytes(path, content):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_bytes(content)
    temporary.replace(path)


def _canonical(value):
    return json.loads(json.dumps(value))


def _config(batch):
    return build_frozen_p42a_config(batch.distance, batch.p, batch.expected_shots,
                                   batch.seed, batch.solver_seed)


def _hash_artifacts(batch_dir):
    result = {}
    for name in REQUIRED_ARTIFACTS:
        path = batch_dir / name
        if not path.is_file():
            raise ValueError(f"missing required artifact: {path}")
        result[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def _require_columns(frame, names, label):
    missing = set(names) - set(frame.columns)
    if missing:
        raise ValueError(f"missing {label} columns: {sorted(missing)}")


def _validate_batch(batch, batch_dir, record):
    hashes = _hash_artifacts(batch_dir)
    if hashes != record.get("artifact_sha256"):
        raise ValueError(f"artifact hash mismatch: {batch.batch_id}")
    metadata = _json(batch_dir / "metadata.json")
    if metadata.get("config") != _canonical(asdict(_config(batch))):
        raise ValueError(f"metadata config mismatch: {batch.batch_id}")
    for key in ("frozen_config_sha256", "experiment_source_manifest_sha256"):
        if metadata.get(key) != record[key]:
            raise ValueError(f"metadata {key} hash mismatch: {batch.batch_id}")
    if metadata.get("exact_enumeration_used") is not False or metadata.get("oracle_used") is not False:
        raise ValueError("confirmation metadata must exclude exact enumeration and oracle use")
    try:
        shots = pd.read_csv(batch_dir / "shot_level.csv", dtype={"syndrome_key": str})
        summary = pd.read_csv(batch_dir / "summary.csv")
        syndromes = pd.read_csv(batch_dir / "syndrome_level.csv", dtype={"syndrome_key": str})
        comparisons = pd.read_csv(batch_dir / "paired_comparisons.csv")
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise ValueError(f"invalid CSV artifact: {batch.batch_id}") from exc
    _require_columns(shots, ["shot", "syndrome_key", "truth_class"] +
                     [f"{kind}_{method}" for method in METHODS for kind in ("pred", "correct")], "shot")
    if len(shots) != batch.expected_shots or shots.shot.tolist() != list(range(batch.expected_shots)):
        raise ValueError(f"shot count/order mismatch: {batch.batch_id}")
    if not shots.truth_class.isin(range(4)).all():
        raise ValueError("invalid truth classes")
    for method in METHODS:
        correct, prediction = shots[f"correct_{method}"], shots[f"pred_{method}"]
        if not correct.isin((0, 1)).all() or not prediction.isin(range(4)).all():
            raise ValueError("invalid binary correctness or prediction")
        if not np.array_equal(correct, (prediction == shots.truth_class).astype(int)):
            raise ValueError("correctness disagrees with paired predictions")
    _require_columns(summary, ("method", "shots", "distance", "p", "logical_success", "logical_error"), "summary")
    if Counter(summary.method) != Counter(METHODS):
        raise ValueError("summary methods missing, duplicate, or unexpected")
    for row in summary.to_dict("records"):
        if row["shots"] != batch.expected_shots or row["distance"] != batch.distance or row["p"] != batch.p:
            raise ValueError("summary config/shot mismatch")
        success = shots[f"correct_{row['method']}"].mean()
        if not np.isclose(row["logical_success"], success, atol=1e-12, rtol=0) or not np.isclose(row["logical_error"], 1-success, atol=1e-12, rtol=0):
            raise ValueError("summary disagrees with shot outcomes")
    _require_columns(syndromes, ("syndrome_key",), "syndrome")
    if syndromes.syndrome_key.duplicated().any() or set(syndromes.syndrome_key) != set(shots.syndrome_key):
        raise ValueError("syndrome coverage mismatch")
    if metadata.get("unique_syndromes") != len(syndromes):
        raise ValueError("metadata syndrome count mismatch")
    _require_columns(comparisons, ("method_a", "method_b", "accuracy_diff_a_minus_b"), "paired")
    if comparisons.empty:
        raise ValueError("missing paired comparisons")
    return shots


def _protocol(plan, base_seed, source_hash):
    rows = "\n".join(f"| {b.batch_id} | {b.role} | {b.family_index} | {b.expected_shots} | {b.seed} | {b.solver_seed} |" for b in plan)
    return f"""# Paper 1 independent confirmation protocol

Predeclared before decoding: six cells, twelve batches, 12000 held-out shots.
Primary: d=5, p=.10/.12, 2 x 1500 shots per cell. Secondary: d=7,
p=.10/.12, 2 x 1000; d=9, p=.10/.12, 2 x 500. Budgets are fixed.
Base seed: {base_seed}. Sample seeds are base + 2*family_index + replicate;
solver seeds add 1000000. The reserved seed namespace starts at 1000000000,
above the retained historical grids. Family indices run from 0 through 5.

frozen_config_sha256: {frozen_config_sha256()}
Generating source_manifest_sha256: {source_hash}

Only build_frozen_p42a_config and run_p42a generate decoder outcomes.
No hyperparameter tuning, exact logical-mass enumeration, or decoder oracle.
The source manifest and all batch identities are registered before the first run.

| Batch | Role | Family index | Expected shots | Sample seed | Solver seed |
|---|---|---:|---:|---:|---:|
{rows}

## Analysis and stopping rule

Pool both replicates within each cell, retaining all shot pairs. Report logical
success/error and Wilson 95% intervals (ci95 refers to success); paired differences
use exact two-sided McNemar/sign tests and conservative simultaneous
Clopper-Pearson difference intervals. Report wins/losses and relative logical
error reduction. Holm and Bonferroni both correct across the same six predeclared
P4_2_Bennett-versus-MWPM cell tests. The primary label does not create a separate
two-test correction or permit selective reporting. All null/negative outcomes stay.
PyMatching/BP_Matching comparisons, if available, are supplementary unadjusted tests
outside the predeclared MWPM family; no confirmatory claim follows from them.

Stopping is allowed only for implementation failure, never statistical outcomes.
Failed/running attempts and their files are retained on resume; validated completed
batches are skipped only with identical source, frozen config, seeds and budgets.
Aggregation refuses incomplete, duplicate, unexpected, mixed-config or mixed-source
families. Historical exploratory data are never pooled into confirmation p-values.

## Paired baseline interface

Use scripts/run_p5_baselines.py --grid <this-directory> --out <this-directory>/baselines,
then run scripts/run_paper1_confirmation.py --out <this-directory> --base-seed {base_seed}
--aggregate-only. The baseline runner verifies the replayed truth and syndrome keys.
The confirmation aggregator additionally checks complete batch/shot identity and
the saved final-decoder predictions before accepting any supplementary comparisons.
The original baseline grid path is provenance only; these content checks also work
after copying the directory or extracting an archive at a different location.

## Published result generations

results_manifest.json is the atomic publication point. Its relative paths and
SHA-256 digests identify an immutable pair of CSVs under aggregates/<generation>/.
Use load_confirmation_results to read and verify one committed generation. Root
pooled_results.csv and paired_results.csv are compatibility copies, replaced from
same-directory temporary files under the output lock. A process crash between
those replacements can leave mixed compatibility copies; the manifest still
selects the previous intact pair until the next generation is fully committed.
"""


def _validated_manifest(outdir, plan, base_seed, current_source=None):
    manifest = _json(outdir / "confirmation_manifest.json")
    if manifest.get("base_seed") != base_seed or manifest.get("plan") != [asdict(batch) for batch in plan]:
        raise ValueError("confirmation manifest plan/base seed mismatch")
    if manifest.get("frozen_config_sha256") != frozen_config_sha256():
        raise ValueError("confirmation frozen config hash mismatch")
    source_hash = manifest.get("source_manifest_sha256")
    if current_source is not None and source_hash != current_source:
        raise ValueError("current generating source hash differs from registered source")
    return manifest


def _load_registered(outdir, plan, base_seed, *, allow_incomplete=False, current_source=None):
    manifest = _validated_manifest(outdir, plan, base_seed, current_source)
    source_hash = manifest["source_manifest_sha256"]
    records = validate_confirmation_registry(_json(outdir / "batch_registry.json"), plan,
                                              allow_incomplete=allow_incomplete)
    for row in records:
        if row["frozen_config_sha256"] != manifest["frozen_config_sha256"] or row["experiment_source_manifest_sha256"] != source_hash:
            raise ValueError("registry source/frozen hash differs from confirmation manifest")
    if (outdir / "CONFIRMATION_PROTOCOL.md").read_text(encoding="utf-8") != _protocol(plan, base_seed, source_hash):
        raise ValueError("registered confirmation protocol mismatch")
    directories = outdir / "batches"
    if directories.exists() and {p.name for p in directories.iterdir()} - {b.batch_id for b in plan}:
        raise ValueError("unexpected batch artifact directory")
    return records


def _initialize_registration(outdir, plan, base_seed, provenance):
    """Finish a registration interrupted before any batch could have started."""
    manifest_path = outdir / "confirmation_manifest.json"
    protocol_path = outdir / "CONFIRMATION_PROTOCOL.md"
    registry_path = outdir / "batch_registry.json"
    names = {path.name for path in outdir.iterdir()}
    allowed = {".confirmation.lock", manifest_path.name, protocol_path.name,
               *(path.name + ".tmp" for path in (manifest_path, protocol_path, registry_path))}
    if names - allowed:
        raise ValueError("unregistered artifacts in confirmation output directory")
    source_hash = provenance["source_manifest_sha256"]
    if manifest_path.exists():
        _validated_manifest(outdir, plan, base_seed, source_hash)
    else:
        if names - {".confirmation.lock", manifest_path.name + ".tmp"}:
            raise ValueError("missing confirmation manifest for partial registration")
        _write_json(manifest_path, {
            "base_seed": base_seed, "plan": [asdict(batch) for batch in plan],
            "frozen_config_sha256": frozen_config_sha256(), **provenance,
        })
    protocol = _protocol(plan, base_seed, source_hash)
    if protocol_path.exists():
        if protocol_path.read_text(encoding="utf-8") != protocol:
            raise ValueError("registered confirmation protocol mismatch")
    else:
        _write_text(protocol_path, protocol)
    records = [{**asdict(batch), "status": "pending", "frozen_config_sha256": frozen_config_sha256(),
                "experiment_source_manifest_sha256": source_hash, "attempt_history": []} for batch in plan]
    _write_json(registry_path, records)
    return records


@contextmanager
def _run_lock(outdir):
    # OS-owned locks release on process exit, so a crash does not block resume.
    # Keep the lock file: unlinking it can let another process lock a new inode.
    with (outdir / ".confirmation.lock").open("a+b") as handle:
        if handle.seek(0, 2) == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        if os.name == "nt":
            import msvcrt
            lock = lambda: msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            unlock = lambda: msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            lock = lambda: fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            unlock = lambda: fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        try:
            lock()
        except OSError as exc:
            raise ValueError("confirmation output is locked by another runner") from exc
        try:
            yield
        finally:
            handle.seek(0)
            unlock()


def run_confirmation(outdir, base_seed=DEFAULT_BASE_SEED):
    """Run fixed batches; stop on implementation failure, retaining every attempt."""
    outdir = Path(outdir)
    build_confirmation_plan(base_seed)
    outdir.mkdir(parents=True, exist_ok=True)
    with _run_lock(outdir):
        return _run_confirmation(outdir, base_seed)


def _run_confirmation(outdir, base_seed):
    outdir = Path(outdir)
    plan = build_confirmation_plan(base_seed)
    outdir.mkdir(parents=True, exist_ok=True)
    provenance = source_manifest()
    source_hash = provenance["source_manifest_sha256"]
    registry_path = outdir / "batch_registry.json"
    if registry_path.exists():
        records = _load_registered(outdir, plan, base_seed, allow_incomplete=True, current_source=source_hash)
    else:
        records = _initialize_registration(outdir, plan, base_seed, provenance)
    by_id = {record["batch_id"]: record for record in records}
    # Validate every completed batch before starting any new work.
    for batch in plan:
        record = by_id[batch.batch_id]
        if record["status"] == "completed":
            _validate_batch(batch, outdir / "batches" / batch.batch_id, record)
    for batch in plan:
        record = by_id[batch.batch_id]
        if record["status"] == "completed":
            continue
        batch_dir = outdir / "batches" / batch.batch_id
        if record["status"] in {"failed", "running"} or batch_dir.exists():
            history = record.setdefault("attempt_history", [])
            attempt_dir = outdir / "attempts" / batch.batch_id / f"attempt-{len(history) + 1:03d}"
            attempt_dir.parent.mkdir(parents=True, exist_ok=True)
            snapshot = {key: value for key, value in record.items() if key != "attempt_history"}
            snapshot["artifact_directory"] = attempt_dir.relative_to(outdir).as_posix()
            if batch_dir.exists():
                batch_dir.rename(attempt_dir)
            else:
                # A crash after rename/mkdir can leave this archive ahead of
                # the registry. With no source directory, finish that same
                # archival transaction; keep its files and append its snapshot.
                attempt_dir.mkdir(exist_ok=True)
            history.append(snapshot)
        record.pop("error", None)
        record.pop("artifact_sha256", None)
        record.pop("finished_at", None)
        record.update(status="running", started_at=datetime.now(timezone.utc).isoformat())
        _write_json(registry_path, records)
        try:
            run_p42a(batch_dir, _config(batch))
            if source_manifest()["source_manifest_sha256"] != source_hash:
                raise ValueError("generating source changed during confirmation run")
            metadata_path = batch_dir / "metadata.json"
            metadata = _json(metadata_path)
            metadata.update(experiment_source_manifest_sha256=source_hash,
                            frozen_config_sha256=frozen_config_sha256(),
                            git_revision=provenance["git_revision"], git_dirty=provenance["git_dirty"])
            _write_json(metadata_path, metadata)
            record["artifact_sha256"] = _hash_artifacts(batch_dir)
            _validate_batch(batch, batch_dir, record)
            record.update(status="completed", finished_at=datetime.now(timezone.utc).isoformat())
        except BaseException as exc:
            record.update(status="failed", error=f"{type(exc).__name__}: {exc}",
                          finished_at=datetime.now(timezone.utc).isoformat())
            _write_json(registry_path, records)
            if not isinstance(exc, Exception):
                raise
            return records
        _write_json(registry_path, records)
    _aggregate_confirmation_unlocked(outdir, base_seed)
    return records


def _load_baselines(outdir, records, frames):
    """Accept the existing baseline script's output only as a complete paired set."""
    directory = outdir / "baselines"
    if not directory.exists():
        return (), None
    metadata = _json(directory / "metadata.json")
    source_hashes = sorted({row["experiment_source_manifest_sha256"] for row in records})
    if metadata.get("source_experiment_manifest_hashes") != source_hashes:
        raise ValueError("baseline generating source does not match confirmation")
    # grid is historical provenance, not identity: archives must remain portable.
    # Identity is checked below against all registered batches, seeds and shots,
    # the generating source hashes and the saved QEC_SB/NetworkX predictions.
    if metadata.get("source_batch_count") != len(records):
        raise ValueError("baseline batch identity mismatch")
    if metadata.get("bp_iterations") != 12 or metadata.get("bp_damping") != .4:
        raise ValueError("baseline fixed BP config mismatch")
    baseline_source = metadata.get("source_manifest_sha256")
    if not isinstance(baseline_source, str) or len(baseline_source) != 64:
        raise ValueError("baseline source hash missing")
    try:
        data = pd.read_csv(directory / "shot_level.csv")
        summary = pd.read_csv(directory / "summary.csv")
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise ValueError("missing or invalid baseline CSV artifacts") from exc
    methods = ("QEC_SB", "NetworkX_MWPM", "PyMatching", "BP_Matching")
    _require_columns(data, ["batch_id", "distance", "p", "seed", "shot", "sample_sha256"] +
                     [f"{kind}_{method}" for method in methods for kind in ("correct", "pred")], "baseline")
    _require_columns(summary, ("batch_id", "method", "status", "shots", "sample_sha256"), "baseline summary")
    ids = {row["batch_id"] for row in records}
    if set(data.batch_id) != ids or set(summary.batch_id) != ids:
        raise ValueError("baseline missing or unexpected batches")
    if data.duplicated(["batch_id", "shot"]).any() or summary.duplicated(["batch_id", "method"]).any():
        raise ValueError("baseline duplicate batch/shot/method")
    for record in records:
        batch_id = record["batch_id"]
        source = frames[batch_id]
        frame = data[data.batch_id == batch_id].sort_values("shot").reset_index(drop=True)
        rows = summary[summary.batch_id == batch_id]
        if len(frame) != len(source) or frame.shot.tolist() != source.shot.tolist():
            raise ValueError("baseline shot coverage mismatch")
        for key in ("distance", "p", "seed"):
            if not frame[key].eq(record[key]).all():
                raise ValueError(f"baseline {key} mismatch")
        hashes = frame.sample_sha256.unique().tolist()
        if len(hashes) != 1 or not isinstance(hashes[0], str) or len(hashes[0]) != 64:
            raise ValueError("baseline sample hash mismatch")
        if Counter(rows.method) != Counter(methods) or not rows.status.eq("completed").all() or not rows.shots.eq(len(source)).all() or not rows.sample_sha256.eq(hashes[0]).all():
            raise ValueError("baseline summary completeness mismatch")
        for method in methods:
            correct, prediction = frame[f"correct_{method}"], frame[f"pred_{method}"]
            if not correct.isin((0, 1)).all() or not prediction.isin(range(4)).all():
                raise ValueError("baseline invalid binary correctness/prediction")
            if not np.array_equal(correct, (prediction == source.truth_class).astype(int)):
                raise ValueError("baseline correctness disagrees with paired truth")
        for method, original in (("QEC_SB", "P4_2_Bennett"), ("NetworkX_MWPM", "MWPM")):
            for kind in ("correct", "pred"):
                if not np.array_equal(frame[f"{kind}_{method}"], source[f"{kind}_{original}"]):
                    raise ValueError("baseline paired source prediction mismatch")
        for method in ("PyMatching", "BP_Matching"):
            source[f"correct_{method}"] = frame[f"correct_{method}"].to_numpy()
    return ("PyMatching", "BP_Matching"), baseline_source


def aggregate_confirmation(outdir, base_seed=DEFAULT_BASE_SEED):
    """Recompute solely from the twelve registered confirmation batch artifacts."""
    outdir = Path(outdir)
    with _run_lock(outdir):
        return _aggregate_confirmation_unlocked(outdir, base_seed)


def _aggregate_confirmation_unlocked(outdir, base_seed):
    plan = build_confirmation_plan(base_seed)
    records = _load_registered(outdir, plan, base_seed)
    by_id = {row["batch_id"]: row for row in records}
    frames = {}
    for batch in plan:
        frames[batch.batch_id] = _validate_batch(batch, outdir / "batches" / batch.batch_id, by_id[batch.batch_id])
    extra_methods, baseline_source = _load_baselines(outdir, records, frames)
    pooled_rows, paired_rows = [], []
    for index in range(6):
        batches = [batch for batch in plan if batch.family_index == index]
        batch = batches[0]
        shots = pd.concat([frames[item.batch_id] for item in batches], ignore_index=True)
        common = {"distance": batch.distance, "p": batch.p, "role": batch.role, "family_index": index,
                  "replicates": 2, "frozen_config_sha256": by_id[batch.batch_id]["frozen_config_sha256"],
                  "experiment_source_manifest_sha256": by_id[batch.batch_id]["experiment_source_manifest_sha256"]}
        for method in (*METHODS, *extra_methods):
            correct = shots[f"correct_{method}"].to_numpy(dtype=int)
            lo, hi = wilson_interval(int(correct.sum()), len(correct))
            pooled_rows.append({**common, "method": method, "shots": len(correct),
                                "successes": int(correct.sum()), "errors": int(len(correct) - correct.sum()),
                                "logical_success": correct.mean(), "logical_error": 1-correct.mean(),
                                "ci95_lo": lo, "ci95_hi": hi, "error_ci95_lo": 1-hi, "error_ci95_hi": 1-lo})
            if method in extra_methods:
                pooled_rows[-1]["baseline_source_manifest_sha256"] = baseline_source
        for baseline in ("MWPM", *extra_methods):
            stats = exact_paired_statistics(shots.correct_P4_2_Bennett, shots[f"correct_{baseline}"])
            baseline_error = 1-stats["baseline_success"]
            paired_rows.append({**common, "method": "P4_2_Bennett", "baseline": baseline, **stats,
                                "relative_logical_error_reduction": stats["paired_difference"]/baseline_error if baseline_error else float("nan"),
                                "correction_family": "six_predeclared_final_vs_MWPM" if baseline == "MWPM" else "supplementary_unadjusted",
                                "correction_family_size": 6 if baseline == "MWPM" else 0,
                                "holm_p_value": float("nan"), "bonferroni_p_value": float("nan")})
            if baseline in extra_methods:
                paired_rows[-1]["baseline_source_manifest_sha256"] = baseline_source
    family = [row for row in paired_rows if row["baseline"] == "MWPM"]
    p_values = [row["exact_p_value"] for row in family]
    for row, holm, bonferroni in zip(family, adjust_p_values(p_values, "holm"), adjust_p_values(p_values, "bonferroni")):
        row.update(holm_p_value=holm, bonferroni_p_value=bonferroni)
    pooled, paired = pd.DataFrame(pooled_rows), pd.DataFrame(paired_rows)
    _publish_results(outdir, pooled, paired)
    return pooled, paired


def _generation_id(digests):
    return hashlib.sha256(json.dumps(digests, sort_keys=True).encode("utf-8")).hexdigest()


def _publish_results(outdir, pooled, paired):
    """Stage both CSVs, then atomically publish their immutable generation pointer."""
    staged, contents, digests = {}, {}, {}
    for name, frame in (("pooled_results.csv", pooled), ("paired_results.csv", paired)):
        staged[name] = outdir / (name + ".tmp")
        frame.to_csv(staged[name], index=False)
        contents[name] = staged[name].read_bytes()
        digests[name] = hashlib.sha256(contents[name]).hexdigest()
    generation = _generation_id(digests)
    directory = outdir / "aggregates" / generation
    directory.mkdir(parents=True, exist_ok=True)
    for name, content in contents.items():
        path = directory / name
        if path.exists():
            if path.read_bytes() != content:
                raise ValueError("immutable aggregate generation artifact mismatch")
        else:
            _write_bytes(path, content)
    manifest = {
        "generation": generation,
        "files": {name: {"path": (directory / name).relative_to(outdir).as_posix(), "sha256": digests[name]}
                  for name in contents},
        "batch_registry_sha256": hashlib.sha256((outdir / "batch_registry.json").read_bytes()).hexdigest(),
    }
    manifest_path = outdir / "results_manifest.json"
    previous = {name: (outdir / name).read_bytes() if (outdir / name).exists() else None for name in contents}
    try:
        for name, temporary in staged.items():
            temporary.replace(outdir / name)
        _write_json(manifest_path, manifest)
    except BaseException:
        # The immutable previous pair remains available even after a hard kill.
        # On a handled write/replace failure, also restore compatibility copies.
        if not manifest_path.exists() or _json(manifest_path) != manifest:
            for name, content in previous.items():
                if content is None:
                    (outdir / name).unlink(missing_ok=True)
                else:
                    _write_bytes(outdir / name, content)
        raise


def load_confirmation_results(outdir):
    """Read one committed CSV pair, even if compatibility copies were interrupted."""
    outdir = Path(outdir)
    manifest = _json(outdir / "results_manifest.json")
    names = ("pooled_results.csv", "paired_results.csv")
    if set(manifest.get("files", {})) != set(names):
        raise ValueError("aggregate generation must contain both result files")
    digests = {name: manifest["files"][name]["sha256"] for name in names}
    generation = _generation_id(digests)
    if manifest.get("generation") != generation:
        raise ValueError("aggregate generation identity mismatch")
    if manifest.get("batch_registry_sha256") != hashlib.sha256((outdir / "batch_registry.json").read_bytes()).hexdigest():
        raise ValueError("aggregate generation registry hash mismatch")
    results = []
    for name in names:
        entry = manifest["files"][name]
        expected = Path("aggregates") / generation / name
        if entry.get("path") != expected.as_posix():
            raise ValueError("aggregate generation artifact path mismatch")
        path = outdir / expected
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError("aggregate generation artifact hash mismatch")
        results.append(pd.read_csv(path))
    return tuple(results)
