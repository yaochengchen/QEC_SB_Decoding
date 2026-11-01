from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json


@dataclass(frozen=True)
class FrozenCPConfig:
    replicas: int = 8
    steps: int = 100
    dt: float = 0.45
    coupling_scale: float = 0.70
    warm_strength: float = 0.08
    random_restarts: int = 1
    warm_restarts: int = 1
    coupling_multipliers: tuple[float, ...] = (1.0,)


@dataclass(frozen=True)
class FrozenImportanceConfig:
    stage0_top_m: int = 1
    stage0_radius: int = 1
    stage0_ratio_threshold: float = 1.05
    top_centers_per_class: int = 2
    core_radius: int = 2
    proposal_rhos: tuple[float, ...] = (0.08, 0.20, 0.40, 0.50)
    uniform_mix: float = 0.05
    batch_size: int = 32
    min_samples: int = 64
    max_samples: int = 160
    z_score: float = 2.58
    prohibit_full_enumeration: bool = True


@dataclass(frozen=True)
class FrozenBennettConfig:
    lambdas: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)
    n_samples: int = 48
    n_chains: int = 4
    burnin: int = 48
    thin: int = 1
    pair_flip_prob: float = 0.20
    bootstrap_reps: int = 20
    refine_ratio_threshold: float = 1.20
    decision_log_margin: float = 0.03
    sign_confidence: float = 0.85
    confirm_on_disagreement: bool = True
    confirm_lambdas: tuple[float, ...] = (
        0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0
    )
    confirm_n_samples: int = 256
    confirm_n_chains: int = 8
    confirm_burnin: int = 128
    confirm_thin: int = 2
    confirm_bootstrap_reps: int = 30
    confirm_sign_confidence: float = 0.95
    confirm_log_margin: float = 0.03


@dataclass(frozen=True)
class FrozenP42AConfig:
    version: str = "P4.2a-frozen"
    distance: int | None = None
    cp: FrozenCPConfig = FrozenCPConfig()
    importance: FrozenImportanceConfig = FrozenImportanceConfig()
    bennett: FrozenBennettConfig = FrozenBennettConfig()

    def to_dict(self) -> dict:
        return asdict(self)


FROZEN_P42A = FrozenP42AConfig()


def frozen_config_sha256() -> str:
    payload = json.dumps(FROZEN_P42A.to_dict(), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def build_frozen_p42a_config(distance, p, shots, seed, solver_seed):
    from ..bridge_ratio import MultiStageBennettConfig
    from ..importance_mass import ImportanceMassConfig
    from ..p42a import P42AConfig
    from ..scalable import ScalableCPConfig

    return P42AConfig(
        distance=int(distance),
        p=float(p),
        shots=int(shots),
        seed=int(seed),
        solver_seed=int(solver_seed),
        cp=ScalableCPConfig(**asdict(FROZEN_P42A.cp)),
        importance=ImportanceMassConfig(**asdict(FROZEN_P42A.importance)),
        bennett=MultiStageBennettConfig(**asdict(FROZEN_P42A.bennett)),
    )


def assert_metadata_uses_frozen_config(config: dict) -> None:
    """Reject archived batches whose algorithm settings differ from the freeze."""
    expected = FROZEN_P42A.to_dict()
    for section in ("cp", "importance", "bennett"):
        actual_value = json.loads(json.dumps(config.get(section), sort_keys=True))
        expected_value = json.loads(json.dumps(expected[section], sort_keys=True))
        if actual_value != expected_value:
            raise ValueError(f"archived batch {section} configuration violates frozen P4.2a contract")

