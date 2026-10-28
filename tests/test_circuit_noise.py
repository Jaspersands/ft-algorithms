from collections import Counter

import pytest
import stabilizer_qec as sq
import stim

from ftalgo.calib import circuit_noise as cn


def _errors(dem: stim.DetectorErrorModel) -> Counter:
    out = Counter()
    for inst in dem.flattened():
        if inst.type == "error":
            out[(round(inst.args_copy()[0], 14), " ".join(sorted(str(t) for t in inst.targets_copy())))] += 1
    return out


@pytest.mark.parametrize("basis", ["z", "x"])
@pytest.mark.parametrize("d,rounds", [(3, 3), (5, 6)])
def test_sd6_reproduces_stabilizer_qec(basis, d, rounds):
    # The noise rules, applied to the noiseless circuit, give stabilizer-qec's own error model.
    ref = stim.Circuit(str(sq.memory_circuit(distance=d, rounds=rounds, p=0.001, basis=basis)))
    mine = cn.memory(d, rounds, basis, "sd6", 0.001)
    assert _errors(mine.detector_error_model()) == _errors(ref.detector_error_model())


@pytest.mark.parametrize("basis", ["z", "x"])
@pytest.mark.parametrize("d", [3, 5])
def test_si1000_cz_memory_is_deterministic_with_distance_d(basis, d):
    c = cn.memory(d, 2 * d, basis, "si1000", 0.001)
    gates = {i.name for i in c.flattened()} - {"DETECTOR", "OBSERVABLE_INCLUDE", "QUBIT_COORDS", "TICK"}
    assert gates == {"H", "CZ", "R", "M", "DEPOLARIZE1", "DEPOLARIZE2", "X_ERROR"}
    c.detector_error_model(decompose_errors=True)  # raises on a non-deterministic detector
    assert len(c.shortest_graphlike_error()) == d


def test_si1000_rates():
    p = 0.002
    c = cn.memory(3, 3, "z", "si1000", p)
    seen = Counter()
    for inst in c.flattened():
        if inst.name in ("X_ERROR", "DEPOLARIZE1", "DEPOLARIZE2"):
            seen[(inst.name, round(inst.gate_args_copy()[0] / p, 9))] += 1
    # Table 3 of Gidney, Newman, Brooks & Jones: CZ p; single-qubit and idle p/10; reset 2p;
    # measurement 5p flip then p depolarizing; qubits idle during measurement or reset 2p.
    assert set(seen) == {("DEPOLARIZE2", 1), ("DEPOLARIZE1", 0.1), ("X_ERROR", 2), ("X_ERROR", 5), ("DEPOLARIZE1", 1), ("DEPOLARIZE1", 2)}


def test_uniform_is_stims_generated_circuit():
    c = cn.memory(3, 3, "x", "uniform", 0.001)
    ref = stim.Circuit.generated("surface_code:rotated_memory_x", distance=3, rounds=3, after_clifford_depolarization=0.001,
                                 before_round_data_depolarization=0.001, before_measure_flip_probability=0.001,
                                 after_reset_flip_probability=0.001)
    assert c == ref
