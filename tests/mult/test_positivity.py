"""Tests for ``schubmult.mult.positivity`` (manifestly positive double structure constants)."""

import itertools
import os
import random
import subprocess
import sys
import threading

import pytest
import sympy

from schubmult import Permutation
from schubmult.abc import beta, x, y, z
from schubmult.mult import positivity
from schubmult.mult.double import schubmult_double, schubmult_double_pair
from schubmult.mult.positivity import _presolve, _resolve_workers, compute_positive_rep, dualcoeff, forwardcoeff, groth_posify, groth_posify_many, posify
from schubmult.symbolic import expand, sympify_sympy
from schubmult.symbolic.poly.schub_poly import schubpoly
from schubmult.utils.schub_lib import will_formula_work

S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
sp = sympify_sympy

POSITIVE_FGL_VALS = [
    0,
    5,
    (z[1] - y[1]) / (1 + beta * y[1]),
    2 * beta * (z[1] - y[1]) ** 2 / ((1 + beta * y[1]) * (1 + beta * z[1])),
    (z[1] - y[2]) * (z[2] - y[1]) + 3 * (z[1] - y[1]) * (z[2] - y[2]),
    (z[2] - y[1]) / ((1 + beta * y[1]) * (1 + beta * z[1])),
    beta * (z[1] - y[1]) * (z[2] - y[1]) / ((1 + beta * y[1]) * (1 + beta * z[2])),
]


@pytest.mark.parametrize("val", POSITIVE_FGL_VALS)
def test_groth_posify_reconstructs_positive_fgl_coefficient(val):
    rep = groth_posify(expand(val), y, z, False)
    assert sp(rep - val).cancel() == 0
    assert rep == val


def test_groth_posify_rejects_nonpositive_coefficient():
    with pytest.raises(ValueError, match="no positive representation found"):
        groth_posify(y[1] - z[1], y, z, False)


@pytest.mark.parametrize("workers", [None, 1, 2, 4])
def test_groth_posify_many_matches_groth_posify(workers):
    vals = [expand(val) for val in POSITIVE_FGL_VALS]
    expected = [groth_posify(val, y, z, False) for val in vals]
    assert list(groth_posify_many(vals, y, z, max_workers=workers)) == expected


@pytest.mark.parametrize("workers", [1, 2])
@pytest.mark.parametrize("bad_at", [0, 1, 2])
def test_groth_posify_many_yields_results_up_to_the_failing_coefficient(workers, bad_at):
    good = (z[1] - y[1]) / (1 + beta * y[1])
    vals = [expand(good)] * 3
    vals[bad_at] = y[1] - z[1]
    results = groth_posify_many(vals, y, z, max_workers=workers)
    for _ in range(bad_at):
        assert next(results) == good
    with pytest.raises(ValueError, match="no positive representation found"):
        next(results)


def test_groth_posify_many_progress_runs_on_the_calling_thread_in_order():
    started = []
    vals = [expand(val) for val in POSITIVE_FGL_VALS]
    list(groth_posify_many(vals, y, z, max_workers=3, progress=lambda i: started.append((i, threading.current_thread()))))
    assert [i for i, _ in started] == list(range(len(vals)))
    assert {thread for _, thread in started} == {threading.current_thread()}


def test_groth_posify_many_checks_max_workers_eagerly():
    for bad in (0, -2):
        with pytest.raises(ValueError, match="max_workers"):
            groth_posify_many([], y, z, max_workers=bad)


def test_resolve_workers():
    assert _resolve_workers(3, 10) == 3
    # never more workers than coefficients, and always at least one
    assert _resolve_workers(8, 3) == 3
    assert _resolve_workers(4, 0) == 1
    assert _resolve_workers(None, 10) == 1


# The tests below drive the scheduler with stand-ins for its two seams: a coefficient is a label whose coroutine asks for
# one integer program named after it, and "bad" has no positive representation.  The solver process of a "hang"
# program is a sleeping subprocess with the working directory on its command line, like the temporary files of PuLP.
def _use_labelled_programs(monkeypatch):
    def steps(val, var2, var3):
        yield ([{val: 1}], {val: 1})
        if val == "bad":
            raise ValueError("no positive representation found")
        return val

    monkeypatch.setattr(positivity, "_posify_steps", steps)


def _hang_in_solver(seen, launched, stop, tmpdir):
    seen.update(tmpdir=tmpdir, solver=subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)", tmpdir]))
    launched.set()
    stop.wait(30)
    raise RuntimeError("integer program cancelled")


def _assert_cleaned_up(seen):
    # a solver that was not killed would still be sleeping and raise TimeoutExpired here
    assert seen["solver"].wait(timeout=10) != 0
    assert not os.path.exists(seen["tmpdir"])


def test_groth_posify_many_interrupt_stops_running_solvers(monkeypatch):
    seen, launched = {}, threading.Event()
    _use_labelled_programs(monkeypatch)
    monkeypatch.setattr(positivity, "_solve_milp", lambda svecs, target, msg, tmpdir, stop: _hang_in_solver(seen, launched, stop, tmpdir))

    def interrupt_while_solving(i):
        if i == 1:
            assert launched.wait(30)
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        list(groth_posify_many(["hang", "hang"], y, z, max_workers=2, progress=interrupt_while_solving))
    _assert_cleaned_up(seen)


def test_groth_posify_many_interrupt_of_an_inline_solve_stops_the_solver(monkeypatch):
    seen = {}
    _use_labelled_programs(monkeypatch)

    def solve(svecs, target, msg, tmpdir, stop):
        seen.update(tmpdir=tmpdir, solver=subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)", tmpdir]))
        raise KeyboardInterrupt

    monkeypatch.setattr(positivity, "_solve_milp", solve)
    with pytest.raises(KeyboardInterrupt):
        list(groth_posify_many(["hang"], y, z, max_workers=1))
    _assert_cleaned_up(seen)


@pytest.mark.parametrize("first", ["fast", "bad"])
def test_groth_posify_many_stops_later_solvers_when_it_ends_early(monkeypatch, first):
    seen, launched = {}, threading.Event()
    _use_labelled_programs(monkeypatch)

    def solve(svecs, target, msg, tmpdir, stop):
        if "hang" in target:
            _hang_in_solver(seen, launched, stop, tmpdir)
        # the first coefficient is only settled once the solver of the second one is running
        assert launched.wait(30)
        return [1]

    monkeypatch.setattr(positivity, "_solve_milp", solve)
    results = groth_posify_many([first, "hang"], y, z, max_workers=2)
    if first == "bad":
        with pytest.raises(ValueError, match="no positive representation found"):
            next(results)
    else:
        assert next(results) == "fast"
        results.close()
    _assert_cleaned_up(seen)


class _SetAfter:
    """Stands in for the stop event and turns set after a given number of polls."""

    def __init__(self, polls):
        self.polls = polls
        self.count = 0

    def is_set(self):
        self.count += 1
        return self.count > self.polls


def test_solve_milp_notices_a_stop_at_any_point_of_its_preparation():
    # nothing is fixed by the presolve, so every candidate takes part in the model
    cands = [{0: 1, 1: 1}] * 10 + [{0: 1, 1: -1}] * 10
    target = {0: 3, 1: 1}
    unset = _SetAfter(10**9)
    values = positivity._solve_milp(cands, target, stop=unset)
    assert sum(a * c[0] for a, c in zip(values, cands)) == 3
    # the presolve, the variables and the constraint matrix are each built candidate by candidate, and a worker thread that
    # keeps building next to an interrupted run makes the interrupt take seconds on the large programs
    assert unset.count >= 3 * len(cands)
    for polls in range(unset.count):
        with pytest.raises(RuntimeError, match="cancelled"):
            positivity._solve_milp(cands, target, stop=_SetAfter(polls))


def test_presolve_keeps_every_nonnegative_integer_solution():
    rng = random.Random(7)
    box = list(itertools.product(range(4), repeat=4))
    for _ in range(300):
        cands = [{r: c for r in range(3) if (c := rng.randint(-2, 2))} for _ in range(4)]
        # planted target: feasible by construction, so presolve must not reject it
        planted = rng.choice(box)
        target = {r: sum(a * cand.get(r, 0) for a, cand in zip(planted, cands)) for r in range(3)}
        pre = _presolve(cands, target)
        assert pre is not None
        fixed, active, resid = pre
        for sol in box:
            if all(sum(a * cand.get(r, 0) for a, cand in zip(sol, cands)) == target[r] for r in range(3)):
                assert all(sol[j] == fixed.get(j, 0) for j in range(4) if j not in active)
                assert all(sum(sol[j] * cands[j].get(r, 0) for j in active) == resid.get(r, 0) for r in range(3))
        # arbitrary target: a rejection means there is no solution at all
        other = {r: rng.randint(-3, 3) for r in range(3)}
        if _presolve(cands, other) is None:
            assert not any(all(sum(a * cand.get(r, 0) for a, cand in zip(sol, cands)) == other[r] for r in range(3)) for sol in box)


def test_compute_positive_rep_reconstructs_known_polynomial():
    val = 2 * (y[1] - z[1]) + 3 * (y[2] - z[1]) * (y[1] - z[2])
    rep = compute_positive_rep(val, y, z)
    assert expand(sp(rep) - sp(val)) == 0


def test_compute_positive_rep_integer_passthrough():
    assert compute_positive_rep(5, y, z) == 5


def test_posify_reconstructs_raw_value_for_s3():
    for u in S3:
        for v in S3:
            coeff_dict = schubmult_double({u: 1}, v, y, z)
            for w, val in coeff_dict.items():
                positive = posify(val, u, v, w, y, z)
                assert expand(sp(positive) - sp(val)) == 0, (u, v, w)


def test_dualcoeff_identity_u_matches_schubpoly():
    u = Permutation([])
    for v in S3:
        coeff_dict = schubmult_double_pair(u, v, y, z)
        for w, val in coeff_dict.items():
            vp = v * (~w)
            if vp.inv == v.inv - w.inv:
                expected = sp(schubpoly(vp, y, z))
                assert expand(sp(dualcoeff(u, v, w, y, z)) - expected) == 0, (v, w)


def test_forwardcoeff_matches_raw_value_when_formula_applies():
    found = 0
    for u in S3:
        for v in S3:
            if not will_formula_work(u, v):
                continue
            coeff_dict = schubmult_double_pair(u, v, y, z)
            for w, val in coeff_dict.items():
                found += 1
                assert expand(sp(forwardcoeff(u, v, w, y, z)) - sp(val)) == 0, (u, v, w)
    assert found > 0


def _manifestly_positive(expr):
    # a minus sign is only allowed directly in front of a z-variable term
    import re

    return re.search(r"-\s*(\d+\*)?[y(]", str(expr)) is None


def test_posify_optimize_controls_lp_fallback():
    u, v, w = Permutation([1, 2, 4, 3]), Permutation([1, 4, 3, 2]), Permutation([2, 3, 4, 1])
    val = schubmult_double({u: 1}, v, y, z)[w]
    assert posify(val, u, v, w, y, z, False, False) == val
    assert _manifestly_positive(posify(val, u, v, w, y, z))


def test_posify_one_dominates_recursion_uses_shifted_alphabet():
    # Regression: u one-dominates w, so posify recurses on (phi1(u), v3, phi1(w)) over the
    # shifted alphabet var2[1:] (a plain tuple). Previously the reduced coefficient was computed
    # in the unshifted alphabet and compute_positive_rep crashed on tuple.index.
    u, v, w = Permutation([2, 4, 3, 1]), Permutation([1, 4, 3, 2]), Permutation([3, 4, 5, 1, 2])
    val = schubmult_double({u: 1}, v, y, z)[w]
    positive = posify(val, u, v, w, y, z)
    assert expand(sp(positive) - sp(val)) == 0
    assert _manifestly_positive(positive)


def test_compute_positive_rep_accepts_plain_tuple_alphabets():
    val = (y[2] - z[1]) + 2 * (y[3] - z[2])
    rep = compute_positive_rep(val, tuple(y[1:5]), tuple(z[1:5]))
    assert expand(sp(rep) - sp(val)) == 0


def test_posify_reconstructs_raw_value_for_s4():
    S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]
    for u in S4:
        for v in S4:
            for w, val in schubmult_double({u: 1}, v, y, z).items():
                positive = posify(val, u, v, w, y, z)
                assert expand(sp(positive) - sp(val)) == 0, (u, v, w)
                assert _manifestly_positive(positive), (u, v, w, positive)
