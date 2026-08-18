"""Trace crossings through Algorithm 3 using the rectification rule.

User's tracing rule, located inside Algorithm 2 (_pieri_rectify in the source):
when an insertion negates a root in a lower row i, that crossing is dumped and a
new crossing is generically inserted in the SAME row i.  The dumped crossing
CORRESPONDS to the newly inserted one.

This gives a labelling of crossings that survives the whole of Algorithm 3:
  - crossings never touched keep their label;
  - the chain crossings deleted in step 7 of Algorithm 3 are dumped and
    reinserted by the Pieri insertion, carrying their labels;
  - a rectification transfers the label of the negated crossing to the crossing
    that replaces it in the same row.

We reimplement the algorithm with labels attached, CHECK that the untraced
output agrees with RCGraph.zero_out_last_row, and then ask what the induced
bijection  crossings(B) -> crossings(Z(B))  preserves: row, letter, or root.
"""

import sys
from collections import Counter
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph, _is_row_root


class Traced:
    """An RCGraph together with a label on each crossing."""

    def __init__(self, rc, labels):
        self.rc = rc
        self.labels = labels

    def add(self, r, c, label):
        rc = self.rc.toggle_ref_at(r, c)
        lab = dict(self.labels)
        lab[(r, c)] = label
        return Traced(rc, lab)

    def remove(self, r, c):
        rc = self.rc.toggle_ref_at(r, c)
        lab = dict(self.labels)
        label = lab.pop((r, c))
        return Traced(rc, lab), label


def insert_row(T, row, descent, dict_by_a, dict_by_b, label_queue):
    """Faithful copy of _pieri_insert_row for backwards=True, left=False."""
    i = 0
    num_done = 0
    num_times = len(label_queue)
    guard = 0
    while num_done < num_times:
        guard += 1
        if guard > 5000:
            raise RuntimeError("insert_row runaway")
        i += 1
        flag = False
        if not T.rc.has_element(row, i):
            a, b = T.rc.right_root_at(row, i)
            if a < b:
                if _is_row_root(descent, (a, b)) and b not in dict_by_b:
                    T = T.add(row, i, label_queue[num_done])
                    dict_by_a.setdefault(a, set()).add(b)
                    dict_by_b[b] = a
                    flag = True
                elif a in dict_by_b and b > descent and b not in dict_by_b:
                    T = T.add(row, i, label_queue[num_done])
                    dict_by_a[dict_by_b[a]].add(b)
                    dict_by_b[b] = dict_by_b[a]
                    flag = True
            if flag:
                num_done += 1
            if row > 1 and not T.rc.is_valid:
                T = rectify(T, row - 1, descent, dict_by_a, dict_by_b)
    return T


def rectify(T, row_below, descent, dict_by_a, dict_by_b):
    """Faithful copy of _pieri_rectify for left=False, with label transfer."""
    if T.rc.is_valid:
        return T
    if row_below == 0:
        return T
    extra = descent if descent is not None else 0
    for j in range(1, T.rc.cols + extra + 5):
        if T.rc.is_valid:
            return T
        if T.rc.has_element(row_below, j):
            a, b = T.rc.right_root_at(row_below, j)
            top, bottom = max(a, b), min(a, b)
            if a < b:
                continue
            flag = False
            if bottom in dict_by_a and top in dict_by_a[bottom]:
                T, label = T.remove(row_below, j)
                dict_by_a[bottom].remove(top)
                if len(dict_by_a[bottom]) == 0:
                    del dict_by_a[bottom]
                del dict_by_b[top]
                flag = True
            elif bottom in dict_by_b and top in dict_by_b and dict_by_b[top] == dict_by_b[bottom]:
                T, label = T.remove(row_below, j)
                dict_by_a[dict_by_b[bottom]].remove(top)
                if len(dict_by_a[dict_by_b[top]]) == 0:
                    del dict_by_a[dict_by_b[top]]
                del dict_by_b[top]
                flag = True
            else:
                raise RuntimeError(f"could not rectify at {(row_below, j)} root {(a, b)}")
            if flag:
                # THE TRACING RULE: the dumped crossing's label is carried by the
                # crossing inserted next in this same row.
                T = insert_row(T, row_below, descent, dict_by_a, dict_by_b, [label])
    return rectify(T, row_below - 1, descent, dict_by_a, dict_by_b)


def traced_zero(B):
    """Algorithm 3 with labels. Returns (Z(B), label map on crossings)."""
    h = len(B)
    if len(B[-1]) != 0:
        raise ValueError("last row not empty")
    if len(B.perm.trimcode) <= h - 1:
        # trivial branch: delete the empty last row
        R = B.rowrange(0, h - 1)
        labels = {(i, L - i + 1): (i, L - i + 1) for i, row in enumerate(B, start=1) for L in row}
        return R, labels

    labels = {(i, L - i + 1): (i, L - i + 1) for i, row in enumerate(B, start=1) for L in row}
    T = Traced(B, labels)

    # step 7: locate and delete the chain, recording rows and labels
    descs = []
    diff_rows = []
    chain_labels = []
    interim = B
    Tcur = T
    while len(interim.perm.trimcode) > h - 1:
        d = len(interim.perm.trimcode)
        descs.append(d)
        nxt, row = interim.exchange_property(d, return_row=True)
        # identify the deleted crossing
        gone = set(_cross(interim)) - set(_cross(nxt))
        if len(gone) != 1:
            raise RuntimeError("exchange deleted != 1 crossing")
        (rr, cc), = gone
        chain_labels.append(Tcur.labels[(rr, cc)])
        lab = dict(Tcur.labels)
        del lab[(rr, cc)]
        Tcur = Traced(nxt, lab)
        diff_rows.append(row)
        interim = nxt

    # step 9: park the chain in the last row (these are placeholders, unlabelled)
    parked = RCGraph([*interim[:-1], tuple(sorted(descs, reverse=True))])
    Tcur = Traced(parked, Tcur.labels)

    # step 10: reinsert
    descent = len(B.perm.trimcode) - 1
    rows_grouping = {}
    for r, lab in zip(diff_rows, chain_labels):
        rows_grouping.setdefault(r, []).append(lab)
    dict_by_a, dict_by_b = {}, {}
    if diff_rows and max(diff_rows) > len(Tcur.rc):
        Tcur = Traced(Tcur.rc.extend(max(diff_rows) - len(Tcur.rc)), Tcur.labels)
    for row in sorted(rows_grouping.keys(), reverse=True):
        queue = rows_grouping[row]
        Tcur = insert_row(Tcur, row, descent, dict_by_a, dict_by_b, queue)
        if row > 1 and not Tcur.rc.is_valid:
            Tcur = rectify(Tcur, row - 1, descent, dict_by_a, dict_by_b)

    R = Tcur.rc.rowrange(0, h - 1)
    final = {pos: lab for pos, lab in Tcur.labels.items() if pos[0] <= h - 1}
    return R, final


def _cross(R):
    return {(i, L - i + 1) for i, row in enumerate(R, start=1) for L in row}


def maxd(p):
    return len(p.trimcode)


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def main(N=5):
    agree = disagree = 0
    total_lab = 0
    bij = 0
    keeps_row = keeps_letter = keeps_root = 0
    per_row = per_letter = per_root = 0
    per_tot = 0
    fails = []

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            h = maxd(w)
            if h < 2 or (w, h) in seen:
                continue
            seen.add((w, h))
            for B in rc0(w, h):
                try:
                    R, labels = traced_zero(B)
                except Exception as e:
                    fails.append((tuple(w[:n]), h, [tuple(r) for r in B], repr(e)))
                    continue
                ref = B.zero_out_last_row()
                if _cross(R) == _cross(ref):
                    agree += 1
                else:
                    disagree += 1
                    continue

                total_lab += 1
                # labels should be a bijection: label set == crossings of B
                srcs = set(labels.values())
                if srcs == _cross(B) and len(labels) == len(_cross(B)):
                    bij += 1
                else:
                    continue

                okr = okl = okrt = True
                for dst, src in labels.items():
                    per_tot += 1
                    if src[0] == dst[0]:
                        per_row += 1
                    else:
                        okr = False
                    if src[0] + src[1] == dst[0] + dst[1]:
                        per_letter += 1
                    else:
                        okl = False
                    if B.right_root_at(*src) == R.right_root_at(*dst):
                        per_root += 1
                    else:
                        okrt = False
                keeps_row += okr
                keeps_letter += okl
                keeps_root += okrt

    print(f"N={N}  (h = maxd(w), nontrivial branch)")
    print(f"  traced algorithm reproduces zero_out_last_row: {agree}, disagreements {disagree}, exceptions {len(fails)}")
    print(f"  labelling is a bijection onto crossings(B): {bij}/{total_lab}")
    print()
    print(f"  graphs where the tracing preserves ROW:    {keeps_row}/{bij}")
    print(f"  graphs where it preserves LETTER:          {keeps_letter}/{bij}")
    print(f"  graphs where it preserves ROOT:            {keeps_root}/{bij}")
    print(f"  per crossing:  row {per_row}/{per_tot}   letter {per_letter}/{per_tot}   root {per_root}/{per_tot}")
    for f in fails[:5]:
        print("   EXCEPTION:", f)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
