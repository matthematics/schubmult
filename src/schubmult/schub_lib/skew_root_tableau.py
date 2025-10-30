from schubmult.schub_lib.perm_lib import Permutation
from schubmult.schub_lib.rc_graph import RCGraph

from .root_tableau import RootTableau


class SkewRootTableau(RootTableau):

    def from_rc_graph(cls, rc: RCGraph, subword_indexes=()):
        rc_hw, raise_seq = rc.to_highest_weight()
        outer_shape = rc_hw.weight_tableau.shape
        base_word = rc.perm_word
        rc_skew = RCGraph([]).resize(len(rc))
        for i in range(len(rc)):
            if i in subword_indexes:
                rc_skew = rc_skew.toggle_ref_at(*rc.left_to_right_inversion_coords(i))
        skew_tableau = RootTableau.from_rc_graph(rc_skew)
        # skew plactic, reverse rectify to outer shape
        while skew_tableau.weight_tableau.shape != outer_shape:
            skew_tableau = skew_tableau.up_jdt_slide