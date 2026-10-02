from schubmult import *

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    for perm in perms:
        graphs = RCGraph.all_rc_graphs(perm)
        forest_invariants = set({graph.forest_invariant: graph for graph in graphs}.values())
        forest_classes_of_words = {invariant: Permutation.forest_class_of(invariant.perm_word) for invariant in forest_invariants}
        forest_classes_of_graphs = {invariant: {graph for graph in graphs if graph.forest_invariant == invariant.forest_invariant} for invariant in forest_invariants}
        for invariant in forest_invariants:
            stack = [invariant]
            results = set()
            while stack:
                current = stack.pop()
                if current in results:
                    continue
                results.add(current)

                cms = current.all_chute_moves()
                #print(cms)
                for ((i, j), (i1, j1)) in cms:
                    if i1 != i + 1:
                        continue
                    new_result = current.toggle_ref_at(i,j).toggle_ref_at(i1, j1)
                    if new_result.forest_invariant == invariant.forest_invariant:
                        stack.append(new_result)
                inv_cms = current.all_inverse_chute_moves()
                #print(inv_cms)
                for ((i, j), (i1, j1)) in inv_cms:
                    if i1 != i - 1:
                        continue
                    new_result = current.toggle_ref_at(i,j).toggle_ref_at(i1, j1)
                    if new_result.forest_invariant == invariant.forest_invariant:
                        stack.append(new_result)
            assert results == forest_classes_of_graphs[invariant], f"Mismatch for invariant {invariant}, {results}, {forest_classes_of_graphs[invariant]}"
        print(f"Barf {perm}")
    print("SAFPJ")
                