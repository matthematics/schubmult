def test_star_import_works():
    namespace = {}
    exec("from schubmult import *", namespace, namespace)
    assert "Sx" in namespace
    assert "DSx" in namespace
    assert "Permutation" in namespace
