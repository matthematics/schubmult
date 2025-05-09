from schubmult.symbolic import expand, sympify


# def sympify(val):
#     try:
#         return symengine.sympify(val)
#     except symengine.SympifyError:
#         # print(f"Bagels {val=}")
#         return sympify(val)


# def expand(val, **kwargs):
#     try:
#         return symengine.expand(val, **kwargs)
#     except Exception:
#         # print(f"Bogels {val=}")
#         return expand(val, **kwargs)
