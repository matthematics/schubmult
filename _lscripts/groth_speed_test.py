from schubmult import *
from schubmult.symbolic import expand

if __name__ == "__main__":
    ring = DGx([]).ring
    #result = ring.chevalley((0,1),uncode([0,1]))
    result = DGx(uncode([0,1]))*DGx(uncode((0,1)))
    #ring.genset[2] 
    polyresult = ring.genset[2] * DGx(uncode([0,1])).as_polynomial()
    print(f"elem: {result}")
    print(f"poly: {(polyresult-result.as_polynomial()).simplify()}")
    