"""Console and helper scripts shipped with schubmult."""

import os

# The CLI never does BLAS work; spinning up numpy's thread pool at import costs ~80 ms wall.
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
