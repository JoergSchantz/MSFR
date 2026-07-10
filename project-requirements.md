# project requirements
- functions and modules must be optimized for long EM algorithm calculations
    - standard `niter = 50000`
- number of studies `S` is assumed to be small 
- number of observations per study `n_s` are assumed to be large
    - with small `n_s = 500` as a reference for a lower bound 
- number of function calls must to be kept at a minimum
- redundancy must be avoided
- direct BLAS calls are fastest for runtime
- if you add new files to `R/`, update `file-structure.md` accordingly