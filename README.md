# MT3DdiverPara
Parallel 3-D MT forward modeling with divergence correction

# Requirements
mpich 3.4.2
PETSc 3.22.0
superlu_dist 8.0.0
parmetis 4.0.3
metis 5.1.0
lapack 3.8.0
openblas

# Make
$ make

# run code
$ mpirun -np 4 ./main -ksp_rtol 1e-10 -max_its 50
