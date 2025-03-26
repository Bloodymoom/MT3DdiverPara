# MT3DdiverPara
Parallel 3-D MT forward modeling with divergence correction

# Requirements
mpich 3.4.2<br>
PETSc 3.22.0<br>
superlu_dist 8.0.0<br>
parmetis 4.0.3<br>
metis 5.1.0<br>
lapack 3.8.0<br>
openblas

# Make
$ make

# run code
$ mpirun -np 4 ./main -ksp_rtol 1e-10 -max_its 50
