# MT3DdiverPara
Parallel 3-D MT forward modeling with divergence correction

# Make
$ make

# run code
$ mpirun -np 4 ./main -ksp_rtol 1e-10 -max_its 50
