![Overall parallel structure of MT3DdiverPara](./kuangjia.png)

# Description
**MT3DdiverPara** implemented using MPI 3.4.2 and PETSc 3.22.0 versions. The overall parallel structure of this algorithm is controlled by three key parameters: N_PROC, MPI_FEM, and MPI_FRE. It governs the total number of processes in the communication domain PETSC_COMM_WORLD, the number of sub-communication domains, and the number of frequencies per subdomain. MPI_FEM = FRE_NUM / MPI_FRE, where FRE_NUM represents the number of frequencies. As shown in the diagram above, 8 processes (N_PROC = 8) are divided into four subdomains (MPI_FEM = 4). Through MPI_Split, the function maps global process ranks to subdomain process ranks, with each subdomain containing two processes that solve two frequencies in parallel (MPI_FRE = 2). Each subdomain stores frequency-dependent data and mesh data including node count, edge count, and element count. Within each subdomain, model initialization and data preprocessing occupy only a small portion of computation time. These steps are executed serially to minimize communication overhead caused by frequent inter-process interactions. Using the PETSc library, the coefficient matrices and solution vectors of the equations are converted to MPI AIJ distributed storage format and stored across multiple processes within subdomains, enabling parallel matrix assembly. During the parallel solving phase, each subdomain directly calls PETSc's parallel iterative solver, using two processes to handle two frequency datasets in parallel. The master process of each subdomain collects solutions from other processes for post-processing and result output.

# Make
`make`

# Requirements
```
mpich 3.4.2
PETSc 3.22.0
superlu_dist 8.0.0
parmetis 4.0.3
metis 5.1.0
lapack 3.8.0
openblas 0.3.29
```

Before running the program, it is necessary to install all the aforementioned dependency libraries. The versions of these libraries should not be arbitrarily changed to avoid potential conflicts. The `PETSc` library must be installed **after** all other dependencies, so please pay close attention to the installation order. Once all dependencies are successfully installed, the program can be compiled using the `make` command. After compilation is complete, the program can be executed using the `mpirun` command to generate the results.

# run code
`mpirun -np 4 ./main -ksp_rtol 1e-10 -max_its 50`

# important parameters
```
np = 4           # number of process
ksp_rtol = 1e-10 # relative tolerance of iterative solver
max_its = 50     # maximum number of iterations for divergence
```

The model data parameters are stored in the `model_data` directory. The program reads the corresponding model parameters from the `model.dat`, while the `freqs.txt` contains the frequency values to be processed and the number of frequencies assigned to each sub-communication domain. These two files, `model.dat` and `freqs.txt`, respectively control the model parameters and the number of frequencies to be processed by the program. By modifying these files, one can control the input of model parameters and frequency data.

# License
MT3DdiverPara is licensed under the MIT License
