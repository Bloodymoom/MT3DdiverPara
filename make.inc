#
#  The name of the libraries to be created/linked to
#

BLAS_DIR     =/thfs1/home/gongchy/whd/installation/openblas
LAPACK_DIR   =/thfs1/software/lapack/3.8.0-gcc9.3.0
PARMETIS_DIR =/thfs1/software/parmetis/4.0.3-gcc9.3.0-mpi-x
METIS_DIR    =/thfs1/software/metis/5.1.0-gcc9.3.0
SUPERLU_DIST =/thfs1/home/gongchy/whd/installation/superlu_dist
PETSC_DIR    =/thfs1/home/gongchy/whd/installation/petsc-3.22.0

LIBS    = -lpetsc
LIBS   += -lopenblas
LIBS   += -llapack -lsuperlu_dist
LIBS   += -lparmetis -lmetis
LIBS   += -lgfortran -ldl -lstdc++ -lX11 -lm

###############################################################

LDFLAGS += -L${BLAS_DIR}/install/lib
LDFLAGS += -L${LAPACK_DIR}/lib
LDFLAGS += -L${METIS_DIR}/lib
LDFLAGS += -L${PARMETIS_DIR}/lib
LDFLAGS += -L${PETSC_DIR}/arch-linux-c-opt/lib -L${PETSC_DIR}/install/lib
LDFLAGS += -L${SUPERLU_DIST}/install/lib

###############################################################

INCLUDES += -I${METIS_DIR}/include
INCLUDES += -I${PARMETIS_DIR}/include
INCLUDES += -I${PETSC_DIR}/include -I${PETSC_DIR}/arch-linux-c-opt/include
INCLUDES += -I${BLAS_DIR}/install/include
INCLUDES += -I${SUPERLU_DIST}/install/include

CC           = /thfs1/software/mpich/mpi-n-gcc9.3.0/bin/mpif90
CFLAGS 	     = -std=c99 -fopenmp -DNDEBUG -fPIC
##-fPIC
CXX          = /thfs1/software/mpich/mpi-n-gcc9.3.0/bin/mpicxx
CXXFLAGS     = -std=c++11 -fopenmp -DNDEBUG -fPIC
