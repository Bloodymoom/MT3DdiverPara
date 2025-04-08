#ifndef V2ASSEMBLY_H
#define V2ASSEMBLY_H
#include <petscmat.h>
#include <petscksp.h>

// Store global impedance matrix data and solution results for source B
typedef struct v2AsEm
{
    PetscScalar *Exx2, *Eyy2, *Ezz2;
    PetscScalar *Hxx2, *Hyy2, *Hzz2;
    PetscScalar *ExB, *EyB, *EzB;

    Mat v2;
} v2AsEm;

// Parallel Solve for Source B
void v2_Solve(v0AsEm *v0Ae, v2AsEm *v2Ae, Divcorre *div, Fmodel *fm, double freq, MPI_Comm curComm);
void freev2Ae(v2AsEm *v2Ae, MPI_Comm curComm);

#endif