#ifndef V1ASSEMBLY_H
#define V1ASSEMBLY_H
#include <petscmat.h>
#include <petscksp.h>

typedef struct v1AsEm
{
    PetscScalar *Exx1, *Eyy1, *Ezz1;
    PetscScalar *Hxx1, *Hyy1, *Hzz1;
    PetscScalar *ExA, *EyA, *EzA;
    // PetscScalar *HxA, *HyA;
    // double *HzA;

    // double rho_a1, rho_a2;

    Mat v1;
} v1AsEm;

void v1_Solve(v0AsEm *v0Ae, v1AsEm *v1Ae, Divcorre *div, Fmodel *fm, double freq, MPI_Comm curComm);
void freev1Ae(v1AsEm *v1Ae, MPI_Comm curComm);

#endif