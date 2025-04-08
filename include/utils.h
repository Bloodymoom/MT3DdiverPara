#ifndef UTILS_H
#define UTILS_H
#include <stdio.h>
#include "fmodel.h"
#include <sys/types.h>
#include <sys/stat.h>

// 
typedef struct MonitorContext {
    PetscReal minRelativeResidual;
    PetscReal rtol;
    PetscReal bnorm;
} MonitorContext;

// Establish conversion relationship between edges and nodes
void edgeTonode(Fmodel *fmodel);
FILE *readModelFile(char *filename, double freq, int isloadB, Fmodel *fmodel);
double getRadian(int degree);
// Read boundary condition values
void fDirBdaries(v1AsEm *v1Ae, v2AsEm *v2Ae, double freq, int polarization, Fmodel *fm);
// Monitor and print iteration solving convergence status for Source A and Source B
PetscErrorCode MyKSPMonitor_xy(KSP ksp, PetscInt n, PetscReal rnorm, void *ctx);
PetscErrorCode MyKSPMonitor_yx(KSP ksp, PetscInt n, PetscReal rnorm, void *ctx);
// Obtain apparent resistivity and phase
void rhoAndpha(v1AsEm *v1Ae, v2AsEm *v2Ae, Fmodel *fm, double freq, MPI_Comm curComm);

#endif