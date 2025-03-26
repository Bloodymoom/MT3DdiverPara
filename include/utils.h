#ifndef UTILS_H
#define UTILS_H
#include <stdio.h>
#include "fmodel.h"
#include <sys/types.h>
#include <sys/stat.h>

typedef struct MonitorContext {
    PetscReal minRelativeResidual;
    PetscReal rtol;
    PetscReal bnorm;
} MonitorContext;

void edgeTonode(Fmodel *fmodel);
FILE *readModelFile(char *filename, double freq, int isloadB, Fmodel *fmodel);
double getRadian(int degree);
void fDirBdaries(v1AsEm *v1Ae, v2AsEm *v2Ae, double freq, int polarization, Fmodel *fm);
PetscErrorCode MyKSPMonitor_xy(KSP ksp, PetscInt n, PetscReal rnorm, void *ctx);
PetscErrorCode MyKSPMonitor_yx(KSP ksp, PetscInt n, PetscReal rnorm, void *ctx);
void rhoAndpha(v1AsEm *v1Ae, v2AsEm *v2Ae, Fmodel *fm, double freq, MPI_Comm curComm);

#endif