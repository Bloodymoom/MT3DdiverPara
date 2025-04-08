#ifndef V0ASSEMBLY_H
#define V0ASSEMBLY_H
#include <petscmat.h>

// Global impedance matrix structure, storing matrix data
typedef struct v0AsEm 
{
    // K11
    double *PxyTyx, *PxzTyx, *PxyTzx, *PxzTzx;
    // K21
    double *PyxTyx, *PyzTyx, *PyxTzx, *PyzTzx;
    // K31
    double *PzxTyx, *PzyTyx, *PzxTzx, *PzyTzx;
    // K12
    double *PxyTxy, *PxzTxy, *PxyTzy, *PxzTzy;
    // K22
    double *PyxTxy, *PyzTxy, *PyxTzy, *PyzTzy;
    // K32
    double *PzxTxy, *PzyTxy, *PzxTzy, *PzyTzy;
    // K13
    double *PxyTxz, *PxzTxz, *PxyTyz, *PxzTyz;
    // K23
    double *PyxTxz, *PyzTxz, *PyxTyz, *PyzTyz;
    // K33
    double *PzxTxz, *PzyTxz, *PzxTyz, *PzyTyz;

    double *k2exx, *k2exy, *k2exz;
    double *k2eyx, *k2eyy, *k2eyz;
    double *k2ezx, *k2ezy, *k2ezz;

    Mat v0;
}v0AsEm;

void getP_xyz(v0AsEm *v0Ae);
// Initialize KE matrix
void initK2eMat(v0AsEm *v0Ae);
// Finite element computation to obtain the global impedance matrix
void v0_FE(double freq, Fmodel *fm, v0AsEm *v0Ae, MPI_Comm curComm);
// Free data
void freev0Ae(v0AsEm *v0Ae);
void freeK1K2(v0AsEm *v0Ae);

#endif