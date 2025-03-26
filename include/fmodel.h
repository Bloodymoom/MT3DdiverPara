#ifndef FMODEL_H
#define FMODEL_H

#define pi_by_18 (double)(PI/18.0)
#define pi_by_9 (double)(PI/9.0)
#define pi_by_6 (double)(PI/6.0)
#define pi_by_3 (double)(PI/3.0)
#define pi_by_2 (double)(PI/2.0)

typedef struct Fmodel
{
    double *A_X, *B_Y, *C_Z;
    int NX, NY, NZ;
    int NE, NP, NL;
    int Nair, Nsea;
    double *rho, *alpha_S, *alpha_D, *alpha_L;
    double *ms, *ms_S, *ms_D, *ms_L;
    int *ME, *EtoN;
}Fmodel;

void init_Fmodel(Fmodel *fmodel, double freq);
void freeModel(Fmodel *fmodel);

#endif