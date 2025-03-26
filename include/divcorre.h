#ifndef DIVCORRE_H
#define DIVCORRE_H
#include "fmodel.h"
#include <petscmat.h>

typedef struct Divcorre 
{
    int *MEs, *MEv;
    double *df1, *df2, *df3,
           *df4, *df5, *df6,
           *df7, *df8, *df9;
    double *tf1, *tf2, *tf3;
    Mat Dv, tv;
}Divcorre;

void getMEs_v(Fmodel *fm, Divcorre *div);
void init_Dfactor(Divcorre *div);
void init_tfactor(Divcorre *div);
void freeDiv(Divcorre *div);
void freedf_tf(Divcorre *div);
void divFE(Fmodel *fm, Divcorre *div, MPI_Comm curComm);

#endif