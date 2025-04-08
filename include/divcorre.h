#ifndef DIVCORRE_H
#define DIVCORRE_H
#include "fmodel.h"
#include <petscmat.h>

// Divergence correction struct, stores data used for divergence correction
typedef struct Divcorre 
{
    int *MEs, *MEv;
    double *df1, *df2, *df3,
           *df4, *df5, *df6,
           *df7, *df8, *df9;
    double *tf1, *tf2, *tf3;
    Mat Dv, tv;
}Divcorre;

// Generate node and edge indices
void getMEs_v(Fmodel *fm, Divcorre *div);
// Initialize derivative coefficient values for differential equations
void init_Dfactor(Divcorre *div);
void init_tfactor(Divcorre *div);
// Free data
void freeDiv(Divcorre *div);
void freedf_tf(Divcorre *div);
// 
void divFE(Fmodel *fm, Divcorre *div, MPI_Comm curComm);

#endif