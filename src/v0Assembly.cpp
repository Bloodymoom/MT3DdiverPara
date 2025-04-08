#include "../include/MT3D.h"

void mergeKmat(double *Ke, double *Kexx, double *Kexy, double *Kexz,
               double *Keyx, double *Keyy, double *Keyz,
               double *Kezx, double *Kezy, double *Kezz)
{
    for (int i = 0, row = 0; i < 4; i++, row++)
    {
        for (int j = 0; j < 4; j++)
        {
            Ke[i * 12 + j] = Kexx[row * 4 + j];
            Ke[i * 12 + (j + 4)] = Kexy[row * 4 + j];
            Ke[i * 12 + (j + 2 * 4)] = Kexz[row * 4 + j];
        }
    }
    for (int i = 4, row = 0; i < 8; i++, row++)
    {
        for (int j = 0; j < 4; j++)
        {
            Ke[i * 12 + j] = Keyx[row * 4 + j];
            Ke[i * 12 + (j + 4)] = Keyy[row * 4 + j];
            Ke[i * 12 + (j + 2 * 4)] = Keyz[row * 4 + j];
        }
    }
    for (int i = 8, row = 0; i < 12; i++, row++)
    {
        for (int j = 0; j < 4; j++)
        {
            Ke[i * 12 + j] = Kezx[row * 4 + j];
            Ke[i * 12 + (j + 4)] = Kezy[row * 4 + j];
            Ke[i * 12 + (j + 2 * 4)] = Kezz[row * 4 + j];
        }
    }
}

void getP_xyz(v0AsEm *v0Ae)
{
    // K123-1
    double tmp_PxyTyx[4 * 4] = {2, 1, -2, -1, 1, 2, -1, -2, -2, -1, 2, 1, -1, -2, 1, 2};
    double tmp_PxzTyx[4 * 4] = {-1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1};
    double tmp_PxyTzx[4 * 4] = {-1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1};
    double tmp_PxzTzx[4 * 4] = {2, -2, 1, -1, -2, 2, -1, 1, 1, -1, 2, -2, -1, 1, -2, 2};

    double tmp_PyxTyx[4 * 4] = {-1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1};
    double tmp_PyzTyx[4 * 4] = {1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1};
    double tmp_PyxTzx[4 * 4] = {1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1};
    double tmp_PyzTzx[4 * 4] = {-2, 2, -1, 1, -1, 1, -2, 2, 2, -2, 1, -1, 1, -1, 2, -2};

    double tmp_PzxTyx[4 * 4] = {1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1};
    double tmp_PzyTyx[4 * 4] = {-2, -1, 2, 1, 2, 1, -2, -1, -1, -2, 1, 2, 1, 2, -1, -2};
    double tmp_PzxTzx[4 * 4] = {-1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1};
    double tmp_PzyTzx[4 * 4] = {1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1};

    // K123-2
    double tmp_PxyTxy[4 * 4] = {-1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1};
    double tmp_PxzTxy[4 * 4] = {1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1};
    double tmp_PxyTzy[4 * 4] = {1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1};
    double tmp_PxzTzy[4 * 4] = {-2, -1, 2, 1, 2, 1, -2, -1, -1, -2, 1, 2, 1, 2, -1, -2};

    double tmp_PyxTxy[4 * 4] = {2, -2, 1, -1, -2, 2, -1, 1, 1, -1, 2, -2, -1, 1, -2, 2};
    double tmp_PyzTxy[4 * 4] = {-1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1};
    double tmp_PyxTzy[4 * 4] = {-1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1};
    double tmp_PyzTzy[4 * 4] = {2, 1, -2, -1, 1, 2, -1, -2, -2, -1, 2, 1, -1, -2, 1, 2};

    double tmp_PzxTxy[4 * 4] = {-2, 2, -1, 1, -1, 1, -2, 2, 2, -2, 1, -1, 1, -1, 2, -2};
    double tmp_PzyTxy[4 * 4] = {1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1};
    double tmp_PzxTzy[4 * 4] = {1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1};
    double tmp_PzyTzy[4 * 4] = {-1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1};

    // K123-3
    double tmp_PxyTxz[4 * 4] = {1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1};
    double tmp_PxzTxz[4 * 4] = {-1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1};
    double tmp_PxyTyz[4 * 4] = {-2, 2, -1, 1, -1, 1, -2, 2, 2, -2, 1, -1, 1, -1, 2, -2};
    double tmp_PxzTyz[4 * 4] = {1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1};

    double tmp_PyxTxz[4 * 4] = {-2, -1, 2, 1, 2, 1, -2, -1, -1, -2, 1, 2, 1, 2, -1, -2};
    double tmp_PyzTxz[4 * 4] = {1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1};
    double tmp_PyxTyz[4 * 4] = {1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1};
    double tmp_PyzTyz[4 * 4] = {-1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1};

    double tmp_PzxTxz[4 * 4] = {2, 1, -2, -1, 1, 2, -1, -2, -2, -1, 2, 1, -1, -2, 1, 2};
    double tmp_PzyTxz[4 * 4] = {-1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1};
    double tmp_PzxTyz[4 * 4] = {-1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1};
    double tmp_PzyTyz[4 * 4] = {2, -2, 1, -1, -2, 2, -1, 1, 1, -1, 2, -2, -1, 1, -2, 2};

    v0Ae->PxyTyx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PxzTyx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PxyTzx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PxzTzx = (double *)malloc(4 * 4 * sizeof(double));

    memcpy(v0Ae->PxyTyx, tmp_PxyTyx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PxzTyx, tmp_PxzTyx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PxyTzx, tmp_PxyTzx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PxzTzx, tmp_PxzTzx, 4 * 4 * sizeof(double));

    v0Ae->PyxTyx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PyzTyx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PyxTzx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PyzTzx = (double *)malloc(4 * 4 * sizeof(double));

    memcpy(v0Ae->PyxTyx, tmp_PyxTyx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PyzTyx, tmp_PyzTyx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PyxTzx, tmp_PyxTzx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PyzTzx, tmp_PyzTzx, 4 * 4 * sizeof(double));

    v0Ae->PzxTyx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PzyTyx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PzxTzx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PzyTzx = (double *)malloc(4 * 4 * sizeof(double));

    memcpy(v0Ae->PzxTyx, tmp_PzxTyx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PzyTyx, tmp_PzyTyx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PzxTzx, tmp_PzxTzx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PzyTzx, tmp_PzyTzx, 4 * 4 * sizeof(double));

    v0Ae->PxyTxy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PxzTxy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PxyTzy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PxzTzy = (double *)malloc(4 * 4 * sizeof(double));

    memcpy(v0Ae->PxyTxy, tmp_PxyTxy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PxzTxy, tmp_PxzTxy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PxyTzy, tmp_PxyTzy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PxzTzy, tmp_PxzTzy, 4 * 4 * sizeof(double));

    v0Ae->PyxTxy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PyzTxy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PyxTzy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PyzTzy = (double *)malloc(4 * 4 * sizeof(double));

    memcpy(v0Ae->PyxTxy, tmp_PyxTxy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PyzTxy, tmp_PyzTxy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PyxTzy, tmp_PyxTzy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PyzTzy, tmp_PyzTzy, 4 * 4 * sizeof(double));

    v0Ae->PzxTxy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PzyTxy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PzxTzy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PzyTzy = (double *)malloc(4 * 4 * sizeof(double));

    memcpy(v0Ae->PzxTxy, tmp_PzxTxy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PzyTxy, tmp_PzyTxy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PzxTzy, tmp_PzxTzy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PzyTzy, tmp_PzyTzy, 4 * 4 * sizeof(double));

    v0Ae->PxyTxz = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PxzTxz = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PxyTyz = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PxzTyz = (double *)malloc(4 * 4 * sizeof(double));

    memcpy(v0Ae->PxyTxz, tmp_PxyTxz, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PxzTxz, tmp_PxzTxz, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PxyTyz, tmp_PxyTyz, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PxzTyz, tmp_PxzTyz, 4 * 4 * sizeof(double));

    v0Ae->PyxTxz = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PyzTxz = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PyxTyz = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PyzTyz = (double *)malloc(4 * 4 * sizeof(double));

    memcpy(v0Ae->PyxTxz, tmp_PyxTxz, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PyzTxz, tmp_PyzTxz, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PyxTyz, tmp_PyxTyz, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PyzTyz, tmp_PyzTyz, 4 * 4 * sizeof(double));

    v0Ae->PzxTxz = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PzyTxz = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PzxTyz = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->PzyTyz = (double *)malloc(4 * 4 * sizeof(double));

    memcpy(v0Ae->PzxTxz, tmp_PzxTxz, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PzyTxz, tmp_PzyTxz, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PzxTyz, tmp_PzxTyz, 4 * 4 * sizeof(double));
    memcpy(v0Ae->PzyTyz, tmp_PzyTyz, 4 * 4 * sizeof(double));
}

void initK2eMat(v0AsEm *v0Ae)
{
    double tmp_K2exx[4 * 4] = {4, 2, 2, 1, 2, 4, 1, 2, 2, 1, 4, 2, 1, 2, 2, 4};
    double tmp_K2exy[4 * 4] = {2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2};
    double tmp_K2exz[4 * 4] = {2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2};

    double tmp_K2eyx[4 * 4] = {2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2};
    double tmp_K2eyy[4 * 4] = {4, 2, 2, 1, 2, 4, 1, 2, 2, 1, 4, 2, 1, 2, 2, 4};
    double tmp_K2eyz[4 * 4] = {2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2};

    double tmp_K2ezx[4 * 4] = {2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2};
    double tmp_K2ezy[4 * 4] = {2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2};
    double tmp_K2ezz[4 * 4] = {4, 2, 2, 1, 2, 4, 1, 2, 2, 1, 4, 2, 1, 2, 2, 4};

    v0Ae->k2exx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->k2exy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->k2exz = (double *)malloc(4 * 4 * sizeof(double));

    v0Ae->k2eyx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->k2eyy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->k2eyz = (double *)malloc(4 * 4 * sizeof(double));

    v0Ae->k2ezx = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->k2ezy = (double *)malloc(4 * 4 * sizeof(double));
    v0Ae->k2ezz = (double *)malloc(4 * 4 * sizeof(double));

    memcpy(v0Ae->k2exx, tmp_K2exx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->k2exy, tmp_K2exy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->k2exz, tmp_K2exz, 4 * 4 * sizeof(double));

    memcpy(v0Ae->k2eyx, tmp_K2eyx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->k2eyy, tmp_K2eyy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->k2eyz, tmp_K2eyz, 4 * 4 * sizeof(double));

    memcpy(v0Ae->k2ezx, tmp_K2ezx, 4 * 4 * sizeof(double));
    memcpy(v0Ae->k2ezy, tmp_K2ezy, 4 * 4 * sizeof(double));
    memcpy(v0Ae->k2ezz, tmp_K2ezz, 4 * 4 * sizeof(double));
}

// Finite element analysis to generate the global impedance matrix
void v0_FE(double freq, Fmodel *fm, v0AsEm *v0Ae, MPI_Comm curComm)
{
    int NX, NY, NE, NL;
    int curRank, curSize;
    PetscLogDouble st_time, ed_time, time;
    NX = fm->NX;
    NY = fm->NY;
    NE = fm->NE;
    NL = fm->NL;
    double w = 2 * PI * freq;
    double prefixVal = -1.0 * w * mu0;
    MPI_Comm_size(curComm, &curSize);
    MPI_Comm_rank(curComm, &curRank);

    MatCreate(curComm, &v0Ae->v0);
    MatSetSizes(v0Ae->v0, PETSC_DECIDE, PETSC_DECIDE, NL, NL);

    double *sig_prm, *RzfuAl_S, *RxfuAl_D, *RzfuAl_L;
    double *RzAl_S, *RxAl_D, *RzAl_L, *sig_tensor;
    double *ms_prm, *RzfuMs_S, *RxfuMs_D, *RzfuMs_L;
    double *RzMs_S, *RxMs_D, *RzMs_L, *ms_tensor;
    double *k1exx, *k1eyx, *k1ezx;
    double *k1exy, *k1eyy, *k1ezy;
    double *k1exz, *k1eyz, *k1ezz;
    double *K1e, *K2e;
    double *T1, *T2, *eye;
    lapack_int *ipiv_m;
    int a_xn, b_yn, c_zn;
    double a, b, c;
    int NJ, NK;
    int h, slice, st_idx, ed_idx, i, j, k;
    std::complex<double> val;
    sig_prm = (double *)malloc(3 * 3 * sizeof(double));
    RzfuAl_S = (double *)malloc(3 * 3 * sizeof(double));
    RxfuAl_D = (double *)malloc(3 * 3 * sizeof(double));
    RzfuAl_L = (double *)malloc(3 * 3 * sizeof(double));
    RzAl_S = (double *)malloc(3 * 3 * sizeof(double));
    RxAl_D = (double *)malloc(3 * 3 * sizeof(double));
    RzAl_L = (double *)malloc(3 * 3 * sizeof(double));
    sig_tensor = (double *)malloc(3 * 3 * sizeof(double));

    ms_prm = (double *)malloc(3 * 3 * sizeof(double));
    RzfuMs_S = (double *)malloc(3 * 3 * sizeof(double));
    RxfuMs_D = (double *)malloc(3 * 3 * sizeof(double));
    RzfuMs_L = (double *)malloc(3 * 3 * sizeof(double));
    RzMs_S = (double *)malloc(3 * 3 * sizeof(double));
    RxMs_D = (double *)malloc(3 * 3 * sizeof(double));
    RzMs_L = (double *)malloc(3 * 3 * sizeof(double));
    ms_tensor = (double *)malloc(3 * 3 * sizeof(double));

    k1exx = (double *)malloc(4 * 4 * sizeof(double));
    k1eyx = (double *)malloc(4 * 4 * sizeof(double));
    k1ezx = (double *)malloc(4 * 4 * sizeof(double));

    k1exy = (double *)malloc(4 * 4 * sizeof(double));
    k1eyy = (double *)malloc(4 * 4 * sizeof(double));
    k1ezy = (double *)malloc(4 * 4 * sizeof(double));

    k1exz = (double *)malloc(4 * 4 * sizeof(double));
    k1eyz = (double *)malloc(4 * 4 * sizeof(double));
    k1ezz = (double *)malloc(4 * 4 * sizeof(double));

    K1e = (double *)malloc(3 * 4 * 4 * 3 * sizeof(double));
    K2e = (double *)malloc(3 * 4 * 4 * 3 * sizeof(double));

    T1 = (double *)malloc(3 * 3 * sizeof(double));
    T2 = (double *)malloc(3 * 3 * sizeof(double));
    eye = (double *)malloc(3 * 3 * sizeof(double));
    ipiv_m = (lapack_int *)malloc(sizeof(lapack_int) * 3);
    memset(eye, 0, 3 * 3 * sizeof(double));
    eye[0 * 3 + 0] = 1.0;
    eye[1 * 3 + 1] = 1.0;
    eye[2 * 3 + 2] = 1.0;

    double sig_xx, sig_xy, sig_xz;
    double sig_yx, sig_yy, sig_yz;
    double sig_zx, sig_zy, sig_zz;
    double mur_xx, mur_xy, mur_xz;
    double mur_yx, mur_yy, mur_yz;
    double mur_zx, mur_zy, mur_zz;
    double r_1, r_2, r_3;
    double a_S, a_D, a_L;
    double ms_1, ms_2, ms_3;
    double m_S, m_D, m_L;

    memset(sig_prm, 0, 3 * 3 * sizeof(double));
    memset(RzfuAl_S, 0, 3 * 3 * sizeof(double));
    memset(RxfuAl_D, 0, 3 * 3 * sizeof(double));
    memset(RzfuAl_L, 0, 3 * 3 * sizeof(double));
    memset(RzAl_S, 0, 3 * 3 * sizeof(double));
    memset(RxAl_D, 0, 3 * 3 * sizeof(double));
    memset(RzAl_L, 0, 3 * 3 * sizeof(double));
    memset(sig_tensor, 0, 3 * 3 * sizeof(double));

    memset(ms_prm, 0, 3 * 3 * sizeof(double));
    memset(RzfuMs_S, 0, 3 * 3 * sizeof(double));
    memset(RxfuMs_D, 0, 3 * 3 * sizeof(double));
    memset(RzfuMs_L, 0, 3 * 3 * sizeof(double));
    memset(RzMs_S, 0, 3 * 3 * sizeof(double));
    memset(RxMs_D, 0, 3 * 3 * sizeof(double));
    memset(RzMs_L, 0, 3 * 3 * sizeof(double));
    memset(ms_tensor, 0, 3 * 3 * sizeof(double));

    slice = floor((double)NE / curSize);
    st_idx = curRank * slice;
    ed_idx = (curRank != curSize - 1) ? (curRank + 1) * slice : NE;
    PetscTime(&st_time);
    for (h = st_idx; h < ed_idx; h++)
    {
        memset(sig_prm, 0, 3 * 3 * sizeof(double));
        memset(ms_prm, 0, 3 * 3 * sizeof(double));
        memset(T1, 0, 3 * 3 * sizeof(double));
        memset(T2, 0, 3 * 3 * sizeof(double));
        memset(k1exx, 0, 4 * 4 * sizeof(double));
        memset(k1exy, 0, 4 * 4 * sizeof(double));
        memset(k1exz, 0, 4 * 4 * sizeof(double));
        memset(k1eyx, 0, 4 * 4 * sizeof(double));
        memset(k1eyy, 0, 4 * 4 * sizeof(double));
        memset(k1eyz, 0, 4 * 4 * sizeof(double));
        memset(k1ezx, 0, 4 * 4 * sizeof(double));
        memset(k1ezy, 0, 4 * 4 * sizeof(double));
        memset(k1ezz, 0, 4 * 4 * sizeof(double));

        // Compute anisotropic conductivity and permeability tensors
        r_1 = fm->rho[h * 3 + 0];
        r_2 = fm->rho[h * 3 + 1];
        r_3 = fm->rho[h * 3 + 2];
        a_S = fm->alpha_S[h];
        a_D = fm->alpha_D[h];
        a_L = fm->alpha_L[h];

        sig_prm[0 * 3 + 0] = 1.0 / (r_1);
        sig_prm[1 * 3 + 1] = 1.0 / (r_2);
        sig_prm[2 * 3 + 2] = 1.0 / (r_3);

        RzfuAl_S[0 * 3 + 0] = cos(-a_S);
        RzfuAl_S[0 * 3 + 1] = sin(-a_S);
        RzfuAl_S[1 * 3 + 0] = -sin(-a_S);
        RzfuAl_S[1 * 3 + 1] = cos(-a_S);
        RzfuAl_S[2 * 3 + 2] = 1;

        RxfuAl_D[0 * 3 + 0] = 1;
        RxfuAl_D[1 * 3 + 1] = cos(-a_D);
        RxfuAl_D[1 * 3 + 2] = sin(-a_D);
        RxfuAl_D[2 * 3 + 1] = -sin(-a_D);
        RxfuAl_D[2 * 3 + 2] = cos(-a_D);

        RzfuAl_L[0 * 3 + 0] = cos(-a_L);
        RzfuAl_L[0 * 3 + 1] = sin(-a_L);
        RzfuAl_L[1 * 3 + 0] = -sin(-a_L);
        RzfuAl_L[1 * 3 + 1] = cos(-a_L);
        RzfuAl_L[2 * 3 + 2] = 1;

        RzAl_S[0 * 3 + 0] = cos(a_S);
        RzAl_S[0 * 3 + 1] = sin(a_S);
        RzAl_S[1 * 3 + 0] = -sin(a_S);
        RzAl_S[1 * 3 + 1] = cos(a_S);
        RzAl_S[2 * 3 + 2] = 1;

        RxAl_D[0 * 3 + 0] = 1;
        RxAl_D[1 * 3 + 1] = cos(a_D);
        RxAl_D[1 * 3 + 2] = sin(a_D);
        RxAl_D[2 * 3 + 1] = -sin(a_D);
        RxAl_D[2 * 3 + 2] = cos(a_D);

        RzAl_L[0 * 3 + 0] = cos(a_L);
        RzAl_L[0 * 3 + 1] = sin(a_L);
        RzAl_L[1 * 3 + 0] = -sin(a_L);
        RzAl_L[1 * 3 + 1] = cos(a_L);
        RzAl_L[2 * 3 + 2] = 1;

        ms_1 = fm->ms[h * 3 + 0];
        ms_2 = fm->ms[h * 3 + 1];
        ms_3 = fm->ms[h * 3 + 2];
        m_S = fm->ms_S[h];
        m_D = fm->ms_D[h];
        m_L = fm->ms_L[h];

        ms_prm[0 * 3 + 0] = ms_1;
        ms_prm[1 * 3 + 1] = ms_2;
        ms_prm[2 * 3 + 2] = ms_3;

        RzfuMs_S[0 * 3 + 0] = cos(-m_S);
        RzfuMs_S[0 * 3 + 1] = sin(-m_S);
        RzfuMs_S[1 * 3 + 0] = -sin(-m_S);
        RzfuMs_S[1 * 3 + 1] = cos(-m_S);
        RzfuMs_S[2 * 3 + 2] = 1;

        RxfuMs_D[0 * 3 + 0] = 1;
        RxfuMs_D[1 * 3 + 1] = cos(-m_D);
        RxfuMs_D[1 * 3 + 2] = sin(-m_D);
        RxfuMs_D[2 * 3 + 1] = -sin(-m_D);
        RxfuMs_D[2 * 3 + 2] = cos(-m_D);

        RzfuMs_L[0 * 3 + 0] = cos(-m_L);
        RzfuMs_L[0 * 3 + 1] = sin(-m_L);
        RzfuMs_L[1 * 3 + 0] = -sin(-m_L);
        RzfuMs_L[1 * 3 + 1] = cos(-m_L);
        RzfuMs_L[2 * 3 + 2] = 1;

        RzMs_S[0 * 3 + 0] = cos(m_S);
        RzMs_S[0 * 3 + 1] = sin(m_S);
        RzMs_S[1 * 3 + 0] = -sin(m_S);
        RzMs_S[1 * 3 + 1] = cos(m_S);
        RzMs_S[2 * 3 + 2] = 1;

        RxMs_D[0 * 3 + 0] = 1;
        RxMs_D[1 * 3 + 1] = cos(m_D);
        RxMs_D[1 * 3 + 2] = sin(m_D);
        RxMs_D[2 * 3 + 1] = -sin(m_D);
        RxMs_D[2 * 3 + 2] = cos(m_D);

        RzMs_L[0 * 3 + 0] = cos(m_L);
        RzMs_L[0 * 3 + 1] = sin(m_L);
        RzMs_L[1 * 3 + 0] = -sin(m_L);
        RzMs_L[1 * 3 + 1] = cos(m_L);
        RzMs_L[2 * 3 + 2] = 1;

        // sigma_tensor=RzfuAlpha_S*RxfuAlpha_D*RzfuAlpha_L*sigma_primary*RzAlpha_L*RxAlpha_D*RzAlpha_S
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    RzfuAl_S, 3, RxfuAl_D, 3, 0.0, T1, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    T1, 3, RzfuAl_L, 3, 0.0, T2, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    T2, 3, sig_prm, 3, 0.0, T1, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    T1, 3, RzAl_S, 3, 0.0, T2, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    T2, 3, RxAl_D, 3, 0.0, T1, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    T1, 3, RzAl_L, 3, 0.0, sig_tensor, 3);

        sig_xx = sig_tensor[0 * 3 + 0];
        sig_xy = sig_tensor[0 * 3 + 1];
        sig_xz = sig_tensor[0 * 3 + 2];
        sig_yx = sig_tensor[1 * 3 + 0];
        sig_yy = sig_tensor[1 * 3 + 1];
        sig_yz = sig_tensor[1 * 3 + 2];
        sig_zx = sig_tensor[2 * 3 + 0];
        sig_zy = sig_tensor[2 * 3 + 1];
        sig_zz = sig_tensor[2 * 3 + 2];

        // ms_tensor=RzfuMs_S*RxfuMs_D*RzfuMs_L*ms_primary*RzMs_L*RxMs_D*RzMs_S;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    RzfuMs_S, 3, RxfuMs_D, 3, 0.0, T1, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    T1, 3, RzfuMs_L, 3, 0.0, T2, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    T2, 3, ms_prm, 3, 0.0, T1, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    T1, 3, RzMs_S, 3, 0.0, T2, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    T2, 3, RxMs_D, 3, 0.0, T1, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    T1, 3, RzMs_L, 3, 0.0, ms_tensor, 3);

        // ms_tensor=(eye(3,3)+ms_tensor);
        cblas_daxpy(3 * 3, 1.0, eye, 1, ms_tensor, 1);

        LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, ms_tensor, 3, ipiv_m);
        LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, ms_tensor, 3, ipiv_m);

        mur_xx = ms_tensor[0 * 3 + 0];
        mur_xy = ms_tensor[0 * 3 + 1];
        mur_xz = ms_tensor[0 * 3 + 2];
        mur_yx = ms_tensor[1 * 3 + 0];
        mur_yy = ms_tensor[1 * 3 + 1];
        mur_yz = ms_tensor[1 * 3 + 2];
        mur_zx = ms_tensor[2 * 3 + 0];
        mur_zy = ms_tensor[2 * 3 + 1];
        mur_zz = ms_tensor[2 * 3 + 2];

        a_xn = (int)fmod(h, NX);
        b_yn = (int)fmod(h - a_xn, NX * NX) / NX;
        c_zn = (int)floor(h / (NX * NY));

        a = fm->A_X[a_xn];
        b = fm->B_Y[b_yn];
        c = fm->C_Z[c_zn];

        getP_xyz(v0Ae);
        cblas_dscal(4 * 4, a * b / 6 / c, v0Ae->PxyTyx, 1);
        cblas_dscal(4 * 4, a / 4, v0Ae->PxzTyx, 1);
        cblas_dscal(4 * 4, a / 4, v0Ae->PxyTzx, 1);
        cblas_dscal(4 * 4, a * c / 6 / b, v0Ae->PxzTzx, 1);

        cblas_dscal(4 * 4, a * b / 4 / c, v0Ae->PyxTyx, 1);
        cblas_dscal(4 * 4, b / 4, v0Ae->PyzTyx, 1);
        cblas_dscal(4 * 4, a / 4, v0Ae->PyxTzx, 1);
        cblas_dscal(4 * 4, c / 6, v0Ae->PyzTzx, 1);

        cblas_dscal(4 * 4, a / 4, v0Ae->PzxTyx, 1);
        cblas_dscal(4 * 4, b / 6, v0Ae->PzyTyx, 1);
        cblas_dscal(4 * 4, a * c / 4 / b, v0Ae->PzxTzx, 1);
        cblas_dscal(4 * 4, c / 4, v0Ae->PzyTzx, 1);

        cblas_dscal(4 * 4, a * b / 4 / c, v0Ae->PxyTxy, 1);
        cblas_dscal(4 * 4, a / 4, v0Ae->PxzTxy, 1);
        cblas_dscal(4 * 4, b / 4, v0Ae->PxyTzy, 1);
        cblas_dscal(4 * 4, c / 6, v0Ae->PxzTzy, 1);

        cblas_dscal(4 * 4, a * b / 6 / c, v0Ae->PyxTxy, 1);
        cblas_dscal(4 * 4, b / 4, v0Ae->PyzTxy, 1);
        cblas_dscal(4 * 4, b / 4, v0Ae->PyxTzy, 1);
        cblas_dscal(4 * 4, b * c / 6 / a, v0Ae->PyzTzy, 1);

        cblas_dscal(4 * 4, a / 6, v0Ae->PzxTxy, 1);
        cblas_dscal(4 * 4, b / 4, v0Ae->PzyTxy, 1);
        cblas_dscal(4 * 4, c / 4, v0Ae->PzxTzy, 1);
        cblas_dscal(4 * 4, b * c / 4 / a, v0Ae->PzyTzy, 1);

        cblas_dscal(4 * 4, a / 4, v0Ae->PxyTxz, 1);
        cblas_dscal(4 * 4, a * c / 4 / b, v0Ae->PxzTxz, 1);
        cblas_dscal(4 * 4, b / 6, v0Ae->PxyTyz, 1);
        cblas_dscal(4 * 4, c / 4, v0Ae->PxzTyz, 1);

        cblas_dscal(4 * 4, a / 6, v0Ae->PyxTxz, 1);
        cblas_dscal(4 * 4, c / 4, v0Ae->PyzTxz, 1);
        cblas_dscal(4 * 4, b / 4, v0Ae->PyxTyz, 1);
        cblas_dscal(4 * 4, b * c / 4 / a, v0Ae->PyzTyz, 1);

        cblas_dscal(4 * 4, a * c / 6 / b, v0Ae->PzxTxz, 1);
        cblas_dscal(4 * 4, c / 4, v0Ae->PzyTxz, 1);
        cblas_dscal(4 * 4, c / 4, v0Ae->PzxTyz, 1);
        cblas_dscal(4 * 4, b * c / 6 / a, v0Ae->PzyTyz, 1);

        cblas_daxpy(4 * 4, mur_yy, v0Ae->PxyTyx, 1, k1exx, 1);
        cblas_daxpy(4 * 4, mur_zy, v0Ae->PxzTyx, 1, k1exx, 1);
        cblas_daxpy(4 * 4, mur_yz, v0Ae->PxyTzx, 1, k1exx, 1);
        cblas_daxpy(4 * 4, mur_zz, v0Ae->PxzTzx, 1, k1exx, 1);

        cblas_daxpy(4 * 4, mur_xy, v0Ae->PyxTyx, 1, k1eyx, 1);
        cblas_daxpy(4 * 4, mur_zy, v0Ae->PyzTyx, 1, k1eyx, 1);
        cblas_daxpy(4 * 4, mur_xz, v0Ae->PyxTzx, 1, k1eyx, 1);
        cblas_daxpy(4 * 4, mur_zz, v0Ae->PyzTzx, 1, k1eyx, 1);

        cblas_daxpy(4 * 4, mur_xy, v0Ae->PzxTyx, 1, k1ezx, 1);
        cblas_daxpy(4 * 4, mur_yy, v0Ae->PzyTyx, 1, k1ezx, 1);
        cblas_daxpy(4 * 4, mur_xz, v0Ae->PzxTzx, 1, k1ezx, 1);
        cblas_daxpy(4 * 4, mur_yz, v0Ae->PzyTzx, 1, k1ezx, 1);

        cblas_daxpy(4 * 4, mur_yx, v0Ae->PxyTxy, 1, k1exy, 1);
        cblas_daxpy(4 * 4, mur_zx, v0Ae->PxzTxy, 1, k1exy, 1);
        cblas_daxpy(4 * 4, mur_yz, v0Ae->PxyTzy, 1, k1exy, 1);
        cblas_daxpy(4 * 4, mur_zz, v0Ae->PxzTzy, 1, k1exy, 1);

        cblas_daxpy(4 * 4, mur_xx, v0Ae->PyxTxy, 1, k1eyy, 1);
        cblas_daxpy(4 * 4, mur_zx, v0Ae->PyzTxy, 1, k1eyy, 1);
        cblas_daxpy(4 * 4, mur_xz, v0Ae->PyxTzy, 1, k1eyy, 1);
        cblas_daxpy(4 * 4, mur_zz, v0Ae->PyzTzy, 1, k1eyy, 1);

        cblas_daxpy(4 * 4, mur_xx, v0Ae->PzxTxy, 1, k1ezy, 1);
        cblas_daxpy(4 * 4, mur_yx, v0Ae->PzyTxy, 1, k1ezy, 1);
        cblas_daxpy(4 * 4, mur_xz, v0Ae->PzxTzy, 1, k1ezy, 1);
        cblas_daxpy(4 * 4, mur_yz, v0Ae->PzyTzy, 1, k1ezy, 1);

        cblas_daxpy(4 * 4, mur_yx, v0Ae->PxyTxz, 1, k1exz, 1);
        cblas_daxpy(4 * 4, mur_zx, v0Ae->PxzTxz, 1, k1exz, 1);
        cblas_daxpy(4 * 4, mur_yy, v0Ae->PxyTyz, 1, k1exz, 1);
        cblas_daxpy(4 * 4, mur_zy, v0Ae->PxzTyz, 1, k1exz, 1);

        cblas_daxpy(4 * 4, mur_xx, v0Ae->PyxTxz, 1, k1eyz, 1);
        cblas_daxpy(4 * 4, mur_zx, v0Ae->PyzTxz, 1, k1eyz, 1);
        cblas_daxpy(4 * 4, mur_xy, v0Ae->PyxTyz, 1, k1eyz, 1);
        cblas_daxpy(4 * 4, mur_zy, v0Ae->PyzTyz, 1, k1eyz, 1);

        cblas_daxpy(4 * 4, mur_xx, v0Ae->PzxTxz, 1, k1ezz, 1);
        cblas_daxpy(4 * 4, mur_yx, v0Ae->PzyTxz, 1, k1ezz, 1);
        cblas_daxpy(4 * 4, mur_xy, v0Ae->PzxTyz, 1, k1ezz, 1);
        cblas_daxpy(4 * 4, mur_yy, v0Ae->PzyTyz, 1, k1ezz, 1);

        mergeKmat(K1e, k1exx, k1exy, k1exz,
                  k1eyx, k1eyy, k1eyz,
                  k1ezx, k1ezy, k1ezz);

        initK2eMat(v0Ae);
        cblas_dscal(4 * 4, sig_xx * prefixVal * a * b * c / 36, v0Ae->k2exx, 1);
        cblas_dscal(4 * 4, sig_xy * prefixVal * a * b * c / 24, v0Ae->k2exy, 1);
        cblas_dscal(4 * 4, sig_xz * prefixVal * a * b * c / 24, v0Ae->k2exz, 1);

        cblas_dscal(4 * 4, sig_yx * prefixVal * a * b * c / 24, v0Ae->k2eyx, 1);
        cblas_dscal(4 * 4, sig_yy * prefixVal * a * b * c / 36, v0Ae->k2eyy, 1);
        cblas_dscal(4 * 4, sig_yz * prefixVal * a * b * c / 24, v0Ae->k2eyz, 1);

        cblas_dscal(4 * 4, sig_zx * prefixVal * a * b * c / 24, v0Ae->k2ezx, 1);
        cblas_dscal(4 * 4, sig_zy * prefixVal * a * b * c / 24, v0Ae->k2ezy, 1);
        cblas_dscal(4 * 4, sig_zz * prefixVal * a * b * c / 36, v0Ae->k2ezz, 1);

        // Global impedance matrix assembly
        mergeKmat(K2e, v0Ae->k2exx, v0Ae->k2exy, v0Ae->k2exz,
                  v0Ae->k2eyx, v0Ae->k2eyy, v0Ae->k2eyz,
                  v0Ae->k2ezx, v0Ae->k2ezy, v0Ae->k2ezz);

        for (j = 0; j < 12; j++)
        {
            NJ = fm->ME[j * NE + h];
            for (k = 0; k < 12; k++)
            {
                NK = fm->ME[k * NE + h];
                val.real(K1e[j * 12 + k]);
                val.imag(K2e[j * 12 + k]);
                MatSetValue(v0Ae->v0, NJ, NK, val, ADD_VALUES);
            }
        }

        freeK1K2(v0Ae);
    }

    MatAssemblyBegin(v0Ae->v0, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(v0Ae->v0, MAT_FINAL_ASSEMBLY);
    PetscTime(&ed_time);
    time = ed_time - st_time;
    PetscPrintf(curComm, "v0Assembly time: %.5lfs\n", time);

    free(sig_prm);
    free(RzfuAl_S);
    free(RxfuAl_D);
    free(RzfuAl_L);
    free(RzAl_S);
    free(RxAl_D);
    free(RzAl_L);
    free(sig_tensor);
    free(ms_prm);
    free(RzfuMs_S);
    free(RxfuMs_D);
    free(RzfuMs_L);
    free(RzMs_S);
    free(RxMs_D);
    free(RzMs_L);
    free(ms_tensor);

    free(k1exx);
    free(k1exy);
    free(k1exz);
    free(k1eyx);
    free(k1eyy);
    free(k1eyz);
    free(k1ezx);
    free(k1ezy);
    free(k1ezz);
    free(K1e);
    free(K2e);
    free(T1);
    free(T2);
    free(eye);
}

void freeK1K2(v0AsEm *v0Ae)
{
    free(v0Ae->PxyTxy);
    free(v0Ae->PxyTxz);
    free(v0Ae->PxyTyx);
    free(v0Ae->PxyTyz);
    free(v0Ae->PxyTzx);
    free(v0Ae->PxyTzy);
    free(v0Ae->PxzTxy);
    free(v0Ae->PxzTxz);
    free(v0Ae->PxzTyx);
    free(v0Ae->PxzTyz);
    free(v0Ae->PxzTzx);
    free(v0Ae->PxzTzy);
    free(v0Ae->PyxTxy);
    free(v0Ae->PyxTxz);
    free(v0Ae->PyxTyx);
    free(v0Ae->PyxTyz);
    free(v0Ae->PyxTzx);
    free(v0Ae->PyxTzy);
    free(v0Ae->PyzTxy);
    free(v0Ae->PyzTxz);
    free(v0Ae->PyzTyx);
    free(v0Ae->PyzTyz);
    free(v0Ae->PyzTzx);
    free(v0Ae->PyzTzy);
    free(v0Ae->PzxTxy);
    free(v0Ae->PzxTxz);
    free(v0Ae->PzxTyx);
    free(v0Ae->PzxTyz);
    free(v0Ae->PzxTzx);
    free(v0Ae->PzxTzy);
    free(v0Ae->PzyTxy);
    free(v0Ae->PzyTxz);
    free(v0Ae->PzyTyx);
    free(v0Ae->PzyTyz);
    free(v0Ae->PzyTzx);
    free(v0Ae->PzyTzy);

    free(v0Ae->k2exx);
    free(v0Ae->k2exy);
    free(v0Ae->k2exz);
    free(v0Ae->k2eyx);
    free(v0Ae->k2eyy);
    free(v0Ae->k2eyz);
    free(v0Ae->k2ezx);
    free(v0Ae->k2ezy);
    free(v0Ae->k2ezz);
}

void freev0Ae(v0AsEm *v0Ae)
{
    MatDestroy(&v0Ae->v0);
}