#include "../include/MT3D.h"

void mergeMatrix(double *merged_te, double *te1, double *te2, double *te3){
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 4; j++) {
            merged_te[i * 12 + j] = te1[i * 4 + j];
            merged_te[i * 12 + (j + 4)] = te2[i * 4 + j];
            merged_te[i * 12 + (j + 2 * 4)] = te3[i * 4 + j];
        }
}

void getMEs_v(Fmodel *fm, Divcorre *div) {
    int NX, NY, NZ, NE;
    int L, M, N, IL, col_index;
    NX = fm->NX; NY = fm->NY; NZ = fm->NZ; NE = fm->NE;
    div->MEs = (int*)malloc(8*NE * sizeof(int));
    div->MEv = (int*)malloc(12*NE * sizeof(int));

    for (L = 0; L < NZ; L++)
        for (M = 0; M < NY; M++)
            for(N = 0; N < NX; N++) {
                IL = L*(NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1));
                col_index = L*NX*NY+M*NX+N;

                div->MEs[0*NE+col_index] = L*(NY+1)*(NX+1)+M*(NX+1)+N;
                div->MEs[1*NE+col_index] = L*(NY+1)*(NX+1)+M*(NX+1)+N+1;
                div->MEs[2*NE+col_index] = L*(NY+1)*(NX+1)+M*(NX+1)+N+NX+2;
                div->MEs[3*NE+col_index] = L*(NY+1)*(NX+1)+M*(NX+1)+N+NX+1;
                div->MEs[4*NE+col_index] = L*(NY+1)*(NX+1)+M*(NX+1)+N+(NX+1)*(NY+1);
                div->MEs[5*NE+col_index] = L*(NY+1)*(NX+1)+M*(NX+1)+N+1+(NX+1)*(NY+1);
                div->MEs[6*NE+col_index] = L*(NY+1)*(NX+1)+M*(NX+1)+N+NX+2+(NX+1)*(NY+1);
                div->MEs[7*NE+col_index] = L*(NY+1)*(NX+1)+M*(NX+1)+N+NX+1+(NX+1)*(NY+1);

                div->MEv[0*NE+col_index] = IL+M*(2*NX+1)+N;
                div->MEv[1*NE+col_index] = IL+M*(2*NX+1)+N+2*NX+1;
                div->MEv[2*NE+col_index] = IL+M*(2*NX+1)+N+NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1);
                div->MEv[3*NE+col_index] = IL+M*(2*NX+1)+N+NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1)+2*NX+1;
                div->MEv[4*NE+col_index] = IL+M*(2*NX+1)+N+NX;
                div->MEv[5*NE+col_index] = IL+M*(2*NX+1)+N+NX+NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1);
                div->MEv[6*NE+col_index] = IL+M*(2*NX+1)+N+NX+1;
                div->MEv[7*NE+col_index] = IL+M*(2*NX+1)+N+NX+1+NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1);
                div->MEv[8*NE+col_index] = IL+NX*(NY+1)+NY*(NX+1)+M*(NX+1)+N;
                div->MEv[9*NE+col_index] = IL+NX*(NY+1)+NY*(NX+1)+M*(NX+1)+N+1;
                div->MEv[10*NE+col_index] = IL+NX*(NY+1)+NY*(NX+1)+M*(NX+1)+N+1+NX;
                div->MEv[11*NE+col_index] = IL+NX*(NY+1)+NY*(NX+1)+M*(NX+1)+N+2+NX;
            }
}

void init_Dfactor(Divcorre *div) {
    double tmp_df1[8*8] = {4, -4, -2,  2,  2, -2, -1,  1,
            -4,  4,  2, -2, -2,  2,  1, -1,
            -2,  2,  4, -4, -1,  1,  2, -2,
            2, -2, -4,  4,  1, -1, -2,  2,
            2, -2, -1,  1,  4, -4, -2,  2,
            -2,  2,  1, -1, -4,  4,  2, -2,
            -1,  1,  2, -2, -2,  2,  4, -4,
            1, -1, -2,  2,  2, -2, -4,  4};

    double tmp_df2[8*8] = {6, -6, -6,  6,  3, -3, -3,  3,
            6, -6, -6,  6,  3, -3, -3,  3,
            -6,  6,  6, -6, -3,  3,  3, -3,
            -6,  6,  6, -6, -3,  3,  3, -3,
            3, -3, -3,  3,  6, -6, -6,  6,
            3, -3, -3,  3,  6, -6, -6,  6,
            -3,  3,  3, -3, -6,  6,  6, -6,
            -3,  3,  3, -3, -6,  6,  6, -6};
    
    double tmp_df3[8*8] = {6, -6, -3,  3,  6, -6, -3,  3,
            6, -6, -3,  3,  6, -6, -3,  3,
            3, -3, -6,  6,  3, -3, -6,  6,
            3, -3, -6,  6,  3, -3, -6,  6,
            -6,  6,  3, -3, -6,  6,  3, -3,
            -6,  6,  3, -3, -6,  6,  3, -3,
            -3,  3,  6, -6, -3,  3,  6, -6,
            -3,  3,  6, -6, -3,  3,  6, -6};

    double tmp_df4[8*8] = {6,  6, -6, -6,  3,  3, -3, -3,
           -6, -6,  6,  6, -3, -3,  3,  3,
           -6, -6,  6,  6, -3, -3,  3,  3,
            6,  6, -6, -6,  3,  3, -3, -3,
            3,  3, -3, -3,  6,  6, -6, -6,
           -3, -3,  3,  3, -6, -6,  6,  6,
           -3, -3,  3,  3, -6, -6,  6,  6,
            3,  3, -3, -3,  6,  6, -6, -6};

    double tmp_df5[8*8] = {4,  2, -2, -4,  2,  1, -1, -2,
            2,  4, -4, -2,  1,  2, -2, -1,
           -2, -4,  4,  2, -1, -2,  2,  1,
           -4, -2,  2,  4, -2, -1,  1,  2,
            2,  1, -1, -2,  4,  2, -2, -4,
            1,  2, -2, -1,  2,  4, -4, -2,
           -1, -2,  2,  1, -2, -4,  4,  2,
           -2, -1,  1,  2, -4, -2,  2,  4};

    double tmp_df6[8*8] = {6,  3, -3, -6,  6,  3, -3, -6,
            3,  6, -6, -3,  3,  6, -6, -3,
            3,  6, -6, -3,  3,  6, -6, -3,
            6,  3, -3, -6,  6,  3, -3, -6,
           -6, -3,  3,  6, -6, -3,  3,  6,
           -3, -6,  6,  3, -3, -6,  6,  3,
           -3, -6,  6,  3, -3, -6,  6,  3,
           -6, -3,  3,  6, -6, -3,  3,  6};

    double tmp_df7[8*8] = {6,  6,  3,  3, -6, -6, -3, -3,
           -6, -6, -3, -3,  6,  6,  3,  3,
           -3, -3, -6, -6,  3,  3,  6,  6,
            3,  3,  6,  6, -3, -3, -6, -6,
            6,  6,  3,  3, -6, -6, -3, -3,
           -6, -6, -3, -3,  6,  6,  3,  3,
           -3, -3, -6, -6,  3,  3,  6,  6,
            3,  3,  6,  6, -3, -3, -6, -6};

    double tmp_df8[8*8] = {6,  3,  3,  6, -6, -3, -3, -6,
            3,  6,  6,  3, -3, -6, -6, -3,
           -3, -6, -6, -3,  3,  6,  6,  3,
           -6, -3, -3, -6,  6,  3,  3,  6,
            6,  3,  3,  6, -6, -3, -3, -6,
            3,  6,  6,  3, -3, -6, -6, -3,
           -3, -6, -6, -3,  3,  6,  6,  3,
           -6, -3, -3, -6,  6,  3,  3,  6};

    double tmp_df9[8*8] = {4,  2,  1,  2, -4, -2, -1, -2,
            2,  4,  2,  1, -2, -4, -2, -1,
            1,  2,  4,  2, -1, -2, -4, -2,
            2,  1,  2,  4, -2, -1, -2, -4,
           -4, -2, -1, -2,  4,  2,  1,  2,
           -2, -4, -2, -1,  2,  4,  2,  1,
           -1, -2, -4, -2,  1,  2,  4,  2,
           -2, -1, -2, -4,  2,  1,  2,  4};
    
    div->df1 = (double*)malloc(8*8 * sizeof(double));
    div->df2 = (double*)malloc(8*8 * sizeof(double));
    div->df3 = (double*)malloc(8*8 * sizeof(double));
    div->df4 = (double*)malloc(8*8 * sizeof(double));
    div->df5 = (double*)malloc(8*8 * sizeof(double));
    div->df6 = (double*)malloc(8*8 * sizeof(double));
    div->df7 = (double*)malloc(8*8 * sizeof(double));
    div->df8 = (double*)malloc(8*8 * sizeof(double));
    div->df9 = (double*)malloc(8*8 * sizeof(double));

    memcpy(div->df1, tmp_df1, 8*8 * sizeof(double));
    memcpy(div->df2, tmp_df2, 8*8 * sizeof(double));
    memcpy(div->df3, tmp_df3, 8*8 * sizeof(double));
    memcpy(div->df4, tmp_df4, 8*8 * sizeof(double));
    memcpy(div->df5, tmp_df5, 8*8 * sizeof(double));
    memcpy(div->df6, tmp_df6, 8*8 * sizeof(double));
    memcpy(div->df7, tmp_df7, 8*8 * sizeof(double));
    memcpy(div->df8, tmp_df8, 8*8 * sizeof(double));
    memcpy(div->df9, tmp_df9, 8*8 * sizeof(double));
}

void init_tfactor (Divcorre *div) {
    double tmp_tf1[8*4] = {-4, -2, -2, -1,
            4,  2,  2,  1,
            2,  4,  1,  2,
           -2, -4, -1, -2,
           -2, -1, -4, -2,
            2,  1,  4,  2,
            1,  2,  2,  4,
           -1, -2, -2, -4};

    double tmp_tf2[8*4] = {-4, -2, -2, -1,
           -2, -1, -4, -2,
            2,  1,  4,  2,
            4,  2,  2,  1,
           -2, -4, -1, -2,
           -1, -2, -2, -4,
            1,  2,  2,  4,
            2,  4,  1,  2};

    double tmp_tf3[8*4] = {-4, -2, -2, -1,
           -2, -4, -1, -2,
           -1, -2, -2, -4,
           -2, -1, -4, -2,
            4,  2,  2,  1,
            2,  4,  1,  2,
            1,  2,  2,  4,
            2,  1,  4,  2};

    div->tf1 = (double*)malloc(8*4 * sizeof(double));
    div->tf2 = (double*)malloc(8*4 * sizeof(double));
    div->tf3 = (double*)malloc(8*4 * sizeof(double));

    memcpy(div->tf1, tmp_tf1, 8*4 * sizeof(double));
    memcpy(div->tf2, tmp_tf2, 8*4 * sizeof(double));
    memcpy(div->tf3, tmp_tf3, 8*4 * sizeof(double));
}

void divFE(Fmodel *fm, Divcorre *div, MPI_Comm curComm) {
    int NX, NY, NE, NP, NL;
    int curRank, curSize;
    NX = fm->NX; NY = fm->NY; NE = fm->NE; 
    NP = fm->NP; NL = fm->NL;
    MPI_Comm_size(curComm, &curSize);
    MPI_Comm_rank(curComm, &curRank);
    getMEs_v(fm, div);

    //Dv tv
    MatCreate(curComm, &div->Dv);
    MatSetSizes(div->Dv, PETSC_DECIDE, PETSC_DECIDE, NP, NP);

    MatCreate(curComm, &div->tv);
    MatSetSizes(div->tv, PETSC_DECIDE, PETSC_DECIDE, NP, NL);

    double *sig_prm, *RzfuAl_S, *RxfuAl_D, *RzfuAl_L;
    double *RzAl_S, *RxAl_D, *RzAl_L, *sig_tensor;
    double *T1, *T2;
    double *De, *te;
    double *te1, *te2, *te3;
    double sig_xx, sig_xy, sig_xz;
    double sig_yx, sig_yy, sig_yz;
    double sig_zx, sig_zy, sig_zz;
    double r_1, r_2, r_3;
    double a_S, a_D, a_L;
    int a_xn, b_yn, c_zn;
    double a, b, c;
    int NJ, NKs, NKv;
    int h, slice, st_idx, ed_idx, i, j, k;
    std::complex<double> val;
    PetscLogDouble st_time, ed_time, time;
    sig_prm = (double*)malloc(3*3 * sizeof(double));
    RzfuAl_S = (double*)malloc(3*3 * sizeof(double));
    RxfuAl_D = (double*)malloc(3*3 * sizeof(double));
    RzfuAl_L = (double*)malloc(3*3 * sizeof(double));
    RzAl_S = (double*)malloc(3*3 * sizeof(double));
    RxAl_D = (double*)malloc(3*3 * sizeof(double));
    RzAl_L = (double*)malloc(3*3 * sizeof(double));
    sig_tensor = (double*)malloc(3*3 * sizeof(double));

    T1 = (double*)malloc(3*3 * sizeof(double));
    T2 = (double*)malloc(3*3 * sizeof(double));

    De = (double*)malloc(8*8 * sizeof(double));
    te = (double*)malloc(8*12 * sizeof(double));
    te1 = (double*)malloc(8*4 * sizeof(double));
    te2 = (double*)malloc(8*4 * sizeof(double));
    te3 = (double*)malloc(8*4 * sizeof(double));

    memset(sig_prm, 0, 3*3*sizeof(double));memset(RzfuAl_S, 0, 3*3*sizeof(double));
    memset(RxfuAl_D, 0, 3*3*sizeof(double));memset(RzfuAl_L, 0, 3*3*sizeof(double));
    memset(RzAl_S, 0, 3*3*sizeof(double));memset(RxAl_D, 0, 3*3*sizeof(double));
    memset(RzAl_L, 0, 3*3*sizeof(double));memset(sig_tensor, 0, 3*3*sizeof(double));

    PetscTime(&st_time);
    slice = floor((double)NE/curSize);
    st_idx = curRank*slice;
    ed_idx = (curRank!=curSize-1) ? (curRank+1)*slice : NE;
    for (h = st_idx; h < ed_idx; h++) {
        memset(T1, 0, 3*3 * sizeof(double));
        memset(T2, 0, 3*3 * sizeof(double));
        memset(sig_tensor, 0, 3*3 * sizeof(double));
        memset(De, 0, 8*8 * sizeof(double));
        memset(te, 0, 8*12 * sizeof(double));
        memset(te1, 0, 8*4 * sizeof(double));
        memset(te2, 0, 8*4 * sizeof(double));
        memset(te3, 0, 8*4 * sizeof(double));

        r_1 = fm->rho[h*3+0];r_2 = fm->rho[h*3+1];r_3 = fm->rho[h*3+2];
        a_S = fm->alpha_S[h];a_D = fm->alpha_D[h];a_L = fm->alpha_L[h];

        sig_prm[0*3+0] = 1.0/(r_1);
        sig_prm[1*3+1] = 1.0/(r_2);
        sig_prm[2*3+2] = 1.0/(r_3);

        RzfuAl_S[0*3+0] = cos(-a_S);
        RzfuAl_S[0*3+1] = sin(-a_S);
        RzfuAl_S[1*3+0] = -sin(-a_S);
        RzfuAl_S[1*3+1] = cos(-a_S);
        RzfuAl_S[2*3+2] = 1;

        RxfuAl_D[0*3+0] = 1;
        RxfuAl_D[1*3+1] = cos(-a_D);
        RxfuAl_D[1*3+2] = sin(-a_D);
        RxfuAl_D[2*3+1] = -sin(-a_D);
        RxfuAl_D[2*3+2] = cos(-a_D);

        RzfuAl_L[0*3+0] = cos(-a_L);
        RzfuAl_L[0*3+1] = sin(-a_L);
        RzfuAl_L[1*3+0] = -sin(-a_L);
        RzfuAl_L[1*3+1] = cos(-a_L);
        RzfuAl_L[2*3+2] = 1;

        RzAl_S[0*3+0] = cos(a_S);
        RzAl_S[0*3+1] = sin(a_S);
        RzAl_S[1*3+0] = -sin(a_S);
        RzAl_S[1*3+1] = cos(a_S);
        RzAl_S[2*3+2] = 1;

        RxAl_D[0*3+0] = 1;
        RxAl_D[1*3+1] = cos(a_D);
        RxAl_D[1*3+2] = sin(a_D);
        RxAl_D[2*3+1] = -sin(a_D);
        RxAl_D[2*3+2] = cos(a_D);

        RzAl_L[0*3+0] = cos(a_L);
        RzAl_L[0*3+1] = sin(a_L);
        RzAl_L[1*3+0] = -sin(a_L);
        RzAl_L[1*3+1] = cos(a_L);
        RzAl_L[2*3+2] = 1;

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    RzfuAl_S, 3, RxfuAl_D, 3, 0.0, T1, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, 
                    T1, 3, RzfuAl_L, 3, 0.0, T2, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, 
                    T2, 3, sig_prm, 3, 0.0, T1, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, 
                    T1, 3, RzAl_L, 3, 0.0, T2, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, 
                    T2, 3, RxAl_D, 3, 0.0, T1, 3);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, 
                    T1, 3, RzAl_S, 3, 0.0, sig_tensor, 3);

        sig_xx=sig_tensor[0*3+0];sig_xy=sig_tensor[0*3+1];sig_xz=sig_tensor[0*3+2];
        sig_yx=sig_tensor[1*3+0];sig_yy=sig_tensor[1*3+1];sig_yz=sig_tensor[1*3+2];
        sig_zx=sig_tensor[2*3+0];sig_zy=sig_tensor[2*3+1];sig_zz=sig_tensor[2*3+2];

        a_xn = (int)fmod(h, NX);
        b_yn = (int)fmod(h-a_xn, NX*NX)/NX;
        c_zn = (int)floor(h/(NX*NY));

        a=fm->A_X[a_xn];b=fm->B_Y[b_yn];c=fm->C_Z[c_zn];

        /*
            b*c/a/36*sigma_xx*Dfactor1+c/72*sigma_yx*Dfactor2+b/72*sigma_zx*Dfactor3+...
            c/72*sigma_xy*Dfactor4+a*c/b/36*sigma_yy*Dfactor5+a/72*sigma_zy*Dfactor6+...
            b/72*sigma_xz*Dfactor7+a/72*sigma_yz*Dfactor8+a*b/c/36*sigma_zz*Dfactor9;
        */
        init_Dfactor(div);
        cblas_daxpy(8*8, b*c/a/36*sig_xx, div->df1, 1, De, 1);
        cblas_daxpy(8*8, c/72*sig_yx, div->df2, 1, De, 1);
        cblas_daxpy(8*8, b/72*sig_zx, div->df3, 1, De, 1);
        cblas_daxpy(8*8, c/72*sig_xy, div->df4, 1, De, 1);
        cblas_daxpy(8*8, a*c/b/36*sig_yy, div->df5, 1, De, 1);
        cblas_daxpy(8*8, a/72*sig_zy, div->df6, 1, De, 1);
        cblas_daxpy(8*8, b/72*sig_xz, div->df7, 1, De, 1);
        cblas_daxpy(8*8, a/72*sig_yz, div->df8, 1, De, 1);
        cblas_daxpy(8*8, a*b/c/36*sig_zz, div->df9, 1, De, 1);

        /*
            sigma_xx*b*c/36*tfactor1+sigma_yx*a*c/36*tfactor2+sigma_zx*a*b/36*tfactor3 ...
            sigma_xy*b*c/36*tfactor1+sigma_yy*a*c/36*tfactor2+sigma_zy*a*b/36*tfactor3 ...
            sigma_xz*b*c/36*tfactor1+sigma_yz*a*c/36*tfactor2+sigma_zz*a*b/36*tfactor3;
        */
        init_tfactor(div);
        cblas_daxpy(8*4, sig_xx*b*c/36, div->tf1, 1, te1, 1);
        cblas_daxpy(8*4, sig_yx*a*c/36, div->tf2, 1, te1, 1);
        cblas_daxpy(8*4, sig_zx*a*b/36, div->tf3, 1, te1, 1);

        cblas_daxpy(8*4, sig_xy*b*c/36, div->tf1, 1, te2, 1);
        cblas_daxpy(8*4, sig_yy*a*c/36, div->tf2, 1, te2, 1);
        cblas_daxpy(8*4, sig_zy*a*b/36, div->tf3, 1, te2, 1);

        cblas_daxpy(8*4, sig_xz*b*c/36, div->tf1, 1, te3, 1);
        cblas_daxpy(8*4, sig_yz*a*c/36, div->tf2, 1, te3, 1);
        cblas_daxpy(8*4, sig_zz*a*b/36, div->tf3, 1, te3, 1);
        mergeMatrix(te, te1, te2, te3);

        for (j = 0; j < 8; j++) {
            NJ = div->MEs[j*NE+h];
            for (k = 0; k < 8; k++) {
                NKs = div->MEs[k*NE+h];
                val.real(De[j*8+k]);
                val.imag(0.0);
                MatSetValue(div->Dv, NJ, NKs, 
                val, ADD_VALUES);
            }
            for (i = 0; i < 12; i++) {
                NKv = div->MEv[i*NE+h];
                val.real(te[j*12+i]);
                val.imag(0.0);
                MatSetValue(div->tv, NJ, NKv,
                val, ADD_VALUES);
            }
        }

        freedf_tf(div);
    }
    MatAssemblyBegin(div->Dv, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(div->tv, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(div->Dv, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(div->tv, MAT_FINAL_ASSEMBLY);
    PetscTime(&ed_time);
    time = ed_time - st_time;
    PetscPrintf(curComm, "Divcorre time: %.5lfs\n", time);

    free(sig_prm);free(RzfuAl_S);free(RxfuAl_D);free(RzfuAl_L);
    free(RzAl_S);free(RxAl_D);free(RzAl_L);free(sig_tensor);
    free(T1);free(T2);free(De);free(te);
    free(te1);free(te2);free(te3);

}

void freedf_tf(Divcorre *div) {
    free(div->df1);free(div->df2);free(div->df3);
    free(div->df4);free(div->df5);free(div->df6);
    free(div->df7);free(div->df8);free(div->df9);
    free(div->tf1);free(div->tf2);free(div->tf3);
}

void freeDiv(Divcorre *div) {
    free(div->MEs);free(div->MEv);
    MatDestroy(&div->Dv);MatDestroy(&div->tv);
}