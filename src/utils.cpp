#include "../include/MT3D.h"

void edgeTonode(Fmodel *fm)
{
    int NX, NY, NZ, NL;
    NX = fm->NX;
    NY = fm->NY;
    NZ = fm->NZ;
    fm->EtoN = (int *)malloc(fm->NL * 3 * sizeof(int));
    int k, j, i;
    int IL, IN1, IN2;

    // 对Ex方向的边线进行循环
    for (k = 0; k < NZ + 1; k++)
        for (j = 0; j < NY + 1; j++)
            for (i = 0; i < NX; i++)
            {
                // 计算边索引IL
                IL = k * (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)) + j * (2 * NX + 1) + i;
                // 根据网格定位计算对应两个节点索引IN1、IN2
                IN1 = k * (NX + 1) * (NY + 1) + j * (NX + 1) + i;
                IN2 = k * (NX + 1) * (NY + 1) + j * (NX + 1) + i + 1;
                // 将节点索引和边长A_X填入EtoN相应位置
                fm->EtoN[IL * 3 + 0] = IN1;
                fm->EtoN[IL * 3 + 1] = IN2;
                fm->EtoN[IL * 3 + 2] = fm->A_X[i];
            }

    // 对Ey方向进行边线循环
    for (k = 0; k < NZ + 1; k++)
        for (j = 0; j < NY; j++)
            for (i = 0; i < NX + 1; i++)
            {
                IL = k * (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)) + j * (2 * NX + 1) + NX + i;
                IN1 = k * (NX + 1) * (NY + 1) + j * (NX + 1) + i;
                IN2 = k * (NX + 1) * (NY + 1) + j * (NX + 1) + i + NX + 1;
                // 将节点索引和边长B_Y填入EtoN相应位置
                fm->EtoN[IL * 3 + 0] = IN1;
                fm->EtoN[IL * 3 + 1] = IN2;
                fm->EtoN[IL * 3 + 2] = fm->B_Y[j];
            }

    // 对Ez方向进行边线循环
    for (k = 0; k < NZ; k++)
        for (j = 0; j < NY + 1; j++)
            for (i = 0; i < NX + 1; i++)
            {
                IL = k * (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)) + NX * (NY + 1) + NY * (NX + 1) + j * (NX + 1) + i;
                IN1 = k * (NX + 1) * (NY + 1) + j * (NX + 1) + i;
                IN2 = k * (NX + 1) * (NY + 1) + j * (NX + 1) + i + (NX + 1) * (NY + 1);
                // 将节点索引和边长C_Z填入EtoN相应位置
                fm->EtoN[IL * 3 + 0] = IN1;
                fm->EtoN[IL * 3 + 1] = IN2;
                fm->EtoN[IL * 3 + 2] = fm->C_Z[k];
            }
    return;
}

FILE *readModelFile(char *filename, double freq, int isloadB, Fmodel *fm)
{
    char prefix[100] = "./model_data/";
    char str_model[100];
    sprintf(str_model, "data-%d%d%d/", fm->NX, fm->NY, fm->NZ);
    strcat(prefix, str_model);
    char str_freq[100];
    if (isloadB)
    {
        sprintf(str_freq, "data-%.4f/", freq);
        strcat(prefix, str_freq);
    }
    strcat(prefix, filename);

    FILE *file;
    file = fopen(prefix, "rb");
    if (file == NULL)
    {
        printf("无法打开文件 %s\n", prefix);
        return NULL;
    }

    return file;
}

double getRadian(int degree)
{
    switch (degree)
    {
    case 0:
        return 0.0;
    case 10:
        return (acos(-1) / 18.0);
    case 20:
        return (acos(-1) / 9.0);
    case 30:
        return (acos(-1) / 6.0);
    case 45:
        return (acos(-1) / 4.0);
    case 60:
        return (acos(-1) / 3.0);
    case 90:
        return (acos(-1) / 2.0);
    default:
        break;
    }

    return 0.0;
}

void load_data(double *real, double *imag, char *filename1, char *filename2, double freq, Fmodel *fmodel)
{
    int count;
    FILE *file;
    char line[200];

    count = 0;
    file = readModelFile(filename1, freq, 1, fmodel);
    while (fgets(line, sizeof(line), file) != NULL)
    {
        real[count++] = strtod(line, NULL);
    }
    fclose(file);

    count = 0;
    file = readModelFile(filename2, freq, 1, fmodel);
    while (fgets(line, sizeof(line), file) != NULL)
    {
        imag[count++] = strtod(line, NULL);
    }
    fclose(file);
}

void fDirBdaries(v1AsEm *v1Ae, v2AsEm *v2Ae, double freq, int polarization, Fmodel *fm)
{
    char line[200];
    char P2REx[20] = "P2REx.dat";
    char P2IEx[20] = "P2IEx.dat";
    char P2REy[20] = "P2REy.dat";
    char P2IEy[20] = "P2IEy.dat";
    char P2REz[20] = "P2REz.dat";
    char P2IEz[20] = "P2IEz.dat";
    char P2RHx[20] = "P2RHx.dat";
    char P2IHx[20] = "P2IHx.dat";
    char P2RHy[20] = "P2RHy.dat";
    char P2IHy[20] = "P2IHy.dat";

    char P1REx[20] = "P1REx.dat";
    char P1IEx[20] = "P1IEx.dat";
    char P1REy[20] = "P1REy.dat";
    char P1IEy[20] = "P1IEy.dat";
    char P1REz[20] = "P1REz.dat";
    char P1IEz[20] = "P1IEz.dat";
    char P1RHx[20] = "P1RHx.dat";
    char P1IHx[20] = "P1IHx.dat";
    char P1RHy[20] = "P1RHy.dat";
    char P1IHy[20] = "P1IHy.dat";

    double *real, *imag;
    int NZ, Nair, Nsea;
    int i;
    std::complex<double> fDval;
    NZ = fm->NZ;
    Nair = fm->Nair;
    Nsea = fm->Nsea;
    real = (double *)malloc(sizeof(double) * ((NZ + 1) + 1));
    imag = (double *)malloc(sizeof(double) * ((NZ + 1) + 1));

    if (polarization == 1)
    {
        v1Ae->ExA = (PetscScalar *)malloc((NZ + 1) * sizeof(PetscScalar));
        v1Ae->EyA = (PetscScalar *)malloc((NZ + 1) * sizeof(PetscScalar));
        v1Ae->EzA = (PetscScalar *)malloc((NZ + 1) * sizeof(PetscScalar));
        // v1Ae->HxA = (PetscScalar *)malloc((NZ + 1) * sizeof(PetscScalar));
        // v1Ae->HyA = (PetscScalar *)malloc((NZ + 1) * sizeof(PetscScalar));
        // v1Ae->HzA = (double *)malloc((NZ + 1) * sizeof(double));

        // ExA
        load_data(real, imag, P2REx, P2IEx, freq, fm);

        for (i = 0; i < NZ + 1; i++)
        {
            fDval.real(real[i]);
            fDval.imag(imag[i]);

            v1Ae->ExA[i] = fDval;
        }

        // EyA
        load_data(real, imag, P2REy, P2IEy, freq, fm);

        for (i = 0; i < NZ + 1; i++)
        {
            fDval.real(real[i]);
            fDval.imag(imag[i]);

            v1Ae->EyA[i] = fDval;
        }

        // EzA
        load_data(real, imag, P2REz, P2IEz, freq, fm);

        for (i = 0; i < NZ + 1; i++)
        {
            fDval.real(real[i]);
            fDval.imag(imag[i]);

            v1Ae->EzA[i] = fDval;
        }
        // file_r.open(prefix + std::string("/P2RHx.dat"));
        // file_i.open(prefix + std::string("/P2IHx.dat"));

        // count = 0;
        // while (std::getline(file_r, line_r) && std::getline(file_i, line_i))
        // {
        //     fDval.real(std::stod(line_r));
        //     fDval.imag(std::stod(line_i));

        //     v1Ae->HxA[count++] = fDval;
        // }

        // file_r.close();
        // file_i.close();

        // file_r.open(prefix + std::string("/P2RHy.dat"));
        // file_i.open(prefix + std::string("/P2IHy.dat"));

        // count = 0;
        // while (std::getline(file_r, line_r) && std::getline(file_i, line_i))
        // {
        //     fDval.real(std::stod(line_r));
        //     fDval.imag(std::stod(line_i));

        //     v1Ae->HyA[count++] = fDval;
        // }

        // file_r.close();
        // file_i.close();

        // for (int i = 0; i < NZ + 1; i++)
        //     v1Ae->HzA[i] = 0.0;

        // v1Ae->rho_a1 = std::pow((std::abs(v1Ae->EyA[Nair + Nsea] / v1Ae->HxA[0])), 2) / w / mu0;
        // v1Ae->rho_a2 = std::pow((std::abs(v1Ae->ExA[Nair + Nsea] / v1Ae->HyA[0])), 2) / w / mu0;
    }
    else if (polarization == 2)
    {
        std::ifstream file_r, file_i;

        v2Ae->ExB = (PetscScalar *)malloc((NZ + 1) * sizeof(PetscScalar));
        v2Ae->EyB = (PetscScalar *)malloc((NZ + 1) * sizeof(PetscScalar));
        v2Ae->EzB = (PetscScalar *)malloc((NZ + 1) * sizeof(PetscScalar));

        // ExB
        load_data(real, imag, P1REx, P1IEx, freq, fm);

        for (i = 0; i < NZ + 1; i++)
        {
            fDval.real(real[i]);
            fDval.imag(imag[i]);

            v2Ae->ExB[i] = fDval;
        }

        // EyB
        load_data(real, imag, P1REy, P1IEy, freq, fm);

        for (i = 0; i < NZ + 1; i++)
        {
            fDval.real(real[i]);
            fDval.imag(imag[i]);

            v2Ae->EyB[i] = fDval;
        }

        // EzB
        load_data(real, imag, P1REz, P1IEz, freq, fm);

        for (i = 0; i < NZ + 1; i++)
        {
            fDval.real(real[i]);
            fDval.imag(imag[i]);

            v2Ae->EzB[i] = fDval;
        }
    }
}

// 自定义监视器函数
PetscErrorCode MyKSPMonitor_xy(KSP ksp, PetscInt n, PetscReal rnorm, void *ctx)
{
    PetscReal relres;
    MonitorContext *context = (MonitorContext*)ctx;
    relres = rnorm;
    if (n == 0) {
        context->minRelativeResidual = relres;
    } else {
        if (relres < context->minRelativeResidual) {
            context->minRelativeResidual = relres;
        }
    }
    // std::ofstream file("res/xymode_iters_div_bcgs.txt", std::ios::app);
    // file << n << " " << relres << "\n";
    // file.close();
    if (relres < context->rtol) {
        PetscPrintf(PETSC_COMM_WORLD, "相对残差 %g 小于容差 %g\n", relres, context->rtol);

        // KSPConvergedReasonSet(ksp, KSP_CONVERGED_RTOL);  // 可选：设置收敛标志
        return 1;  // 强制终止迭代
    }
    return 0;
}

// 自定义监视器函数
PetscErrorCode MyKSPMonitor_yx(KSP ksp, PetscInt n, PetscReal rnorm, void *ctx)
{
    PetscReal relres;
    MonitorContext *context = (MonitorContext*)ctx;
    relres = rnorm;
    if (n == 0) {
        context->minRelativeResidual = relres;
    } else {
        if (relres < context->minRelativeResidual) {
            context->minRelativeResidual = relres;
        }
    }
    std::ofstream file("res/0.0100_yxmode_iters_bcgs.txt", std::ios::app);
    file << n << " " << relres << "\n";
    file.close();
    if (relres < context->rtol) {
        PetscPrintf(PETSC_COMM_WORLD, "相对残差 %g 小于容差 %g\n", relres, context->rtol);

        // KSPConvergedReasonSet(ksp, KSP_CONVERGED_RTOL);  // 可选：设置收敛标志
        return 1;  // 强制终止迭代
    }
    return 0;
}

void isExit(char* prefix){
    struct stat st;
    if(!(stat(prefix, &st) == 0) || !S_ISDIR(st.st_mode))
        mkdir(prefix, 0777);
}

void rhoAndpha(v1AsEm *v1Ae, v2AsEm *v2Ae, Fmodel *fm, double freq, MPI_Comm curComm)
{
    int NX, NY;
    NX = fm->NX;
    NY = fm->NY;
    int curRank;
    MPI_Comm_rank(curComm, &curRank);
    PetscLogDouble st_time, ed_time, time;
    if (curRank == 0)
    {
        PetscTime(&st_time);

        int i, j;
        double w = 2 * PI * freq;
        std::complex<double> tmpVal, sigVal;
        sigVal.real(0.0);
        sigVal.imag(1.0);
        PetscScalar *zxx, *zxy, *zyx, *zyy;
        zxx = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        zxy = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        zyx = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        zyy = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));

        double PhaseXX[NX][NY], PhaseXY[NX][NY], PhaseYX[NX][NY], PhaseYY[NX][NY];
        double ResXX[NX][NY], ResXY[NX][NY], ResYX[NX][NY], ResYY[NX][NY];

        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = v1Ae->Hxx1[i * NY + j] * v2Ae->Hyy2[i * NY + j] - v1Ae->Hyy1[i * NY + j] * v2Ae->Hxx2[i * NY + j];
                zxx[i * NY + j] = (v1Ae->Exx1[i * NY + j] * v2Ae->Hyy2[i * NY + j] - v2Ae->Exx2[i * NY + j] * v1Ae->Hyy1[i * NY + j]) / tmpVal;
                zxy[i * NY + j] = (v2Ae->Exx2[i * NY + j] * v1Ae->Hxx1[i * NY + j] - v1Ae->Exx1[i * NY + j] * v2Ae->Hxx2[i * NY + j]) / tmpVal;
                zyx[i * NY + j] = (v1Ae->Eyy1[i * NY + j] * v2Ae->Hyy2[i * NY + j] - v2Ae->Eyy2[i * NY + j] * v1Ae->Hyy1[i * NY + j]) / tmpVal;
                zyy[i * NY + j] = (v2Ae->Eyy2[i * NY + j] * v1Ae->Hxx1[i * NY + j] - v1Ae->Eyy1[i * NY + j] * v2Ae->Hxx2[i * NY + j]) / tmpVal;

                PhaseXX[i][j] = -atan(zxx[i * NY + j].imag() / zxx[i * NY + j].real()) * 180 / PI;
                PhaseXY[i][j] = -atan(zxy[i * NY + j].imag() / zxy[i * NY + j].real()) * 180 / PI;
                PhaseYX[i][j] = -atan(zyx[i * NY + j].imag() / zyx[i * NY + j].real()) * 180 / PI;
                PhaseYY[i][j] = -atan(zyy[i * NY + j].imag() / zyy[i * NY + j].real()) * 180 / PI;
                ResXX[i][j] = std::abs(std::pow(zxx[i * NY + j], 2) * sigVal / (w * mu0));
                ResXY[i][j] = std::abs(std::pow(zxy[i * NY + j], 2) * sigVal / (w * mu0));
                ResYX[i][j] = std::abs(std::pow(zyx[i * NY + j], 2) * sigVal / (w * mu0));
                ResYY[i][j] = std::abs(std::pow(zyy[i * NY + j], 2) * sigVal / (w * mu0));
            }
        }

        int midx = floor(NX / 2), midy = floor(NY / 2);
        PetscPrintf(curComm, "PhaseXX mid: %g\n", PhaseXX[midx][midy]);
        PetscPrintf(curComm, "PhaseXY mid: %g\n", PhaseXY[midx][midy]);
        PetscPrintf(curComm, "PhaseYX mid: %g\n", PhaseYX[midx][midy]);
        PetscPrintf(curComm, "PhaseYY mid: %g\n", PhaseYY[midx][midy]);

        PetscPrintf(curComm, "RseXX mid: %g\n", ResXX[midx][midy]);
        PetscPrintf(curComm, "RseXY mid: %g\n", ResXY[midx][midy]);
        PetscPrintf(curComm, "RseYX mid: %g\n", ResYX[midx][midy]);
        PetscPrintf(curComm, "RseYY mid: %g\n", ResYY[midx][midy]);

        char fPhaXX[20] = "PhaseXX.txt";
        char fPhaXY[20] = "PhaseXY.txt";
        char fPhaYX[20] = "PhaseYX.txt";
        char fPhaYY[20] = "PhaseYY.txt";

        char fResXX[20] = "ResXX.txt";
        char fResXY[20] = "ResXY.txt";
        char fResYX[20] = "ResYX.txt";
        char fResYY[20] = "ResYY.txt";

        char prefix[50] = "res/";
        char str_freq[20], tmp_prefix[50];
        sprintf(str_freq, "%.4f_bcgs/", freq);
        strcat(prefix, str_freq);
        isExit(prefix);

        FILE *filePXX, *filePXY, *filePYX, *filePYY;
        FILE *fileRXX, *fileRXY, *fileRYX, *fileRYY;

        strcpy(tmp_prefix, prefix);
        strcat(tmp_prefix, fPhaXX);
        filePXX = fopen(tmp_prefix, "wb");

        strcpy(tmp_prefix, prefix);
        strcat(tmp_prefix, fPhaXY);
        filePXY = fopen(tmp_prefix, "wb");

        strcpy(tmp_prefix, prefix);
        strcat(tmp_prefix, fPhaYX);
        filePYX = fopen(tmp_prefix, "wb");

        strcpy(tmp_prefix, prefix);
        strcat(tmp_prefix, fPhaYY);
        filePYY = fopen(tmp_prefix, "wb");

        strcpy(tmp_prefix, prefix);
        strcat(tmp_prefix, fResXX);
        fileRXX = fopen(tmp_prefix, "wb");

        strcpy(tmp_prefix, prefix);
        strcat(tmp_prefix, fResXY);
        fileRXY = fopen(tmp_prefix, "wb");

        strcpy(tmp_prefix, prefix);
        strcat(tmp_prefix, fResYX);
        fileRYX = fopen(tmp_prefix, "wb");

        strcpy(tmp_prefix, prefix);
        strcat(tmp_prefix, fResYY);
        fileRYY = fopen(tmp_prefix, "wb");

        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                fprintf(filePXX, "%lf\t", PhaseXX[i][j]);
                fprintf(filePXY, "%lf\t", PhaseXY[i][j]);
                fprintf(filePYX, "%lf\t", PhaseYX[i][j]);
                fprintf(filePYY, "%lf\t", PhaseYY[i][j]);

                fprintf(fileRXX, "%lf\t", ResXX[i][j]);
                fprintf(fileRXY, "%lf\t", ResXY[i][j]);
                fprintf(fileRYX, "%lf\t", ResYX[i][j]);
                fprintf(fileRYY, "%lf\t", ResYY[i][j]);
            }
            fprintf(filePXX, "\n");
            fprintf(filePXY, "\n");
            fprintf(filePYX, "\n");
            fprintf(filePYY, "\n");

            fprintf(fileRXX, "\n");
            fprintf(fileRXY, "\n");
            fprintf(fileRYX, "\n");
            fprintf(fileRYY, "\n");
        }
        fclose(filePXX);
        fclose(filePXY);
        fclose(filePYX);
        fclose(filePYY);

        fclose(fileRXX);
        fclose(fileRXY);
        fclose(fileRYX);
        fclose(fileRYY);
    }
}