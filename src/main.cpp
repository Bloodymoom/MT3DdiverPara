#include "../include/MT3D.h"

int main(int argc, char **args)
{
    PetscErrorCode ierr;
    int N_PROC, rank, MPI_FEM, MPI_FRE, FRE_NUM;
    int Subdomain_Id = 0, key;
    int counter = 0;
    int fnn, f_start, f_end;
    MPI_Comm curComm;
    PetscLogDouble total_st, total_ed, total,
        st1, ed1, total1,
        st2, ed2, total2;
    Fmodel *fm;
    Divcorre *div;
    v0AsEm *v0Ae;
    v1AsEm *v1Ae;
    v2AsEm *v2Ae;

    double *freqs;
    char line[1024];
    char *token;
    FILE *freqs_file;
    freqs = (double *)malloc(sizeof(double));
    // 0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1.0,5.0,10.0,100.0,200.0,500.0,1000.0,2000.0
    freqs_file = fopen("./model_data/freqs.txt", "rb");
    if (freqs_file == NULL)
    {
        std::cerr << "file is not exist." << std::endl;
        return 0;
    }

    if (fgets(line, 1024, freqs_file) != NULL)
    {
        token = strtok(line, ",");
        while (token != NULL)
        {
            counter++;
            freqs = (double *)realloc(freqs, sizeof(double) * counter);
            freqs[counter - 1] = strtod(token, NULL);
            token = strtok(NULL, ","); 
        }
    }
    FRE_NUM = counter - 1;

    fclose(freqs_file);

    ierr = PetscInitialize(&argc, &args, NULL, "Program help message or NULL");
    if (ierr)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Error during initialization: %d\n", ierr);
        return ierr;
    }

    PetscTime(&total_st);

    // 划分子域
    MPI_Comm_size(PETSC_COMM_WORLD, &N_PROC);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_FRE = (int)freqs[FRE_NUM]; // 每个子域处理频率数据个数
    if (N_PROC == 1)
    {
        MPI_FRE = FRE_NUM;
        MPI_FEM = FRE_NUM / MPI_FRE;
    }
    else if (FRE_NUM % MPI_FRE == 0)
    {
        MPI_FEM = FRE_NUM / MPI_FRE;
        if (N_PROC < MPI_FEM)
        {
            PetscPrintf(PETSC_COMM_WORLD, "error process number.\n");
            return 0;
        }
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD, "MPI Split wrong.\n");
        return 0;
    }
    Subdomain_Id = rank % MPI_FEM;
    key = rank;
    MPI_Comm_split(PETSC_COMM_WORLD, Subdomain_Id, key, &curComm);

    PetscTime(&st1);
    fm = (Fmodel *)malloc(sizeof(Fmodel));
    init_Fmodel(fm, freqs[MPI_FRE * Subdomain_Id]);
    edgeTonode(fm);
    PetscTime(&ed1);
    total1 = ed1 - st1;
    PetscPrintf(curComm, "Fmodel time: %.5lfs\n", total1);

    div = (Divcorre *)malloc(sizeof(Divcorre));
    divFE(fm, div, curComm);

    v0Ae = (v0AsEm *)malloc(sizeof(v0AsEm));
    v1Ae = (v1AsEm *)malloc(sizeof(v1AsEm));
    v2Ae = (v2AsEm *)malloc(sizeof(v2AsEm));

    f_start = Subdomain_Id * MPI_FRE;
    f_end = f_start + MPI_FRE;
    PetscPrintf(PETSC_COMM_WORLD, "total %d Source\n", FRE_NUM);
    for (fnn = f_start; fnn < f_end; fnn++)
    {
        PetscTime(&st2);
        PetscPrintf(curComm, "第%d个频点: %f\n", fnn + 1, freqs[fnn]);
        v0_FE(freqs[fnn], fm, v0Ae, curComm);
        v1_Solve(v0Ae, v1Ae, div, fm, freqs[fnn], curComm);
        // v2_Solve(v0Ae, v2Ae, div, fm, freqs[fnn], curComm);
        // rhoAndpha(v1Ae, v2Ae, fm, freqs[fnn], curComm);
        PetscTime(&ed2);
        total2 = ed2 - st2;
        PetscPrintf(curComm, "第%d个频点 time: %.5lfs\n", fnn + 1, total2);

        freev0Ae(v0Ae);
        freev1Ae(v1Ae, curComm);
        // freev2Ae(v2Ae, curComm);
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    freeDiv(div);
    freeModel(fm);

    PetscTime(&total_ed);
    total = total_ed - total_st;
    PetscPrintf(PETSC_COMM_WORLD, "total time: %.5lfs\n", total);

    ierr = PetscFinalize();
    if (ierr)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Error during finalization: %d\n", ierr);
        return ierr;
    }

    return 0;
}