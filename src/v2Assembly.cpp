#include "../include/MT3D.h"

void v2_Solve(v0AsEm *v0Ae, v2AsEm *v2Ae, Divcorre *div, Fmodel *fm, double freq, MPI_Comm curComm)
{
    int NX, NY, NZ, NP, NL, Nair, Nsea;
    NX = fm->NX;
    NY = fm->NY;
    NZ = fm->NZ;
    NP = fm->NP;
    NL = fm->NL;
    Nair = fm->Nair;
    Nsea = fm->Nsea;
    int curRank, curSize;
    PetscScalar *Ezx2, *Exz2;
    PetscScalar *Eyx2, *Exy2;
    PetscScalar *Ezy2, *Eyz2;
    PetscScalar *Bp;
    int *Bid, *Bidphi;
    int polarization = 2;
    int i, j, k, h, row, col;
    int idn, idn1, idn2, idn3, idn4, idn5, idn6, idn7;
    int BDN = ((NX + 1) * NY + (NY + 1) * NY) * 2 + NZ * (NX + 1 + NX + 1 + NY - 1 + NY - 1 + NX + NX + NY + NY) - NX - NX - NY - NY;
    Vec P, tvE, Xphi, Xe, deltaPhi;
    PetscLogDouble st_time, ed_time, time;
    int slice, st_idx, ed_idx, iters;
    double w = 2 * PI * freq;

    MPI_Comm_size(curComm, &curSize);
    MPI_Comm_rank(curComm, &curRank);
    MatCreate(curComm, &v2Ae->v2);
    MatGetSize(v0Ae->v0, &row, &col);
    MatSetSizes(v2Ae->v2, PETSC_DECIDE, PETSC_DECIDE, row, col);
    MatDuplicate(v0Ae->v0, MAT_COPY_VALUES, &v2Ae->v2);

    VecCreate(curComm, &P);
    VecSetSizes(P, PETSC_DECIDE, NL);
    VecSetUp(P);

    VecCreate(curComm, &tvE);
    VecSetSizes(tvE, PETSC_DECIDE, NP);
    VecSetUp(tvE);

    VecCreate(curComm, &Xphi);
    VecSetSizes(Xphi, PETSC_DECIDE, NP);
    VecSetUp(Xphi);

    VecCreate(curComm, &Xe);
    VecSetSizes(Xe, PETSC_DECIDE, NL);
    VecSetUp(Xe);

    VecCreate(curComm, &deltaPhi);
    VecSetSizes(deltaPhi, PETSC_DECIDE, NL);
    VecSetUp(deltaPhi);

    Bid = (int *)malloc(BDN * sizeof(int));
    Bp = (PetscScalar *)malloc(BDN * sizeof(PetscScalar));

    if (curRank == 0)
    {
        v2Ae->Exx2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        v2Ae->Eyy2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        v2Ae->Ezz2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));

        v2Ae->Hxx2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        v2Ae->Hyy2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        v2Ae->Hzz2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));

        Ezx2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        Exz2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        Eyx2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        Exy2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        Ezy2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        Eyz2 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
    }

    PetscTime(&st_time);
    fDirBdaries(NULL, v2Ae, freq, polarization, fm);

    idn = 0;
    for (j = 0; j < NY + 1; j++)
    {
        for (i = 0; i < NX; i++)
        {
            Bid[idn++] = j * (NX + NX + 1) + i;
        }
    }

    idn1 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn1);
    for (k = 0; k < NZ - 1; k++)
    {
        for (i = 0; i < NX; i++)
        {
            Bid[idn++] = (k + 1) * (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)) + i;
        }
        for (i = 0; i < NX; i++)
        {
            Bid[idn++] = (k + 1) * (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)) + NY * (NX + 1) + NX * NY + i;
        }
    }

    idn2 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn2);
    for (j = 0; j < NY + 1; j++)
    {
        for (i = 0; i < NX; i++)
        {
            Bid[idn++] = NZ * (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)) + j * (NX + NX + 1) + i;
        }
    }

    idn3 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn3);
    for (j = 0; j < NY; j++)
    {
        for (i = 0; i < NX + 1; i++)
        {
            Bid[idn++] = j * (NX + NX + 1) + NX + i;
        }
    }

    idn4 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn4);
    for (k = 0; k < NZ - 1; k++)
    {
        for (j = 0; j < NY; j++)
        {
            Bid[idn++] = (k + 1) * (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)) + j * (NX + NX + 1) + NX;

            Bid[idn++] = (k + 1) * (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)) + j * (NX + NX + 1) + NX + NX;
        }
    }

    idn5 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn5);
    for (j = 0; j < NY; j++)
    {
        for (i = 0; i < NX + 1; i++)
        {
            Bid[idn++] = NZ * (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)) + j * (NX + NX + 1) + NX + i;
        }
    }

    idn6 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn6);
    for (k = 0; k < NZ; k++)
    {
        for (j = 0; j < 1; j++)
        {
            for (i = 0; i < NX + 1; i++)
            {
                Bid[idn++] = k * ((NY + 1) * NX + (NX + 1) * NY + (NX + 1) * (NY + 1)) + (NY + 1) * NX + (NX + 1) * NY + i;
            }
        }

        for (j = 1; j < NY; j++)
        {
            Bid[idn++] = k * ((NY + 1) * NX + (NX + 1) * NY + (NX + 1) * (NY + 1)) + (NY + 1) * NX + (NX + 1) * NY + j * (NX + 1);

            Bid[idn++] = k * ((NY + 1) * NX + (NX + 1) * NY + (NX + 1) * (NY + 1)) + (NY + 1) * NX + (NX + 1) * NY + (j + 1) * (NX + 1) - 1;
        }

        for (j = NY; j < NY + 1; j++)
        {
            for (i = 0; i < NX + 1; i++)
            {
                Bid[idn++] = k * ((NY + 1) * NX + (NX + 1) * NY + (NX + 1) * (NY + 1)) + (NY + 1) * NX + (NX + 1) * NY + (NX + 1) * NY + i;
            }
        }
    }

    idn7 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn7);
    PetscPrintf(PETSC_COMM_WORLD, "第一类边界棱边总数: %d\n", idn7);
    PetscPrintf(PETSC_COMM_WORLD, "A 2\n");

    // 设置边界条件
    idn = 0;
    for (i = 0; i < idn1; i++)
    {
        Bp[idn++] = v2Ae->ExB[0];
    }

    for (i = idn; i < idn2; i++)
    {
        h = floor(Bid[idn] / (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)));
        Bp[idn++] = v2Ae->ExB[h];
    }

    for (i = idn2; i < idn3; i++)
    {
        Bp[idn++] = v2Ae->ExB[NZ];
    }

    for (i = idn3; i < idn4; i++)
    {
        Bp[idn++] = v2Ae->EyB[0];
    }

    for (i = idn4; i < idn5; i++)
    {
        h = floor(Bid[idn] / (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)));
        Bp[idn++] = v2Ae->EyB[h];
    }

    for (i = idn5; i < idn6; i++)
    {
        Bp[idn++] = v2Ae->EyB[NZ];
    }

    for (i = idn6; i < idn7; i++)
    {
        h = floor(Bid[idn] / (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)));
        Bp[idn++] = v2Ae->EzB[h];
    }

    slice = floor(idn7 / curSize);
    st_idx = curRank * slice;
    ed_idx = (curRank != curSize - 1) ? (curRank + 1) * slice : idn7;
    PetscPrintf(PETSC_COMM_WORLD, "A 3\n");
    for (i = st_idx; i < ed_idx; i++)
        VecSetValue(P, Bid[i], Bp[i] * pow(10, 12), INSERT_VALUES);

    VecAssemblyBegin(P);
    VecAssemblyEnd(P);

    PetscPrintf(PETSC_COMM_WORLD, "A 4\n");
    PetscPrintf(PETSC_COMM_WORLD, "A 5\n");
    for (i = st_idx; i < ed_idx; i++)
        MatSetValue(v2Ae->v2, Bid[i], Bid[i], pow(10, 12), INSERT_VALUES);

    MatAssemblyBegin(v2Ae->v2, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(v2Ae->v2, MAT_FINAL_ASSEMBLY);

    PetscPrintf(PETSC_COMM_WORLD, "A change v\n");
    PetscPrintf(PETSC_COMM_WORLD, "A 6\n");

    // 计算表面节点号
    int IDNphi = 0;
    Bidphi = (int *)malloc(sizeof(int));
    for (k = 0; k < NZ + 1; k++)
    {
        for (j = 0; j < NY + 1; j++)
        {
            for (i = 0; i < NX + 1; i++)
            {
                if (k == 0 || k == NZ || j == 0 || j == NY || i == 0 || i == NX)
                {
                    Bidphi = (int *)realloc(Bidphi, sizeof(int) * (IDNphi + 1));
                    Bidphi[IDNphi++] = k * (NX + 1) * (NY + 1) + j * (NX + 1) + i;
                }
            }
        }
    }

    // PetscViewer viewer;
    // PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    // PetscViewerSetType(viewer, PETSCVIEWERASCII);

    // PetscViewerASCIIOpen(PETSC_COMM_WORLD, "res/v2.mtx", &viewer);
    // PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    // MatView(v2Ae->v2, viewer);

    // PetscViewerDestroy(&viewer);

    // int size;
    // VecGetSize(P, &size);
    // FILE *fout = fopen("res/b2.mtx", "w");
    // fprintf(fout, "%%MatrixMarket matrix array complex general\n");
    // fprintf(fout, "%d %d\n", size, 2);
    // for (int i=0; i<size; i++)
    // {
    //     std::complex<double> val;
    //     VecGetValues(P, 1, &i, &val);
    //     fprintf(fout, "%17.16e %17.16e\n", val.real(), val.imag());
    // }
    // fclose(fout);

    VecGetLocalSize(Xe, &slice);
    st_idx = (curRank != curSize - 1) ? curRank * slice : NL - slice;
    ed_idx = (curRank != curSize - 1) ? (curRank + 1) * slice : NL;
    // 节点元、边索引映射
    PetscInt *nodes1, *nodes2, len = 0;
    nodes1 = (PetscInt *)malloc(sizeof(PetscInt) * slice);
    nodes2 = (PetscInt *)malloc(sizeof(PetscInt) * slice);
    Vec sub_Xphi1, sub_Xphi2, distances, tmp_result;
    IS is_node1, is_node2;
    VecCreate(curComm, &distances);
    VecSetSizes(distances, PETSC_DECIDE, NL);
    VecSetUp(distances);

    VecCreate(curComm, &tmp_result);
    VecSetSizes(tmp_result, PETSC_DECIDE, NL);
    VecSetUp(tmp_result);
    VecSet(tmp_result, 0.0);
    for (j = st_idx; j < ed_idx; j++)
    {
        nodes2[len] = fm->EtoN[j * 3 + 1];
        nodes1[len] = fm->EtoN[j * 3 + 0];

        VecSetValue(distances, j, fm->EtoN[j * 3 + 2], INSERT_VALUES);

        len++;
    }
    VecAssemblyBegin(distances);
    VecAssemblyEnd(distances);

    VecReciprocal(distances);

    ISCreateGeneral(curComm, len, nodes1, PETSC_COPY_VALUES, &is_node1);
    ISCreateGeneral(curComm, len, nodes2, PETSC_COPY_VALUES, &is_node2);

    int maxsteps = 100, its;
    double tol = 1e-9;
    PetscReal bnorm, rnorm, relres;
    KSP ksp_x2, ksp_div2;
    PC pc_x2, pc_div2;
    MonitorContext context;

    VecNorm(P, NORM_2, &bnorm);

    // init
    KSPCreate(curComm, &ksp_x2);
    KSPCreate(curComm, &ksp_div2);
    KSPSetOperators(ksp_x2, v2Ae->v2, v2Ae->v2);
    KSPSetOperators(ksp_div2, div->Dv, div->Dv);
    KSPGetPC(ksp_x2, &pc_x2);
    KSPGetPC(ksp_div2, &pc_div2);
    KSPSetType(ksp_x2, KSPBCGS);
    KSPSetType(ksp_div2, KSPQMRCGS);
    PCSetType(pc_x2, PCJACOBI);
    // PCFactorSetMatSolverType(pc_x2, MATSOLVERSUPERLU_DIST);
    PCSetType(pc_div2, PCJACOBI);
    KSPSetTolerances(ksp_x2, tol, PETSC_DEFAULT, PETSC_DEFAULT, maxsteps);
    KSPSetTolerances(ksp_div2, tol, PETSC_DEFAULT, PETSC_DEFAULT, maxsteps);
    KSPSetFromOptions(ksp_x2);
    // KSPSetFromOptions(ksp_div2);

    KSPGetTolerances(ksp_x2, &tol, NULL, NULL, &maxsteps);
    context.rtol = tol;
    context.bnorm = bnorm;
    int total_iters = 2000;
    its = (int)ceil((double)total_iters / maxsteps);
    // Moniter
    // context.minRelativeResidual = PETSC_INFINITY; // 初始化为无穷大
    // KSPMonitorSet(ksp_x2, MyKSPMonitor_yx, &context, NULL);
    // KSPSolve(ksp_x2, P, Xe);
    for (iters = 0; iters < its; iters++)
    {
        KSPSetInitialGuessNonzero(ksp_x2, PETSC_TRUE);
        // Moniter
        // context.minRelativeResidual = PETSC_INFINITY; // 初始化为无穷大
        // KSPMonitorSet(ksp_x2, MyKSPMonitor_yx, &context, NULL);
        KSPSolve(ksp_x2, P, Xe);

        KSPGetResidualNorm(ksp_x2, &rnorm);
        rnorm = context.minRelativeResidual;

        if (rnorm < tol)
        {
            PetscPrintf(curComm, "相对残差 %g 小于容差 %g\n", rnorm, tol);
            break;
        }
        else
            PetscPrintf(curComm, "相对残差为: %g\n", rnorm);

        // tvE = tv*Xe
        MatMult(div->tv, Xe, tvE);
        slice = floor(IDNphi / curSize);
        st_idx = curRank * slice;
        ed_idx = (curRank != curSize - 1) ? (curRank + 1) * slice : IDNphi;
        for (i = st_idx; i < ed_idx; i++)
            VecSetValue(tvE, Bidphi[i], 0, INSERT_VALUES);

        VecAssemblyBegin(tvE);
        VecAssemblyEnd(tvE);

        KSPSolve(ksp_div2, tvE, Xphi);

        VecGetSubVector(Xphi, is_node1, &sub_Xphi1);
        VecGetSubVector(Xphi, is_node2, &sub_Xphi2);

        VecWAXPY(tmp_result, -1.0, sub_Xphi1, sub_Xphi2);

        VecPointwiseMult(deltaPhi, tmp_result, distances);

        VecAXPY(Xe, -1.0, deltaPhi);
    }

    VecDestroy(&sub_Xphi1);
    VecDestroy(&sub_Xphi2);
    free(nodes1);
    free(nodes2);
    VecDestroy(&tmp_result);
    VecDestroy(&distances);
    ISDestroy(&is_node1);
    ISDestroy(&is_node2);
    PCDestroy(&pc_x2);
    PCDestroy(&pc_div2);

    PetscInt Xe_s;
    PetscScalar *Xe_a;

    VecGetLocalSize(Xe, &Xe_s);
    VecGetArray(Xe, &Xe_a);
    PetscInt *recive_size = (PetscInt *)malloc(sizeof(PetscInt) * curSize), displs[curSize];
    PetscScalar *XeVal = (PetscScalar *)malloc(sizeof(PetscScalar) * NL);

    // 收集每个进程的数据大小 size
    MPI_Gather(&Xe_s, 1, MPIU_INT, recive_size, 1, MPIU_INT, 0, curComm);
    // 计算偏移量
    if (!curRank)
    {
        displs[0] = 0;
        for (i = 1; i < curSize; i++)
            displs[i] = displs[i - 1] + recive_size[i - 1];
    }

    // 根据偏移量收集每个进程的Xe数据到0号进程
    MPI_Gatherv(Xe_a, Xe_s, MPIU_SCALAR, XeVal, recive_size, displs, MPIU_SCALAR, 0, curComm);
    free(recive_size);
    VecRestoreArray(Xe, &Xe_a);

    if (curRank == 0)
    {
        PetscScalar Ex2[NX][NY + 1][NZ + 1];
        PetscScalar Ey2[NX + 1][NY][NZ + 1];
        PetscScalar Ez2[NX + 1][NY + 1][NZ];
        int MEij;
        PetscScalar tmpVal;

        // Ex
        for (k = 0; k < NZ + 1; k++)
        {
            for (j = 0; j < NY + 1; j++)
            {
                for (i = 0; i < NX; i++)
                {
                    MEij = k * (NX * (NY + 1) + NY * (NX + 1) + (NX + 1) * (NY + 1)) + j * (2 * NX + 1) + i;
                    tmpVal = XeVal[MEij];
                    Ex2[i][j][k] = tmpVal;
                }
            }
        }
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ex2[i][j + 1][Nair + Nsea] + Ex2[i][j][Nair + Nsea]) / 2.0;
                tmpVal += (Ex2[i][j + 1][Nair + Nsea + 1] + Ex2[i][j][Nair + Nsea + 1]) / 2.0;
                tmpVal = tmpVal / 2.0;

                v2Ae->Exx2[i * NY + j] = tmpVal;
            }
        }

        // Ey
        for (k = 0; k < NZ + 1; k++)
        {
            for (i = 0; i < NX + 1; i++)
            {
                for (j = 0; j < NY; j++)
                {
                    MEij = k * (NX * (NY + 1) + NY * (NX + 1) + (NX + 1) * (NY + 1)) + j * (2 * NX + 1) + NX + i;
                    tmpVal = XeVal[MEij];
                    Ey2[i][j][k] = tmpVal;
                }
            }
        }
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ey2[i][j][Nair + Nsea] + Ey2[i + 1][j][Nair + Nsea]) / 2.0;
                tmpVal += (Ey2[i][j][Nair + Nsea + 1] + Ey2[i + 1][j][Nair + Nsea + 1]) / 2.0;
                tmpVal = tmpVal / 2.0;

                v2Ae->Eyy2[i * NY + j] = tmpVal;
            }
        }

        // Ez
        for (k = 0; k < NZ; k++)
        {
            for (i = 0; i < NX + 1; i++)
            {
                for (j = 0; j < NY + 1; j++)
                {
                    MEij = k * (NX * (NY + 1) + NY * (NX + 1) + (NX + 1) * (NY + 1)) + NX * (NY + 1) + NY * (NX + 1) + j * (NX + 1) + i;
                    tmpVal = XeVal[MEij];
                    Ez2[i][j][k] = tmpVal;
                }
            }
        }
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ez2[i][j][Nair + Nsea] + Ez2[i][j + 1][Nair + Nsea]) / 2.0;
                tmpVal += (Ez2[i + 1][j][Nair + Nsea] + Ez2[i + 1][j + 1][Nair + Nsea]) / 2.0;
                tmpVal = tmpVal / 2.0;

                v2Ae->Ezz2[i * NY + j] = tmpVal;
            }
        }

        // Exz Ezx
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ez2[i + 1][j][Nair + Nsea] - Ez2[i][j][Nair + Nsea]);
                tmpVal += (Ez2[i + 1][j + 1][Nair + Nsea] - Ez2[i][j + 1][Nair + Nsea]);
                tmpVal = tmpVal / 2.0 / fm->A_X[i];

                Ezx2[i * NY + j] = tmpVal;

                tmpVal = (Ex2[i][j][Nair + Nsea + 1] - Ex2[i][j][Nair + Nsea]);
                tmpVal += (Ex2[i][j + 1][Nair + Nsea + 1] - Ex2[i][j + 1][Nair + Nsea]);
                tmpVal = tmpVal / 2.0 / fm->C_Z[Nair + Nsea];

                Exz2[i * NY + j] = tmpVal;
            }
        }

        // Eyx Exy
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ey2[i + 1][j][Nair + Nsea] - Ey2[i][j][Nair + Nsea]);
                tmpVal += (Ey2[i + 1][j][Nair + Nsea + 1] - Ey2[i][j][Nair + Nsea + 1]);
                tmpVal = tmpVal / 2.0 / fm->A_X[i];

                Eyx2[i * NY + j] = tmpVal;

                tmpVal = (Ex2[i][j + 1][Nair + Nsea] - Ex2[i][j][Nair + Nsea]);
                tmpVal += (Ex2[i][j + 1][Nair + Nsea + 1] - Ex2[i][j][Nair + Nsea + 1]);
                tmpVal = tmpVal / 2.0 / fm->B_Y[j];

                Exy2[i * NY + j] = tmpVal;
            }
        }

        // Ezy Eyz
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ez2[i][j + 1][Nair + Nsea] - Ez2[i][j][Nair + Nsea]);
                tmpVal += (Ez2[i + 1][j + 1][Nair + Nsea] - Ez2[i + 1][j][Nair + Nsea]);
                tmpVal = tmpVal / 2.0 / fm->B_Y[j];

                Ezy2[i * NY + j] = tmpVal;

                tmpVal = (Ey2[i + 1][j][Nair + Nsea + 1] - Ey2[i + 1][j][Nair + Nsea]);
                tmpVal += (Ey2[i][j][Nair + Nsea + 1] - Ey2[i][j][Nair + Nsea]);
                tmpVal = tmpVal / 2.0 / fm->C_Z[Nair + Nsea];

                Eyz2[i * NY + j] = tmpVal;
            }
        }

        // Hxyz
        std::complex<double> w_mu0(0, w * mu0);

        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = 1.0 / w_mu0 * (Ezy2[i * NY + j] - Eyz2[i * NY + j]);
                v2Ae->Hxx2[i * NY + j] = tmpVal;

                tmpVal = 1.0 / w_mu0 * (Exz2[i * NY + j] - Ezx2[i * NY + j]);
                v2Ae->Hyy2[i * NY + j] = tmpVal;

                v2Ae->Hzz2[i * NY + j] = 0.0;
            }
        }
    }

    PetscTime(&ed_time);
    time = ed_time - st_time;
    PetscPrintf(curComm, "v2Assembly time: %.5lfs\n", time);

    if (curRank == 0)
    {
        free(Ezx2);
        free(Exz2);
        free(Eyx2);
        free(Exy2);
        free(Ezy2);
        free(Eyz2);
    }

    free(Bid);
    free(Bp);

    VecDestroy(&P);
    VecDestroy(&tvE);
    VecDestroy(&Xphi);
    VecDestroy(&Xe);
    VecDestroy(&deltaPhi);
}

void freev2Ae(v2AsEm *v2Ae, MPI_Comm curComm)
{
    free(v2Ae->ExB);
    free(v2Ae->EyB);
    free(v2Ae->EzB);
    int curRank;
    MPI_Comm_rank(curComm, &curRank);
    if (curRank == 0)
    {
        free(v2Ae->Exx2);
        free(v2Ae->Eyy2);
        free(v2Ae->Ezz2);
        free(v2Ae->Hxx2);
        free(v2Ae->Hyy2);
        free(v2Ae->Hzz2);
    }

    MatDestroy(&v2Ae->v2);
}