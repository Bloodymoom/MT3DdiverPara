#include "../include/MT3D.h"

void v1_Solve(v0AsEm *v0Ae, v1AsEm *v1Ae, Divcorre *div, Fmodel *fm, double freq, MPI_Comm curComm)
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
    PetscScalar *Ezx1, *Exz1;
    PetscScalar *Eyx1, *Exy1;
    PetscScalar *Ezy1, *Eyz1;
    PetscScalar *Bp;
    int *Bid, *Bidphi;
    int polarization = 1;
    int i, j, k, h, row, col;
    int idn, idn1, idn2, idn3, idn4, idn5, idn6, idn7;
    int BDN = ((NX + 1) * NY + (NY + 1) * NY) * 2 + NZ * (NX + 1 + NX + 1 + NY - 1 + NY - 1 + NX + NX + NY + NY) - NX - NX - NY - NY;
    Vec P, tvE, Xphi, Xe, deltaPhi;
    PetscLogDouble st_time, ed_time, time;
    int slice, st_idx, ed_idx, iters;
    double w = 2 * PI * freq;

    MPI_Comm_size(curComm, &curSize);
    MPI_Comm_rank(curComm, &curRank);
    MatCreate(curComm, &v1Ae->v1);
    MatGetSize(v0Ae->v0, &row, &col);
    MatSetSizes(v1Ae->v1, PETSC_DECIDE, PETSC_DECIDE, row, col);
    MatDuplicate(v0Ae->v0, MAT_COPY_VALUES, &v1Ae->v1);

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
        v1Ae->Exx1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        v1Ae->Eyy1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        v1Ae->Ezz1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));

        v1Ae->Hxx1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        v1Ae->Hyy1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        v1Ae->Hzz1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));

        Ezx1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        Exz1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        Eyx1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        Exy1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        Ezy1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
        Eyz1 = (PetscScalar *)malloc(NX * NY * sizeof(PetscScalar));
    }

    PetscTime(&st_time);
    fDirBdaries(v1Ae, NULL, freq, polarization, fm);

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
        Bp[idn++] = v1Ae->ExA[0];
    }

    for (i = idn; i < idn2; i++)
    {
        h = floor(Bid[idn] / (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)));
        Bp[idn++] = v1Ae->ExA[h];
    }

    for (i = idn2; i < idn3; i++)
    {
        Bp[idn++] = v1Ae->ExA[NZ];
    }

    for (i = idn3; i < idn4; i++)
    {
        Bp[idn++] = v1Ae->EyA[0];
    }

    for (i = idn4; i < idn5; i++)
    {
        h = floor(Bid[idn] / (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)));
        Bp[idn++] = v1Ae->EyA[h];
    }

    for (i = idn5; i < idn6; i++)
    {
        Bp[idn++] = v1Ae->EyA[NZ];
    }

    for (i = idn6; i < idn7; i++)
    {
        h = floor(Bid[idn] / (NX * (NY + 1) + (NX + 1) * NY + (NX + 1) * (NY + 1)));
        Bp[idn++] = v1Ae->EzA[h];
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
        MatSetValue(v1Ae->v1, Bid[i], Bid[i], pow(10, 12), INSERT_VALUES);

    MatAssemblyBegin(v1Ae->v1, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(v1Ae->v1, MAT_FINAL_ASSEMBLY);

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

    // 加入外边界
    slice = floor(IDNphi / curSize);
    st_idx = curRank * slice;
    ed_idx = (curRank != curSize - 1) ? (curRank + 1) * slice : IDNphi;
    for (i = st_idx; i < ed_idx; i++)
        MatSetValue(div->Dv, Bidphi[i], Bidphi[i], pow(10, 10), INSERT_VALUES);

    MatAssemblyBegin(div->Dv, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(div->Dv, MAT_FINAL_ASSEMBLY);

    // PetscViewer viewer;
    // PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    // PetscViewerSetType(viewer, PETSCVIEWERASCII);

    // // PetscViewerASCIIOpen(PETSC_COMM_WORLD, "res/Dv.mtx", &viewer);
    // // PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    // // MatView(div->Dv, viewer);

    // // PetscViewerASCIIOpen(PETSC_COMM_WORLD, "res/tv.mtx", &viewer);
    // // PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    // // MatView(div->tv, viewer);

    // PetscViewerASCIIOpen(PETSC_COMM_WORLD, "res/v1.mtx", &viewer);
    // PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    // MatView(v1Ae->v1, viewer);

    // PetscViewerDestroy(&viewer);

    // int size;
    // VecGetSize(P, &size);
    // FILE *fout = fopen("res/b1.mtx", "w");
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
    KSP ksp_x1, ksp_div1;
    PC pc_x1, pc_div1;
    MonitorContext context;

    VecNorm(P, NORM_2, &bnorm);

    // init
    KSPCreate(curComm, &ksp_x1);
    KSPCreate(curComm, &ksp_div1);
    KSPSetOperators(ksp_x1, v1Ae->v1, v1Ae->v1);
    KSPSetOperators(ksp_div1, div->Dv, div->Dv);
    KSPGetPC(ksp_x1, &pc_x1);
    KSPGetPC(ksp_div1, &pc_div1);
    KSPSetType(ksp_x1, KSPBCGS);
    KSPSetType(ksp_div1, KSPQMRCGS);
    PCSetType(pc_x1, PCJACOBI);
    // PCFactorSetMatSolverType(pc_x1, MATSOLVERSUPERLU_DIST);
    PCSetType(pc_div1, PCJACOBI);
    KSPSetTolerances(ksp_x1, tol, PETSC_DEFAULT, PETSC_DEFAULT, maxsteps);
    KSPSetTolerances(ksp_div1, tol, PETSC_DEFAULT, PETSC_DEFAULT, maxsteps);
    KSPSetFromOptions(ksp_x1);
    // KSPSetFromOptions(ksp_div1);

    KSPGetTolerances(ksp_x1, &tol, NULL, NULL, &maxsteps);
    context.rtol = tol;
    context.bnorm = bnorm;
    int total_iters = 2000;
    its = (int)ceil((double)total_iters / maxsteps);
    // context.minRelativeResidual = PETSC_INFINITY; // 初始化为无穷大
    // KSPMonitorSet(ksp_x1, MyKSPMonitor_xy, &context, NULL);
    // KSPSolve(ksp_x1, P, Xe);
    for (iters = 0; iters < its; iters++)
    {
        KSPSetInitialGuessNonzero(ksp_x1, PETSC_TRUE);
        // Moniter
        context.minRelativeResidual = PETSC_INFINITY; // 初始化为无穷大
        KSPMonitorSet(ksp_x1, MyKSPMonitor_xy, &context, NULL);
        KSPSolve(ksp_x1, P, Xe);

        KSPGetResidualNorm(ksp_x1, &rnorm);
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

        KSPSolve(ksp_div1, tvE, Xphi);

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
    PCDestroy(&pc_x1);
    PCDestroy(&pc_div1);

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
        PetscScalar Ex1[NX][NY + 1][NZ + 1];
        PetscScalar Ey1[NX + 1][NY][NZ + 1];
        PetscScalar Ez1[NX + 1][NY + 1][NZ];
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
                    Ex1[i][j][k] = tmpVal;
                }
            }
        }

        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ex1[i][j + 1][Nair + Nsea] + Ex1[i][j][Nair + Nsea]) / 2.0;
                tmpVal += (Ex1[i][j + 1][Nair + Nsea + 1] + Ex1[i][j][Nair + Nsea + 1]) / 2.0;
                tmpVal = tmpVal / 2.0;

                v1Ae->Exx1[i * NY + j] = tmpVal;
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
                    Ey1[i][j][k] = tmpVal;
                }
            }
        }
        std::ofstream fout("res/div-sourceA-45.txt");
        for (i=0; i<NY; i++)
            fout << Ey1[33][i][50].real() << ", ";
        fout.close();

        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ey1[i][j][Nair + Nsea] + Ey1[i + 1][j][Nair + Nsea]) / 2.0;
                tmpVal += (Ey1[i][j][Nair + Nsea + 1] + Ey1[i + 1][j][Nair + Nsea + 1]) / 2.0;
                tmpVal = tmpVal / 2.0;

                v1Ae->Eyy1[i * NY + j] = tmpVal;
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
                    Ez1[i][j][k] = tmpVal;
                }
            }
        }
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ez1[i][j][Nair + Nsea] + Ez1[i][j + 1][Nair + Nsea]) / 2.0;
                tmpVal += (Ez1[i + 1][j][Nair + Nsea] + Ez1[i + 1][j + 1][Nair + Nsea]) / 2.0;
                tmpVal = tmpVal / 2.0;

                v1Ae->Ezz1[i * NY + j] = tmpVal;
            }
        }

        // Exz Ezx
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ez1[i + 1][j][Nair + Nsea] - Ez1[i][j][Nair + Nsea]);
                tmpVal += (Ez1[i + 1][j + 1][Nair + Nsea] - Ez1[i][j + 1][Nair + Nsea]);
                tmpVal = tmpVal / 2.0 / fm->A_X[i];

                Ezx1[i * NY + j] = tmpVal;

                tmpVal = (Ex1[i][j][Nair + Nsea + 1] - Ex1[i][j][Nair + Nsea]);
                tmpVal += (Ex1[i][j + 1][Nair + Nsea + 1] - Ex1[i][j + 1][Nair + Nsea]);
                tmpVal = tmpVal / 2.0 / fm->C_Z[Nair + Nsea];

                Exz1[i * NY + j] = tmpVal;
            }
        }

        // Eyx Exy
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ey1[i + 1][j][Nair + Nsea] - Ey1[i][j][Nair + Nsea]);
                tmpVal += (Ey1[i + 1][j][Nair + Nsea + 1] - Ey1[i][j][Nair + Nsea + 1]);
                tmpVal = tmpVal / 2.0 / fm->A_X[i];

                Eyx1[i * NY + j] = tmpVal;

                tmpVal = (Ex1[i][j + 1][Nair + Nsea] - Ex1[i][j][Nair + Nsea]);
                tmpVal += (Ex1[i][j + 1][Nair + Nsea + 1] - Ex1[i][j][Nair + Nsea + 1]);
                tmpVal = tmpVal / 2.0 / fm->B_Y[j];

                Exy1[i * NY + j] = tmpVal;
            }
        }

        // Ezy Eyz
        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = (Ez1[i][j + 1][Nair + Nsea] - Ez1[i][j][Nair + Nsea]);
                tmpVal += (Ez1[i + 1][j + 1][Nair + Nsea] - Ez1[i + 1][j][Nair + Nsea]);
                tmpVal = tmpVal / 2.0 / fm->B_Y[j];

                Ezy1[i * NY + j] = tmpVal;

                tmpVal = (Ey1[i + 1][j][Nair + Nsea + 1] - Ey1[i + 1][j][Nair + Nsea]);
                tmpVal += (Ey1[i][j][Nair + Nsea + 1] - Ey1[i][j][Nair + Nsea]);
                tmpVal = tmpVal / 2.0 / fm->C_Z[Nair + Nsea];

                Eyz1[i * NY + j] = tmpVal;
            }
        }

        // Hxyz
        std::complex<double> w_mu0(0, w * mu0);

        for (i = 0; i < NX; i++)
        {
            for (j = 0; j < NY; j++)
            {
                tmpVal = 1.0 / w_mu0 * (Ezy1[i * NY + j] - Eyz1[i * NY + j]);
                v1Ae->Hxx1[i * NY + j] = tmpVal;

                tmpVal = 1.0 / w_mu0 * (Exz1[i * NY + j] - Ezx1[i * NY + j]);
                v1Ae->Hyy1[i * NY + j] = tmpVal;

                v1Ae->Hzz1[i * NY + j] = 0.0;
            }
        }
    }

    PetscTime(&ed_time);
    time = ed_time - st_time;
    PetscPrintf(curComm, "v1Assembly time: %.5lfs\n", time);

    if (curRank == 0)
    {
        free(Ezx1);
        free(Exz1);
        free(Eyx1);
        free(Exy1);
        free(Ezy1);
        free(Eyz1);
    }

    free(Bid);
    free(Bp);

    VecDestroy(&P);
    VecDestroy(&tvE);
    VecDestroy(&Xphi);
    VecDestroy(&Xe);
    VecDestroy(&deltaPhi);
}

void freev1Ae(v1AsEm *v1Ae, MPI_Comm curComm)
{
    free(v1Ae->ExA);
    free(v1Ae->EyA);
    free(v1Ae->EzA);
    // free(v1Ae->HxA);
    // free(v1Ae->HyA);
    // free(v1Ae->HzA);
    int curRank;
    MPI_Comm_rank(curComm, &curRank);
    if (curRank == 0)
    {
        free(v1Ae->Exx1);
        free(v1Ae->Eyy1);
        free(v1Ae->Ezz1);
        free(v1Ae->Hxx1);
        free(v1Ae->Hyy1);
        free(v1Ae->Hzz1);
    }

    MatDestroy(&v1Ae->v1);
}