#include "../include/MT3D.h"

void init_Fmodel(Fmodel *fm, double freq)
{
    FILE *file;
    char line[200];
    char model[20] = "model.dat";
    char prefix[100] = "./model_data/";
    char A_X[20] = "A_X.dat";
    char B_Y[20] = "B_Y.dat";
    char C_Z[20] = "C_Z.dat";
    char *token;
    int h, idx, x, y, z, a_bday, s_bday, counter = 0;
    double rho_air, rho_sea, rho_under;
    int anomalies, L, M, N;
    int *z_d, *y_d, *x_d, *alpha_d, *alpha_ms_d;
    double *rho_d, *ms_d;
    int IL, col_index;

    strcat(prefix, model);
    file = fopen(prefix, "rb");
    if (file == NULL)
    {
        printf("file can not open %s.\n", prefix);
        return;
    }
    idx = 1;
    while (fgets(line, sizeof(line), file) != NULL)
    {
        switch (idx)
        {
        case 1:
            sscanf(line, "%d %d %d %d %d", &fm->NX, &fm->NY, &fm->NZ, &fm->Nair, &fm->Nsea);
            break;
        case 2:
            sscanf(line, "%lf %lf %lf", &rho_air, &rho_sea, &rho_under);
            break;
        case 3:
            sscanf(line, "%d", &anomalies);
            z_d = (int *)malloc((2 * anomalies) * sizeof(int));
            y_d = (int *)malloc((2 * anomalies) * sizeof(int));
            x_d = (int *)malloc((2 * anomalies) * sizeof(int));
            alpha_d = (int *)malloc((3 * anomalies) * sizeof(int));
            alpha_ms_d = (int *)malloc((3 * anomalies) * sizeof(int));
            rho_d = (double *)malloc((3 * anomalies) * sizeof(double));
            ms_d = (double *)malloc((3 * anomalies) * sizeof(double));
            break;
        case 4:
            counter = 0;
            token = strtok(line, " ");
            while (token != NULL && strcmp(token, "#"))
            {
                z_d[counter++] = atoi(token);
                token = strtok(NULL, " ");
            }
            break;
        case 5:
            counter = 0;
            token = strtok(line, " ");
            while (token != NULL && strcmp(token, "#"))
            {
                y_d[counter++] = atoi(token);
                token = strtok(NULL, " ");
            }
            break;
        case 6:
            counter = 0;
            token = strtok(line, " ");
            while (token != NULL && strcmp(token, "#"))
            {
                x_d[counter++] = atoi(token);
                token = strtok(NULL, " ");
            }
            break;
        case 7:
            counter = 0;
            token = strtok(line, " ");
            while (token != NULL && strcmp(token, "#"))
            {
                rho_d[counter++] = strtod(token, NULL);
                token = strtok(NULL, " ");
            }
            break;
        case 8:
            counter = 0;
            token = strtok(line, " ");
            while (token != NULL && strcmp(token, "#"))
            {
                alpha_d[counter++] = atoi(token);
                token = strtok(NULL, " ");
            }
            break;
        case 9:
            counter = 0;
            token = strtok(line, " ");
            while (token != NULL && strcmp(token, "#"))
            {
                ms_d[counter++] = strtod(token, NULL);
                token = strtok(NULL, " ");
            }
            break;
        case 10:
            counter = 0;
            token = strtok(line, " ");
            while (token != NULL && strcmp(token, "#"))
            {
                alpha_ms_d[counter++] = atoi(token);
                token = strtok(NULL, " ");
            }
            break;
        default:
            break;
        }
        idx++;
    }

    fclose(file);

    x = fm->NX;
    y = fm->NY;
    z = fm->NZ;
    fm->NE = x * y * z;
    fm->NP = (x + 1) * (y + 1) * (z + 1);
    fm->NL = (x * (y + 1) * (z + 1) + (x + 1) * y * (z + 1) + (x + 1) * (y + 1) * z);
    // ######################
    a_bday = fm->Nair * x * y;
    s_bday = (fm->Nair + fm->Nsea) * x * y;
    fm->rho = (double *)malloc(fm->NE * 3 * sizeof(double));
    fm->ms = (double *)malloc(fm->NE * 3 * sizeof(double));
    fm->alpha_S = (double *)malloc(fm->NE * sizeof(double));
    fm->alpha_D = (double *)malloc(fm->NE * sizeof(double));
    fm->alpha_L = (double *)malloc(fm->NE * sizeof(double));
    fm->ms_S = (double *)malloc(fm->NE * sizeof(double));
    fm->ms_D = (double *)malloc(fm->NE * sizeof(double));
    fm->ms_L = (double *)malloc(fm->NE * sizeof(double));
    fm->ME = (int *)malloc(fm->NE * 12 * sizeof(int));
    memset(fm->alpha_S, 0, fm->NE * sizeof(double));
    memset(fm->alpha_D, 0, fm->NE * sizeof(double));
    memset(fm->alpha_L, 0, fm->NE * sizeof(double));
    memset(fm->ms_S, 0, fm->NE * sizeof(double));
    memset(fm->ms_D, 0, fm->NE * sizeof(double));
    memset(fm->ms_L, 0, fm->NE * sizeof(double));
    memset(fm->ms, 0, fm->NE * 3 * sizeof(double));

    for (h = 0; h < fm->NE; h++)
    {
        if (h < a_bday)
        {
            fm->rho[h * 3 + 0] = rho_air;
            fm->rho[h * 3 + 1] = rho_air;
            fm->rho[h * 3 + 2] = rho_air;
        }
        else if (fm->Nsea != 0 && h < s_bday && (h > a_bday | h == a_bday))
        {
            fm->rho[h * 3 + 0] = rho_sea;
            fm->rho[h * 3 + 1] = rho_sea;
            fm->rho[h * 3 + 2] = rho_sea;
        }
        else
        {
            fm->rho[h * 3 + 0] = rho_under;
            fm->rho[h * 3 + 1] = rho_under;
            fm->rho[h * 3 + 2] = rho_under;
        }
    }

    idx = 0;
    counter = 0;
    while (anomalies--)
    {
        for (L = z_d[idx]; L <= z_d[idx + 1]; L++)
        {
            for (M = y_d[idx]; M <= y_d[idx + 1]; M++)
            {
                for (N = x_d[idx]; N <= x_d[idx + 1]; N++)
                {
                    h = (L - 1) * x * y + (M - 1) * x + N - 1;
                    fm->rho[h * 3 + 0] = rho_d[counter + 0],
                                    fm->rho[h * 3 + 1] = rho_d[counter + 1],
                                    fm->rho[h * 3 + 2] = rho_d[counter + 2];
                    fm->alpha_S[h] = getRadian(alpha_d[counter + 0]),
                    fm->alpha_D[h] = getRadian(alpha_d[counter + 1]),
                    fm->alpha_L[h] = getRadian(alpha_d[counter + 2]);
                    fm->ms[h * 3 + 0] = ms_d[counter + 0],
                                   fm->ms[h * 3 + 1] = ms_d[counter + 1],
                                   fm->ms[h * 3 + 2] = ms_d[counter + 2];
                    fm->ms_S[h] = getRadian(alpha_ms_d[counter + 0]),
                    fm->ms_D[h] = getRadian(alpha_ms_d[counter + 1]),
                    fm->ms_L[h] = getRadian(alpha_ms_d[counter + 2]);
                }
            }
        }
        idx += 2;
        counter += 3;
    }

    for (L = 0; L < z; L++)
        for (M = 0; M < y; M++)
            for (N = 0; N < x; N++)
            {
                // 索引中间变量
                IL = L * (x * (y + 1) + y * (x + 1) + (x + 1) * (y + 1));
                col_index = L * x * y + M * x + N;
                // 为ME矩阵的每一行元素设置索引值
                fm->ME[0 * fm->NE + col_index] = IL + M * (2 * x + 1) + N;
                fm->ME[1 * fm->NE + col_index] = IL + M * (2 * x + 1) + N + 2 * x + 1;
                fm->ME[2 * fm->NE + col_index] = IL + M * (2 * x + 1) + N + x * (y + 1) + y * (x + 1) + (x + 1) * (y + 1);
                fm->ME[3 * fm->NE + col_index] = IL + M * (2 * x + 1) + N + x * (y + 1) + y * (x + 1) + (x + 1) * (y + 1) + 2 * x + 1;
                fm->ME[4 * fm->NE + col_index] = IL + M * (2 * x + 1) + N + x;
                fm->ME[5 * fm->NE + col_index] = IL + M * (2 * x + 1) + N + x + x * (y + 1) + y * (x + 1) + (x + 1) * (y + 1);
                fm->ME[6 * fm->NE + col_index] = IL + M * (2 * x + 1) + N + x + 1;
                fm->ME[7 * fm->NE + col_index] = IL + M * (2 * x + 1) + N + x + 1 + x * (y + 1) + y * (x + 1) + (x + 1) * (y + 1);
                fm->ME[8 * fm->NE + col_index] = IL + x * (y + 1) + y * (x + 1) + M * (x + 1) + N;
                fm->ME[9 * fm->NE + col_index] = IL + x * (y + 1) + y * (x + 1) + M * (x + 1) + N + 1;
                fm->ME[10 * fm->NE + col_index] = IL + x * (y + 1) + y * (x + 1) + M * (x + 1) + N + 1 + x;
                fm->ME[11 * fm->NE + col_index] = IL + x * (y + 1) + y * (x + 1) + M * (x + 1) + N + 2 + x;
            }

    fm->A_X = (double *)malloc(x * sizeof(double));
    fm->B_Y = (double *)malloc(y * sizeof(double));
    fm->C_Z = (double *)malloc(z * sizeof(double));

    idx = 0;
    file = readModelFile(A_X, freq, 0, fm);
    while (fgets(line, sizeof(line), file) != NULL)
    {
        fm->A_X[idx] = strtod(line, NULL);
        idx++;
    }
    fclose(file);

    idx = 0;
    file = readModelFile(B_Y, freq, 0, fm);
    while (fgets(line, sizeof(line), file) != NULL)
    {
        fm->B_Y[idx] = strtod(line, NULL);
        idx++;
    }
    fclose(file);

    idx = 0;
    file = readModelFile(C_Z, freq, 0, fm);
    while (fgets(line, sizeof(line), file) != NULL)
    {
        fm->C_Z[idx] = strtod(line, NULL);
        idx++;
    }
    fclose(file);

    free(z_d);
    free(y_d);
    free(x_d);
    free(rho_d);
    free(ms_d);
    free(alpha_d);
    free(alpha_ms_d);
}

void freeModel(Fmodel *fm)
{
    free(fm->A_X);
    free(fm->B_Y);
    free(fm->C_Z);
    free(fm->alpha_S);
    free(fm->alpha_D);
    free(fm->alpha_L);
    free(fm->ms_S);
    free(fm->ms_D);
    free(fm->ms_L);
    free(fm->rho);
    free(fm->ME);
    free(fm->EtoN);
}