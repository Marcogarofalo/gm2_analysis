#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cmath>     /* erf */ 
#include <ranges>

#include "global.hpp"
#include "resampling.hpp"
#include "read.hpp"
// #include "m_eff.hpp"
// #include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
#include "functions_gm2_analysis.hpp"
// #include "correlators_analysis.hpp"
// #include "eigensystem.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
#include "fit_all.hpp"
#include "resampling_new.hpp"
#include "global.hpp"

#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>
#include <vector>
#include <map>


enum enum_ensembles {
    B72_64,
    B72_96,
    C06,
    D54,
    A53,
    A40,
    A30,
    E112,
    C112
};
constexpr double Mpi_MeV = 135;
constexpr double Mpi_MeV_err = 0.2;

generic_header read_header(FILE* stream) {
    generic_header header;
    int ir = 0;
    ir += fread(&header.T, sizeof(int), 1, stream);
    ir += fread(&header.L, sizeof(int), 1, stream);
    int s;
    ir += fread(&s, sizeof(int), 1, stream);
    header.mus = std::vector<double>(s);
    for (int i = 0; i < s; i++) {
        ir += fread(&header.mus[i], sizeof(double), 1, stream);
    }
    ir += fread(&s, sizeof(int), 1, stream);
    header.thetas = std::vector<double>(s);
    for (int i = 0; i < s; i++) {
        ir += fread(&header.thetas[i], sizeof(double), 1, stream);
    }

    ir += fread(&header.Njack, sizeof(int), 1, stream);
    header.struct_size = ftell(stream);
    return header;


}


double read_single_Nobs(FILE* stream, int header_size, int Njack) {
    int Nobs;
    long int tmp;
    int s = header_size;

    // size_t i = fread(&Njack, sizeof(int), 1, stream);


    fseek(stream, 0, SEEK_END);
    tmp = ftell(stream);
    tmp -= header_size;

    s = Njack;

    Nobs = (tmp) / ((s) * sizeof(double));

    fseek(stream, header_size, SEEK_SET);

    return Nobs;

}

data_single read_single_dataj(FILE* stream) {

    int Njack;
    int Nobs;

    //read_single_Njack_Nobs(stream, header.header_size, Njack, Nobs);
    // data_single dj(Nobs,Njack);
    data_single dj;
    dj.header = read_header(stream);
    dj.Nobs = read_single_Nobs(stream, dj.header.struct_size, dj.header.Njack);
    dj.Njack = dj.header.Njack;
    dj.jack = double_malloc_2(dj.Nobs, dj.Njack);

    //
    size_t i = 0;
    for (int obs = 0; obs < dj.Nobs; obs++) {
        i += fread(dj.jack[obs], sizeof(double), dj.Njack, stream);
    }
    return dj;

}

data_all read_all_the_files(std::vector<std::string> files, const char* resampling) {
    data_all jackall;
    jackall.resampling = resampling;
    //jackall->en = (data_single*)malloc(sizeof(data_single) * files.size());
    jackall.en = new data_single[files.size()];
    jackall.ens = files.size();
    int count = 0;
    for (std::string s : files) {
        FILE* f = open_file(s.c_str(), "r");

        // read_single_dataj(f, params, &(jackall->en[count]));
        jackall.en[count] = read_single_dataj(f);
        jackall.en[count].resampling = resampling;
        count++;
        fclose(f);
    }
    return jackall;

}

double lhs_sum(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    // if (n == 0)        r = gjack.en[e].jack[fit_info.corr_id[0]][j] + gjack.en[e].jack[fit_info.corr_id[2]][j] + gjack.en[e].jack[fit_info.corr_id[4]][j];
    // else if (n == 1)   r = gjack.en[e].jack[fit_info.corr_id[1]][j] + gjack.en[e].jack[fit_info.corr_id[3]][j] + gjack.en[e].jack[fit_info.corr_id[5]][j];
    // r = gjack.en[e].jack[fit_info.corr_id[0 + n]][j] + gjack.en[e].jack[fit_info.corr_id[2 + n]][j] + gjack.en[e].jack[fit_info.corr_id[4 + n]][j];
    r = 0;
    for (int i = 0;i < fit_info.corr_id.size();i += 2)
        r += gjack.en[e].jack[fit_info.corr_id[i + n]][j];
    return r;
}

void sum_lsc(data_all in, const char* outpath, const char* filename) {
    int N = in.Nfits;
    int Njack = in.fits[0].Njack;
    char name[NAMESIZE];
    mysprintf(name, NAMESIZE, "%s/%s", outpath, filename);
    mysprintf(name, NAMESIZE, "%s/%s", outpath, filename);
    FILE* f = open_file(name, "w+");
    printf("writing: %s\n", name);
    double* ave = (double*)calloc(Njack, sizeof(double*));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < Njack; j++) {
            ave[j] += in.fits[i].P[0][j];
        }
        fprintf(f, "%s %g  %g\n", in.fits[i].name, in.fits[i].P[0][Njack - 1], error_jackboot(in.resampling.c_str(), Njack, in.fits[i].P[0]));
    }
    fprintf(f, "total: %g  %g\n", ave[Njack - 1], error_jackboot(in.resampling.c_str(), Njack, ave));

    double** y = double_malloc_2(N, Njack);
    for (int j = 0; j < Njack; j++) {
        for (int i = 0; i < N; i++) {
            y[i][j] = in.fits[i].P[0][j];
        }
    }
    double** cov = covariance(in.resampling.c_str(), N, Njack, y);
    fprintf(f, "\n#covariance:\n");
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            fprintf(f, "%g\t", cov[i][k]);
        }
        fprintf(f, "\n");
    }
    fprintf(f, "\n#correlation:\n");
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            fprintf(f, "%g\t", cov[i][k] / sqrt(cov[i][i] * cov[k][k]));
        }
        fprintf(f, "\n");
    }
    fclose(f);

}

double comp_pool_variable(double* j1, double* j2) {
    double err1 = myres->comp_error(j1);
    double err2 = myres->comp_error(j2);
    int  Nj = myres->Njack;
    return (j1[Nj - 1] - j2[Nj - 1]) / sqrt(err1 * err1 + err2 * err2);
}
double comp_error_pool(double* j1, double* j2) {

    int  Nj = myres->Njack;
    double P = comp_pool_variable(j1, j2);

    return fabs(j1[Nj - 1] - j2[Nj - 1]) * erf(fabs(P) / sqrt(2.0));
}



double comp_error_pool_func(data_all gjack, fit_type fit_info, double (*lhs_fun)(int, int, int, data_all, struct fit_type)) {

    int  Nj = myres->Njack;
    double* j1 = (double*)malloc(sizeof(double) * Nj);
    double* j2 = (double*)malloc(sizeof(double) * Nj);
    for (int j = 0;j < Nj;j++) {
        j1[j] = lhs_fun(0, fit_info.myen[0], j, gjack, fit_info);
        j2[j] = lhs_fun(0, fit_info.myen[1], j, gjack, fit_info);
    }

    return comp_error_pool(j1, j2);
}



int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");
    char namefile[NAMESIZE];


    std::vector<std::string> files;
    mysprintf(namefile, NAMESIZE, "%s/%s_cB.72.64_mu.0.000720", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB.72.96_mu.0.000720", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC.06.80_mu.0.000600", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD.54.96_mu.0.000540", argv[2], argv[1]);
    files.emplace_back(namefile);

    std::vector<int> myen(files.size());
    for (int e = 0; e < files.size(); e++) {
        myen[e] = e;
        printf("file:  %s\n", files[e].c_str());
        //printf("Nobs=%d  Njack=%d    mus=", jackall.en[e].Nobs, jackall.en[e].Njack);
        //for (double mu : jackall.en[e].header.mus)
        //    printf("%g   ", mu);
        //printf("\n");
    }

    mysprintf(namefile, NAMESIZE, "%s/%s_cA.53.24_mu.0.005300", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cA.40.24_mu.0.004000", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cA.30.32_mu.0.003000", argv[2], argv[1]);
    files.emplace_back(namefile);

    mysprintf(namefile, NAMESIZE, "%s/%s_cE.44.112_mu.0.000440", argv[2], argv[1]);
    files.emplace_back(namefile);

    mysprintf(namefile, NAMESIZE, "%s/%s_cC.06.112_mu.0.000600", argv[2], argv[1]);
    files.emplace_back(namefile);


    data_all jackall = read_all_the_files(files, argv[1]);
    jackall.create_generalised_resampling();


    std::vector<int> myen_full(jackall.ens);
    for (int e = 0; e < jackall.ens; e++) {
        myen_full[e] = e;
        printf("file:  %s\n", files[e].c_str());
        printf("Nobs=%d  Njack=%d    mus=", jackall.en[e].Nobs, jackall.en[e].Njack);
        for (double mu : jackall.en[e].header.mus)
            printf("%g   ", mu);
        printf("\n");
    }
    int Njack = jackall.en[0].Njack;
    std::vector<int> myen_charm = { 0, 2, 3, 4, 5, 6 };

    if (strcmp(argv[1], "jack") == 0) {
        myres = new resampling_jack(Njack - 1);
    }
    else if (strcmp(argv[1], "boot") == 0) {
        myres = new resampling_boot(Njack - 1);
    }

    double* jack_Mpi_MeV_exp = fake_sampling(argv[1], Mpi_MeV, Mpi_MeV_err, Njack, 1003);
    //////////////////////////////////////////////////////////////
    // new jackall
    //////////////////////////////////////////////////////////////
    data_all jackextra;
    jackextra.resampling = jackall.resampling;

    jackextra.en = new data_single[jackall.ens];
    jackextra.ens = jackall.ens;
    for (int count = 0; count < jackextra.ens; count++) {

        data_single dj;
        dj.header = jackall.en[count].header;
        dj.Nobs = jackall.en[count].Nobs + 1 + 50;
        dj.Njack = dj.header.Njack;
        dj.jack = double_malloc_2(dj.Nobs, dj.Njack);

        //
        size_t i = 0;
        for (int obs = 0; obs < jackall.en[count].Nobs; obs++) {

            for (int j = 0;j < Njack;j++)
                dj.jack[obs][j] = jackall.en[count].jack[obs][j];
        }
        dj.resampling = jackall.en[count].resampling;

        jackextra.en[count] = dj;

    }



    int obs = jackall.en[0].Nobs;
    std::vector<int> id_SD = { 167, 168 };
    std::vector<int> id_W = { 169, 170 };
    std::vector<int> id_LD = { 175, 176 };
    std::vector<int> id_full = { 146, 147 };
    std::vector<int> id_fulltree = { 181, 182 };
    std::vector<int> id_SDtmin0 = { 211, 212 };
    std::vector<int> id_SDtmin1 = { 213, 214 };
    std::vector<int> id_SDtmin2 = { 215, 216 };
    std::vector<int> id_SDtmin3 = { 217, 218 };
    std::vector<int> id_SDtmin4 = { 219, 220 };
    std::vector<int> id_SDeta = { 37, 40 };
    std::vector<int> id_Weta = { 54, 57 };
    std::vector<int> id_LDeta = { 232, 233 };



    std::vector<int> id_SD_cor = { obs + 1, obs + 2 };
    std::vector<int> id_W_cor = { obs + 3, obs + 4 };
    std::vector<int> id_LD_cor = { obs + 5, obs + 6 };
    std::vector<int> id_full_cor = { obs + 7, obs + 8 };

    std::vector<int> id_SD_FVE = { obs + 9, obs + 10 };
    std::vector<int> id_W_FVE = { obs + 11, obs + 12 };
    std::vector<int> id_LD_FVE = { obs + 13, obs + 14 };
    std::vector<int> id_full_FVE = { obs + 15, obs + 16 };

    std::vector<int> ensemble_to_correct = { B72_64, B72_96, C06, C112 ,D54, E112 };
    std::vector<double*> damu_SD(files.size());
    std::vector<double*> damu_W(files.size());
    std::vector<double*> damu_LD(files.size());
    std::vector<double*> damu_full(files.size());

    double scale = 1e-10;
    for (int e = 0; e < files.size(); e++) {
        if (e == B72_64 || e == B72_96) {
            damu_SD[e] = myres->create_fake(4.4159e-05 * scale, 0.000392418 * scale, 2000 + B72_64);
            damu_W[e] = myres->create_fake(-0.00688788 * scale, 0.00925944 * scale, 2000 + B72_64 + files.size());
            damu_full[e] = myres->create_fake(-0.0327914 * scale, 0.0385565 * scale, 2000 + B72_64 + 2 * files.size());
            damu_LD[e] = myres->create_fake(-0.0259477 * scale, 0.0297455 * scale, 2000 + B72_64 + 3 * files.size());
        }
        else if (e == C06 || e == C112) {
            damu_SD[e] = myres->create_fake(0.0008153 * scale, 0.0013851 * scale, 2000 + C06);
            damu_W[e] = myres->create_fake(0.0284585 * scale, 0.0159039 * scale, 2000 + C06 + files.size());
            damu_full[e] = myres->create_fake(0.09633 * scale, 0.0580504 * scale, 2000 + C06 + 2 * files.size()); // FULl HVP and then subtract
            damu_LD[e] = myres->create_fake(0.0670562 * scale, 0.0436722 * scale, 2000 + C06 + 3 * files.size());
        }
        else if (e == D54) {
            damu_SD[e] = myres->create_fake(0.000635405 * scale, 0.000967111 * scale, 2000 + D54);
            damu_W[e] = myres->create_fake(0.0104084 * scale, 0.0315747 * scale, 2000 + D54 + files.size());
            damu_full[e] = myres->create_fake(0.0321713 * scale, 0.13259 * scale, 2000 + D54 + 2 * files.size()); // FULl HVP and then subtract
            damu_LD[e] = myres->create_fake(0.0211275 * scale, 0.101898 * scale, 1);
        }
        else if (e == E112) {
            damu_SD[e] = myres->create_fake(0.00230056 * scale, 0.00218935 * scale, 2000 + E112);
            damu_W[e] = myres->create_fake(0.0196878 * scale, 0.0205063 * scale, 2000 + E112 + files.size());
            damu_full[e] = myres->create_fake(0.0893201 * scale, 0.0523576 * scale, 2000 + E112 + 2 * files.size()); // FULl HVP and then subtract
            damu_LD[e] = myres->create_fake(0.0673317 * scale, 0.0338949 * scale, 2000 + E112 + 3 * files.size());
        }
        else {
            damu_SD[e] = myres->create_fake(0.0, 1e-16, 2000 + D54);
            damu_W[e] = myres->create_fake(0.0, 1e-16, 2000 + D54 + files.size());
            damu_full[e] = myres->create_fake(0 * scale, 1e-16 * scale, 2000 + E112 + 2 * files.size()); // FULl HVP and then subtract
            damu_LD[e] = myres->create_fake(0.0, 1e-16, 2000 + D54 + 2 * files.size()); // FULl HVP and then subtract

        }

    }

    for (int e : ensemble_to_correct) {
        for (int j = 0;j < Njack;j++) {
            for (int tm = 0;tm < 2;tm++) {
                jackextra.en[e].jack[id_SD_cor[tm]][j] = damu_SD[e][j];
                jackextra.en[e].jack[id_W_cor[tm]][j] = damu_W[e][j];
                jackextra.en[e].jack[id_LD_cor[tm]][j] = damu_LD[e][j];
                jackextra.en[e].jack[id_full_cor[tm]][j] = damu_full[e][j];


                jackextra.en[e].jack[id_SD_FVE[tm]][j] = jackall.en[e].jack[id_SDeta[tm]][j];
                jackextra.en[e].jack[id_W_FVE[tm]][j] = jackall.en[e].jack[id_Weta[tm]][j];
                jackextra.en[e].jack[id_LD_FVE[tm]][j] = jackall.en[e].jack[id_LDeta[tm]][j];
                // jackextra.en[e].jack[id_full_FVE[tm]][j] = jackall.en[e].jack[id_fulleta[tm]][j];

            }
        }
    }
    {
        int e = B72_64;
        int eadd = B72_96;
        for (int j = 0;j < Njack;j++) {
            for (int tm = 0;tm < 2;tm++) {
                jackextra.en[e].jack[id_SD_FVE[tm]][j] += jackextra.en[eadd].jack[id_SD_FVE[tm]][j];
                jackextra.en[e].jack[id_W_FVE[tm]][j] += jackextra.en[eadd].jack[id_W_FVE[tm]][j];
                jackextra.en[e].jack[id_LD_FVE[tm]][j] += jackextra.en[eadd].jack[id_LD_FVE[tm]][j];
                // jackextra.en[e].jack[id_full_FVE[tm]][j] += jackextra.en[eadd].jack[id_full_FVE[tm]][j];
                jackextra.en[e].jack[id_SD_FVE[tm]][j] /= 2.0;
                jackextra.en[e].jack[id_W_FVE[tm]][j] /= 2.0;
                jackextra.en[e].jack[id_LD_FVE[tm]][j] /= 2.0;
                // jackextra.en[e].jack[id_full_FVE[tm]][j] /= 2.0;
            }
        }
        e = C06;
        eadd = C112;
        for (int j = 0;j < Njack;j++) {
            for (int tm = 0;tm < 2;tm++) {
                jackextra.en[e].jack[id_SD_FVE[tm]][j] += jackextra.en[eadd].jack[id_SD_FVE[tm]][j];
                jackextra.en[e].jack[id_W_FVE[tm]][j] += jackextra.en[eadd].jack[id_W_FVE[tm]][j];
                jackextra.en[e].jack[id_LD_FVE[tm]][j] += jackextra.en[eadd].jack[id_LD_FVE[tm]][j];
                // jackextra.en[e].jack[id_full_FVE[tm]][j] += jackextra.en[eadd].jack[id_full_FVE[tm]][j];
                jackextra.en[e].jack[id_SD_FVE[tm]][j] /= 2.0;
                jackextra.en[e].jack[id_W_FVE[tm]][j] /= 2.0;
                jackextra.en[e].jack[id_LD_FVE[tm]][j] /= 2.0;
                // jackextra.en[e].jack[id_full_FVE[tm]][j] /= 2.0;
            }
        }
    }



    std::vector<std::string>   interpolations;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // fits 
    /////////////////////////////////////////////////////////////////////////////////////////////////
    fit_type fit_info;
    data_all  syst_amu_SD_l;
    syst_amu_SD_l.resampling = argv[1];
    data_all  sum_amu_SD;
    sum_amu_SD.resampling = argv[1];
    data_all  sum_amu_W;
    sum_amu_W.resampling = argv[1];

    int count = 0;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.restore_default();
    std::string namefit;


    for (int iW = 0;iW < 28;iW++) {
        for (int ie = 0;ie < 14;ie++) {

            std::vector<int> fi_list;
            if (ie < 7)
                fi_list = { 0,1,2,3,4,5,6,   10,11,12,13,14,15 };
            else if (ie >= 7 && ie < 13) // only one regularization
                fi_list = { 7,8,9, 16,17 };
            else if (ie > 12) // with FVE
                fi_list = { 0,1,2,3,4,5,6,   10,11,12,13,14,15 };



            for (int fi : fi_list) {

                namefit = "amu";

                switch (iW) {
                case 0:
                    namefit = namefit + "_SD";
                    fit_info.corr_id = { 167, 168 };
                    break;
                case 1:
                    namefit = namefit + "_W";
                    fit_info.corr_id = { 169, 170 };
                    break;
                case 2:
                    namefit = namefit + "_LD";
                    fit_info.corr_id = { 175, 176 };
                    break;
                case 3:
                    namefit = namefit + "_full";
                    fit_info.corr_id = { 146, 147 };
                    break;
                case 4:
                    namefit = namefit + "_fulltree";
                    fit_info.corr_id = { 181, 182 };
                    break;
                case 5:
                    namefit = namefit + "_SDpWpLD";
                    fit_info.corr_id = { 167, 168,169, 170, 175, 176 };
                    break;
                case 6:
                    namefit = namefit + "_SDtmin0";
                    fit_info.corr_id = { 211, 212 , id_SD_cor[0], id_SD_cor[1] }; // SD tmin 0
                    break;
                case 7:
                    namefit = namefit + "_SDtmin1";
                    fit_info.corr_id = { 213, 214 , id_SD_cor[0], id_SD_cor[1] }; // SD tmin 1
                    break;
                case 8:
                    namefit = namefit + "_SDtmin2";
                    fit_info.corr_id = { 215, 216 , id_SD_cor[0], id_SD_cor[1] }; // SD tmin 2
                    break;
                case 9:
                    namefit = namefit + "_SDtmin3";
                    fit_info.corr_id = { 217, 218 , id_SD_cor[0], id_SD_cor[1] }; // SD tmin 3
                    break;
                case 10:
                    namefit = namefit + "_SDtmin4";
                    fit_info.corr_id = { 219, 220 , id_SD_cor[0], id_SD_cor[1] }; // SD tmin 4
                    break;
                    //////////////////  corrected
                case 11:
                    namefit = namefit + "_SDcor";
                    fit_info.corr_id = { id_SD[0], id_SD[1],id_SD_cor[0], id_SD_cor[1] };
                    break;
                case 12:
                    namefit = namefit + "_Wcor";
                    fit_info.corr_id = { id_W[0], id_W[1], id_W_cor[0], id_W_cor[1] };
                    break;
                case 13:
                    namefit = namefit + "_LDcor";
                    fit_info.corr_id = { id_LD[0], id_LD[1], id_LD_cor[0], id_LD_cor[1] };
                    break;
                case 14:
                    namefit = namefit + "_SDpWpLDcor";
                    fit_info.corr_id = { id_SD[0], id_SD[1],id_W[0], id_W[1], id_LD[0], id_LD[1], id_full_cor[0], id_full_cor[1] };
                    // fit_info.corr_id = { id_SD[0], id_SD[1], id_SD_cor[0], id_SD_cor[1],id_W[0], id_W[1], id_W_cor[0], id_W_cor[1], id_LD[0], id_LD[1], id_LD_cor[0], id_LD_cor[1] };
                    // fit_info.corr_id = {id_LD_cor[0], id_LD_cor[1]      , id_LD[0], id_LD[1], id_LD[0], id_LD[1], id_W_cor[0], id_W_cor[1], id_SD_cor[1],id_W[0], id_W[1],
                    //  id_SD[0], id_SD[1], id_SD_cor[0]  };
                    break;
                case 15:
                    namefit = namefit + "_SDetas";
                    fit_info.corr_id = { id_SDeta[0], id_SDeta[1], id_SD_cor[0], id_SD_cor[1] };
                    break;
                case 16:
                    namefit = namefit + "_Wetas";
                    fit_info.corr_id = { id_Weta[0], id_Weta[1], id_W_cor[0], id_W_cor[1] };
                    break;
                case 17:
                    namefit = namefit + "_LDetas";
                    fit_info.corr_id = { id_LDeta[0], id_LDeta[1], id_LD_cor[0], id_LD_cor[1] };
                    break;
                case 18:
                    namefit = namefit + "_SDpWpLDetas";
                    fit_info.corr_id = { id_SDeta[0], id_SDeta[1],id_Weta[0], id_Weta[1], id_LDeta[0], id_LDeta[1] ,id_full_cor[0], id_full_cor[1] };
                    break;
                case 19:
                    namefit = namefit + "_SDetasFVE";
                    fit_info.corr_id = { id_SD_FVE[0], id_SD_FVE[1], id_SD_cor[0], id_SD_cor[1] };
                    break;
                case 20:
                    namefit = namefit + "_WetasFVE";
                    fit_info.corr_id = { id_W_FVE[0], id_W_FVE[1], id_W_cor[0], id_W_cor[1] };
                    break;
                case 21:
                    namefit = namefit + "_LDetasFVE";
                    fit_info.corr_id = { id_LD_FVE[0], id_LD_FVE[1], id_LD_cor[0], id_LD_cor[1] };
                    break;
                case 22:
                    namefit = namefit + "_SDpWpLDetasFVE";
                    fit_info.corr_id = { id_SD_FVE[0], id_SD_FVE[1],id_W_FVE[0], id_W_FVE[1], id_LD_FVE[0], id_LD_FVE[1] ,id_full_cor[0], id_full_cor[1] };
                    break;
                case 23:
                    namefit = namefit + "_SDetasNoCor";
                    fit_info.corr_id = { id_SDeta[0], id_SDeta[1] };
                    break;
                case 24:
                    namefit = namefit + "_WetasNoCor";
                    fit_info.corr_id = { id_Weta[0], id_Weta[1] };
                    break;
                case 25:
                    namefit = namefit + "_LDetasNoCor";
                    fit_info.corr_id = { id_LDeta[0], id_LDeta[1] };
                    break;
                case 26:
                    namefit = namefit + "_SDpWpLDetasNoCor";
                    fit_info.corr_id = { id_SDeta[0], id_SDeta[1],id_Weta[0], id_Weta[1], id_LDeta[0], id_LDeta[1] };
                    break;
                case 27:
                    namefit = namefit + "_SDetasFVENoCor";
                    fit_info.corr_id = { id_SD_FVE[0], id_SD_FVE[1] };
                    break;
                    // case 15:
                    //     namefit = namefit + "_SDtmin0cor";
                    //     fit_info.corr_id = { id_SDtmin0_cor[0], id_SDtmin0_cor[1] }; // SD tmin 0
                    //     break;
                    // case 16:
                    //     namefit = namefit + "_SDtmin1cor";
                    //     fit_info.corr_id = { id_SDtmin1_cor[0], id_SDtmin1_cor[1] };
                    //     break;
                    // case 17:
                    //     namefit = namefit + "_SDtmin2cor";
                    //     fit_info.corr_id = { id_SDtmin2_cor[0], id_SDtmin2_cor[1] };
                    //     break;
                    // case 18:
                    //     namefit = namefit + "_SDtmin3cor";
                    //     fit_info.corr_id = { id_SDtmin3_cor[0], id_SDtmin3_cor[1] };
                    //     break;
                    // case 19:
                    //     namefit = namefit + "_SDtmin4cor";
                    //     fit_info.corr_id = { id_SDtmin4_cor[0], id_SDtmin4_cor[1] };
                    //     break;
                default: break;
                }


                // if fititng only the TM we need to put the id only in the even slots so that 
                // the function lhs_sum sum them when n=0
                if (ie == 8 || ie == 10 || ie == 12) {
                    for (int i = 0;i < fit_info.corr_id.size() / 2;i++)
                        fit_info.corr_id[i * 2] = fit_info.corr_id[i * 2 + 1];
                }

                // if fitting the sum SD+W+LD we need to select lhs function to sum
                double (*lhs_fun)(int, int, int, data_all, struct fit_type) = lhs_sum;
                // if (iW == 5 || iW == 14 || iW == 18) lhs_fun = lhs_sum;
                // else lhs_fun = lhs_amu;


                switch (ie) {
                case 0:
                    namefit = namefit + "_3b";
                    fit_info.Nxen = { {  C06 ,D54, E112},
                                      { C06, D54, E112} };
                    break;
                case 1:
                    namefit = namefit + "_3b_BOS";
                    fit_info.Nxen = { { B72_64, C06 ,D54, E112},
                                      { C06, D54, E112} };
                    break;
                case 2:
                    namefit = namefit + "_3b_BTM";
                    fit_info.Nxen = { {  C06 ,D54, E112},
                                      {B72_64, C06, D54, E112} };
                    break;
                case 3:
                    namefit = namefit + "_3b_BOS_BTM";
                    fit_info.Nxen = { { B72_64, C06 ,D54, E112},
                                      {B72_64, C06, D54, E112} };
                    break;
                case 4:
                    namefit = namefit + "_3b_noC";
                    fit_info.Nxen = { {  B72_64 ,D54, E112},
                                      { B72_64, D54, E112} };
                    break;
                case 5:
                    namefit = namefit + "_3b_noC_BOS";
                    fit_info.Nxen = { { B72_64, C06 ,D54, E112},
                                      { B72_64, D54, E112} };
                    break;
                case 6:
                    namefit = namefit + "_3b_noC_BTM";
                    fit_info.Nxen = { {  B72_64 ,D54, E112},
                                      {B72_64, C06, D54, E112} };
                    break;
                case 7:
                    namefit = namefit + "_3b_onlyOS";
                    fit_info.Nxen = { {  C06 ,D54, E112} };
                    break;
                case 8:
                    namefit = namefit + "_3b_onlyTM";
                    fit_info.Nxen = { {  C06 ,D54, E112} };
                    break;
                case 9:
                    namefit = namefit + "_4b_onlyOS";
                    fit_info.Nxen = { { B72_64, C06 ,D54, E112} };
                    break;
                case 10:
                    namefit = namefit + "_4b_onlyTM";
                    fit_info.Nxen = { { B72_64, C06 ,D54, E112} };
                    break;
                case 11:
                    namefit = namefit + "_3b_noC_onlyOS";
                    fit_info.Nxen = { {  B72_64 ,D54, E112} };
                    break;
                case 12:
                    namefit = namefit + "_3b_noC_onlyTM";
                    fit_info.Nxen = { {  B72_64 ,D54, E112} };
                    break;
                case 13:
                    namefit = namefit + "_3b_BOS_BTM_FVE";
                    fit_info.Nxen = { { B72_64, B72_96, C06, C112 ,D54, E112},
                                      { B72_64, B72_96, C06, C112, D54, E112} };
                    break;
                default:
                    break;
                }
                fit_info.N = fit_info.Nxen.size();


                switch (fi) {
                case 0:
                    namefit = namefit + "";
                    fit_info.Npar = 3;
                    fit_info.function = rhs_amu_common;
                    break;
                case 1:
                    namefit = namefit + "_a4OS";
                    fit_info.Npar = 4;
                    fit_info.function = rhs_amu_a4OS_common;
                    break;
                case 2:
                    namefit = namefit + "_a4TM";
                    fit_info.Npar = 4;
                    fit_info.function = rhs_amu_a4TM_common;
                    break;
                case 3:
                    namefit = namefit + "_a4OS_a4TM";
                    fit_info.Npar = 5;
                    fit_info.function = rhs_amu_a4OS_a4TM_common;
                    break;
                case 4:
                    namefit = namefit + "_alogOS";
                    fit_info.Npar = 4;
                    fit_info.function = rhs_amu_alogOS_common;
                    break;
                case 5:
                    namefit = namefit + "_alogTM";
                    fit_info.Npar = 4;
                    fit_info.function = rhs_amu_alogTM_common;
                    break;
                case 6:
                    namefit = namefit + "_alogOS_alogTM";
                    fit_info.Npar = 5;
                    fit_info.function = rhs_amu_alogOS_alogTM_common;
                    break;
                case 7:
                    namefit = namefit + "";
                    fit_info.Npar = 2;
                    fit_info.function = rhs_amu_onlyOSTM;
                    break;
                case 8:
                    namefit = namefit + "_a4";
                    fit_info.Npar = 3;
                    fit_info.function = rhs_amu_a4_onlyOSTM;
                    break;
                case 9:
                    namefit = namefit + "_alog";
                    fit_info.Npar = 3;
                    fit_info.function = rhs_amu_alog_onlyOSTM;
                    break;
                case 10:
                    namefit = namefit + "_alog2OS";
                    fit_info.Npar = 4;
                    fit_info.function = rhs_amu_alog2OS_common;
                    break;
                case 11:
                    namefit = namefit + "_alog2TM";
                    fit_info.Npar = 4;
                    fit_info.function = rhs_amu_alog2TM_common;
                    break;
                case 12:
                    namefit = namefit + "_alog2OS_alog2TM";
                    fit_info.Npar = 5;
                    fit_info.function = rhs_amu_alog2OS_alog2TM_common;
                    break;
                case 13:
                    namefit = namefit + "_alog3OS";
                    fit_info.Npar = 4;
                    fit_info.function = rhs_amu_alog3OS_common;
                    break;
                case 14:
                    namefit = namefit + "_alog3TM";
                    fit_info.Npar = 4;
                    fit_info.function = rhs_amu_alog3TM_common;
                    break;
                case 15:
                    namefit = namefit + "_alog3OS_alog3TM";
                    fit_info.Npar = 5;
                    fit_info.function = rhs_amu_alog3OS_alog3TM_common;
                    break;
                case 16:
                    namefit = namefit + "_alog2";
                    fit_info.Npar = 3;
                    fit_info.function = rhs_amu_alog2_onlyOSTM;
                    break;
                case 17:
                    namefit = namefit + "_alog3";
                    fit_info.Npar = 3;
                    fit_info.function = rhs_amu_alog3_onlyOSTM;
                    break;
                default:
                    break;
                }




                fit_info.Nvar = 8;
                fit_info.Njack = Njack;
                fit_info.init_N_etot_form_Nxen();
                fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);

                if (fit_info.entot <= fit_info.Npar) continue;

                count = 0;
                for (int n = 0;n < fit_info.N;n++) {
                    for (int e : fit_info.Nxen[n]) {
                        for (int j = 0;j < Njack;j++) {
                            fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
                            // fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
                            fit_info.x[1][count][j] = jackall.en[e].jack[58][j];  // Delta_FV_GS
                            fit_info.x[2][count][j] = jackall.en[e].jack[1][j];  //Mpi
                            fit_info.x[3][count][j] = jack_Mpi_MeV_exp[j];
                            fit_info.x[4][count][j] = 0/* l */ + 1e-6;
                            fit_info.x[5][count][j] = 0/* a */ + 1e-6;
                            fit_info.x[6][count][j] = 0 + 1e-6;
                            fit_info.x[7][count][j] = 0/* w */ + 1e-6;
                        }
                        count++;
                    }
                }

                fit_info.linear_fit = true;
                // fit_info.acc= 1e-6;
                // fit_info.chi2_gap_jackboot=0.1;
                // fit_info.guess_per_jack=5;
                // fit_info.repeat_start=5;
                fit_info.verbosity = 0;
                fit_info.covariancey = true;
                fit_info.compute_cov_fit(argv, jackextra, lhs_fun);
                int ide = 0, ide1 = 0;
                for (int n = 0;n < fit_info.Nxen.size();n++) {
                    for (int e : fit_info.Nxen[n]) {
                        ide1 = 0;
                        for (int n1 = 0;n1 < fit_info.Nxen.size();n1++) {
                            for (int e1 : fit_info.Nxen[n1]) {
                                if (e != e1)   fit_info.cov[ide][ide1] = 0;
                                // printf("%-12.5g ", fit_info.cov[ide][ide1]);
                                ide1++;
                            }
                        }
                        // printf("\n");
                        ide++;
                    }
                }
                fit_info.compute_cov1_fit();
                fit_result amu_SD_l_common_a4 = fit_all_data(argv, jackextra, lhs_fun, fit_info, namefit.c_str());
                fit_info.band_range = { 0,0.0081 };
                std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
                     fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };


                //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
                print_fit_band(argv, jackextra, fit_info, fit_info, namefit.c_str(), "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0002, xcont);

                free_fit_result(fit_info, amu_SD_l_common_a4);
                // if (namefit.compare())
            }
        }
    }

    //////////////////////////////////////////////////////////////
    // all ens
    //////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////
    // volume pool variable
    //////////////////////////////////////////////////////////////
    {
        std::vector<double> P(4);
        std::vector<int> id;
        id = id_SDeta;
        printf("%g  %g\n", jackextra.en[B72_64].jack[id[0]][Njack - 1], myres->comp_error(jackextra.en[B72_64].jack[id[0]]));
        printf("%g  %g\n", jackextra.en[B72_96].jack[id[0]][Njack - 1], myres->comp_error(jackextra.en[B72_96].jack[id[0]]));
        P[0] = comp_error_pool(jackextra.en[B72_64].jack[id[0]], jackextra.en[B72_96].jack[id[0]]);
        P[1] = comp_error_pool(jackextra.en[B72_64].jack[id[1]], jackextra.en[B72_96].jack[id[1]]);
        P[2] = comp_error_pool(jackextra.en[C06].jack[id[0]], jackextra.en[C112].jack[id[0]]);
        P[3] = comp_error_pool(jackextra.en[C06].jack[id[1]], jackextra.en[C112].jack[id[1]]);
        auto delta = std::ranges::max_element(P.begin(), P.end());
        for (double p : P) printf("%g\n", p);
        printf("SD volume error FVE to add Delta=%g\n", *delta);
        id = id_Weta;
        P[0] = comp_error_pool(jackextra.en[B72_64].jack[id[0]], jackextra.en[B72_96].jack[id[0]]);
        P[1] = comp_error_pool(jackextra.en[B72_64].jack[id[1]], jackextra.en[B72_96].jack[id[1]]);
        P[2] = comp_error_pool(jackextra.en[C06].jack[id[0]], jackextra.en[C112].jack[id[0]]);
        P[3] = comp_error_pool(jackextra.en[C06].jack[id[1]], jackextra.en[C112].jack[id[1]]);
        delta = std::ranges::max_element(P.begin(), P.end());
        for (double p : P) printf("%g\n", p);
        printf("W volume error FVE to add Delta=%g\n", *delta);
        id = id_LDeta;
        P[0] = comp_error_pool(jackextra.en[B72_64].jack[id[0]], jackextra.en[B72_96].jack[id[0]]);
        P[1] = comp_error_pool(jackextra.en[B72_64].jack[id[1]], jackextra.en[B72_96].jack[id[1]]);
        P[2] = comp_error_pool(jackextra.en[C06].jack[id[0]], jackextra.en[C112].jack[id[0]]);
        P[3] = comp_error_pool(jackextra.en[C06].jack[id[1]], jackextra.en[C112].jack[id[1]]);
        // for (int i = 0;i < P.size();i++) printf("Pi=%g\n", P[i]);
        delta = std::ranges::max_element(P.begin(), P.end());
        for (double p : P) printf("%g\n", p);
        printf("LD volume error FVE to add Delta=%g\n", *delta);

        /// full
        fit_info.myen = { B72_64, B72_96 };
        fit_info.corr_id = { id_SDeta[0], -1, id_Weta[0], -1, id_LDeta[0] };
        P[0] = comp_error_pool_func(jackextra, fit_info, lhs_sum);
        fit_info.corr_id = { id_SDeta[1], -1, id_Weta[1],-1, id_LDeta[1] }; // TM
        P[1] = comp_error_pool_func(jackextra, fit_info, lhs_sum);
        fit_info.myen = { C06,C112 };
        fit_info.corr_id = { id_SDeta[0], -1, id_Weta[0],-1, id_LDeta[0] };
        P[2] = comp_error_pool_func(jackextra, fit_info, lhs_sum);
        fit_info.corr_id = { id_SDeta[1], -1, id_Weta[1],-1, id_LDeta[1] }; // TM
        P[3] = comp_error_pool_func(jackextra, fit_info, lhs_sum);
        delta = std::ranges::max_element(P.begin(), P.end());
        for (double p : P) printf("%g\n", p);
        printf("SDpWpLD volume error FVE to add Delta=%g\n", *delta);


        // id = id_full;
        // P[0] = comp_error_pool(jackextra.en[B72_64].jack[id[0]], jackextra.en[B72_96].jack[id[0]]);
        // P[1] = comp_error_pool(jackextra.en[B72_64].jack[id[1]], jackextra.en[B72_96].jack[id[1]]);
        // P[2] = comp_error_pool(jackextra.en[C06].jack[id[0]], jackextra.en[C112].jack[id[0]]);
        // P[3] = comp_error_pool(jackextra.en[C06].jack[id[1]], jackextra.en[C112].jack[id[1]]);
        // delta = std::ranges::max_element(P.begin(), P.end());
        // printf("SDpWpLD volume error FVE to add Delta=%g\n", *delta);
    }


    /// fit data with a polynomila   data(a^2) = lamb + Poly(a^2)
    // namefit = "amu";
    // namefit = namefit + "_Wcor";
    // fit_info.corr_id = { id_W[0],  id_W_cor[0] };
    // double (*lhs_fun)(int, int, int, data_all, struct fit_type) = lhs_sum;
    // namefit = namefit + "_3b_BOS";
    // fit_info.Nxen = { { B72_64, C06 ,D54, E112}};
}