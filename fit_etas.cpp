#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

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

int get_en_in_lhs_Nxen(int n, int e, struct fit_type fit_info) {
    int count = 0;
    for (int in = 0;in < n;in++) {
        for (int ie : fit_info.Nxen[in]) {
            count++;
        }
    }
    for (int ee : fit_info.Nxen[n]) {
        if (ee == e) {
            break;
        }
        count++;
    }
    return count;
}

double lhs_Meta_MeV(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
    double r;
    int ie = get_en_in_lhs_Nxen(n, e, fit_info);
    r = gjack.en[e].jack[fit_info.corr_id[0]][j] *hbarc/ sqrt(fit_info.x[0][ie][j]);
    return r;
}


double rhs_a2_a4(int n, int Nvar, double* x, int Npar, double* P) {
    double r;
    double a = x[0];
    r = P[0] + a * P[1] + a * a * P[2];
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


    for (int iW = 0;iW < 1;iW++) {
        for (int ie = 0;ie < 2;ie++) {

            std::vector<int> fi_list;
            fi_list = { 0,1 };

            for (int fi : fi_list) {

                namefit = "Metas";

                switch (iW) {
                case 0:
                    namefit = namefit + "_TM";
                    fit_info.corr_id = { 231 };
                    break;
                default:
                    break;
                }



                double (*lhs_fun)(int, int, int, data_all, struct fit_type);
                lhs_fun = lhs_Meta_MeV;


                switch (ie) {
                case 0:
                    namefit = namefit + "_4b";
                    fit_info.Nxen = { { B72_64, C06 ,D54, E112} };
                    break;
                case 1:
                    namefit = namefit + "_3b";
                    fit_info.Nxen = { {  C06 ,D54, E112} };
                    break;
                default:
                    break;
                }
                fit_info.N = fit_info.Nxen.size();


                switch (fi) {
                case 0:
                    namefit = namefit + "_a2";
                    fit_info.Npar = 2;
                    fit_info.function = rhs_a2;
                    break;
                case 1:
                    namefit = namefit + "_a4";
                    fit_info.Npar = 3;
                    fit_info.function = rhs_a2_a4;
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
                fit_info.compute_cov_fit(argv, jackall, lhs_fun);
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
                fit_result amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_fun, fit_info, namefit.c_str());
                fit_info.band_range = { 0,0.0081 };
                std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
                     fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };


                //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
                print_fit_band(argv, jackall, fit_info, fit_info, namefit.c_str(), "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0002, xcont);

                free_fit_result(fit_info, amu_SD_l_common_a4);
                // if (namefit.compare())
            }
        }
    }

    //////////////////////////////////////////////////////////////
    // all ens
    //////////////////////////////////////////////////////////////
}