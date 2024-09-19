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
#include "fve.hpp"
// #include "do_analysis_charm.hpp"

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
    C112,
    B14_64,
    B25_48,
    C20
};


constexpr double fpi_MeV = 130.5;
constexpr double fpi_MeV_err = 0.04;

// constexpr double Mpi_MeV = 135;
// constexpr double Mpi_MeV_err = 0.2;

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
        std::cout << "reading:  " << s << "\n";
        FILE* f = open_file(s.c_str(), "r");

        // read_single_dataj(f, params, &(jackall->en[count]));
        jackall.en[count] = read_single_dataj(f);
        jackall.en[count].resampling = resampling;
        count++;
        fclose(f);
    }
    return jackall;

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


void compute_amul_print_res(char** argv, char* namefit, fit_type fit_info, int Njack, fit_result fit_aMpi2_over_afpi2, double* Mpi2_fpi2_phys, std::vector<std::string> lattices) {
    char file_amul[NAMESIZE];
    mysprintf(file_amul, NAMESIZE, "%s/%s_amul_res.txt", argv[3], namefit);
    FILE* famul = open_file(file_amul, "w+");
    fprintf(famul, "# Lattname amul   damul   \n");

    double** tif = swap_indices(fit_info.Npar, Njack, fit_aMpi2_over_afpi2.P);
    std::vector<double> swapped_x(fit_info.Nvar);
    std::vector<double> amu(Njack);

    int count = 0;
    for (int l = 0; l < lattices.size(); l++) {
        for (size_t j = 0; j < Njack; j++) {
            for (int i = 0; i < fit_info.Nvar; i++) {
                swapped_x[i] = fit_info.x[i][count][j];
            }
            amu[j] = rtbis_func_eq_input(fit_info.function, l /*n*/, fit_info.Nvar, swapped_x.data(), fit_info.Npar, tif[j], 0, Mpi2_fpi2_phys[j], 1e-6, 2, 1e-10, 2);
        }
        printf("amu(%s)=%-20.12g %-20.12g\n",lattices[l].c_str(), amu[Njack - 1], myres->comp_error(amu.data()));
        fprintf(famul, "%s %-20.12g %-20.12g\n",lattices[l].c_str(), amu[Njack - 1], myres->comp_error(amu.data()));

        char file_amul_jack[NAMESIZE];
        mysprintf(file_amul_jack, NAMESIZE, "%s/%s_amul_jack_%s.txt", argv[3], namefit,lattices[l].c_str());
        myres->write_jack_in_file(amu.data(), file_amul_jack);

        count += fit_info.Nxen[l].size();
    }
}



int main(int argc, char** argv) {
    error(argc != 4, 1, "main ",
        "usage:./fit_all_phi4  jack/boot   path_to_jack   output_dir");
    char namefile[NAMESIZE];
    char namefit[NAMESIZE];


    std::vector<std::string> files;
    mysprintf(namefile, NAMESIZE, "%s/%s_cB.72.64_mu.0.000720", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cB.72.96_mu.0.000720", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cC.06.80_mu.0.000600", argv[2], argv[1]);
    files.emplace_back(namefile);
    mysprintf(namefile, NAMESIZE, "%s/%s_cD.54.96_mu.0.000540", argv[2], argv[1]);
    files.emplace_back(namefile);


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

    mysprintf(namefile, NAMESIZE, "%s/%s_cB.14.64_mu.0.001400", argv[2], argv[1]);
    files.emplace_back(namefile);

    mysprintf(namefile, NAMESIZE, "%s/%s_cB.25.48_mu.0.002500", argv[2], argv[1]);
    files.emplace_back(namefile);

    mysprintf(namefile, NAMESIZE, "%s/%s_cC.20.48_mu.0.002000", argv[2], argv[1]);
    files.emplace_back(namefile);

    std::vector<int> myen(files.size());
    for (int e = 0; e < files.size(); e++) {
        myen[e] = e;
        printf("file:  %s\n", files[e].c_str());
    }

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
    double* jack_fpi_phys_MeV = myres->create_fake(fpi_MeV, fpi_MeV_err, 1);
    double* jack_Mpi_phys_MeV = myres->create_fake(Mpi_MeV, Mpi_MeV_err, 2);
    double* Mpi2_fpi2_phys = myres->create_copy(jack_Mpi_phys_MeV);
    myres->div(Mpi2_fpi2_phys, jack_Mpi_phys_MeV, jack_fpi_phys_MeV);
    myres->mult(Mpi2_fpi2_phys, Mpi2_fpi2_phys, Mpi2_fpi2_phys);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.restore_default();
    fit_info.N = 8;
    fit_info.Nvar = 8;
    fit_info.Npar = 7;
    fit_info.Njack = Njack;
    fit_info.Nxen = { {A53, A40, A30}, {B25_48, B14_64, B72_64, B72_96},
                        {C20, C06, C112}, {D54},
                        {E112},  {B72_64},
                        {C06}, {D54} };
    fit_info.init_N_etot_form_Nxen();

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.Nxen.size();n++) {
        for (int e : fit_info.Nxen[n]) {
            for (int j = 0;j < Njack;j++) {
                double my_mu, my_M, my_fpi;
                if (n < 5) {
                    my_mu = jackall.en[e].jack[165][j];
                    my_M = jackall.en[e].jack[1][j];
                    my_fpi = jackall.en[e].jack[163][j];
                }
                else if (n >= 5) {
                    my_mu = jackall.en[e].jack[166][j];
                    my_M = jackall.en[e].jack[123][j];
                    my_fpi = jackall.en[e].jack[164][j];
                }

                fit_info.x[0][count][j] = my_mu; // 
                fit_info.x[1][count][j] = my_M;  // 
                fit_info.x[2][count][j] = my_fpi;  //
                fit_info.x[3][count][j] = jackall.en[e].header.L;

                double xi = my_M / (4 * M_PI * my_fpi);
                xi *= xi;
                double delta_FVE = FVE_GL_Mpi(jackall.en[e].header.L, xi, my_fpi);
                xi *= (1 + delta_FVE) * (1 + delta_FVE) / (1 - 0.25 * delta_FVE) * (1 - 0.25 * delta_FVE);
                fit_info.x[4][count][j] = xi;

                fit_info.x[5][count][j] = jack_Mpi_phys_MeV[j] / hbarc;
                fit_info.x[6][count][j] = jack_fpi_phys_MeV[j] / hbarc;


                fit_info.x[7][count][j] = jack_Mpi_phys_MeV[j] / (4 * M_PI * jack_fpi_phys_MeV[j]);
                fit_info.x[7][count][j] *= fit_info.x[7][count][j];

                // fit_info.x[3][count][j] = jack_Mpi_MeV_exp[j];
                // fit_info.x[4][count][j] = l + 1e-6;
                // fit_info.x[5][count][j] = a + 1e-6;
                // fit_info.x[6][count][j] = 0 + 1e-6;
                // fit_info.x[7][count][j] = w + 1.0;
            }
            count++;
        }
    }

    // fit_info.corr_id = { 1, 123, 163, 164 }; // Mpi(mu1), Mpi(mu2), fpi(mu1), fpi(mu2)
    fit_info.corr_id = { 163, 163, 163, 163, 163, 164, 164, 164 }; //  fpi(mu1), fpi(mu2)
    fit_info.function = rhs_afpi;
    fit_info.linear_fit = false;
    // fit_info.covariancey = false;
    // fit_info.acc= 1e-10;
    // fit_info.h=1e-7;
    // fit_info.maxiter=500;
    // fit_info.NM = true;
    // // fit_info.chi2_gap_jackboot=0.1;
    // // fit_info.guess_per_jack=5;
    // // fit_info.repeat_start=100;
    // fit_info.verbosity = 1;
    // fit_info.compute_cov_fit(argv, jackall, lhs_afpi);
    // int ide = 0, ide1 = 0;
    // for (int n = 0;n < fit_info.Nxen.size();n++) {
    //     for (int e : fit_info.Nxen[n]) {
    //         ide1 = 0;
    //         for (int n1 = 0;n1 < fit_info.Nxen.size();n1++) {
    //             for (int e1 : fit_info.Nxen[n1]) {
    //                 if (e != e1)   fit_info.cov[ide][ide1] = 0;
    //                 ide1++;
    //             }
    //         }
    //         ide++;
    //     }
    // }
    // fit_info.compute_cov1_fit();
    fit_info.guess = { 0.0908026, 0.07951 , 0.06816, 0.05688, 0.04891 ,-8.75 ,-500 };
    mysprintf(namefit, NAMESIZE, "afpi_cov");
    fit_result fit_afpi = fit_all_data(argv, jackall, lhs_afpi, fit_info, namefit);
    fit_info.band_range = { 0.005,0.035 };
    // std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
    //      fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };


    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    std::vector<double> xcont = {};
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "xi", fit_afpi, fit_afpi, 4, 0, 0.0005, xcont);

    myres->write_jack_in_file(fit_afpi.P[0], "../../g-2_new_stat/out/a_fm_A.txt");
    myres->write_jack_in_file(fit_afpi.P[1], "../../g-2_new_stat/out/a_fm_B.txt");
    myres->write_jack_in_file(fit_afpi.P[2], "../../g-2_new_stat/out/a_fm_C.txt");
    myres->write_jack_in_file(fit_afpi.P[3], "../../g-2_new_stat/out/a_fm_D.txt");
    myres->write_jack_in_file(fit_afpi.P[4], "../../g-2_new_stat/out/a_fm_E.txt");

    // fit_afpi.clear();


    //////////////////////////////////////////////////////////////
    // fit only physical point
    //////////////////////////////////////////////////////////////

    fit_info.restore_default();
    fit_info.N = 7;
    fit_info.Nvar = 8;
    fit_info.Npar = 6;
    fit_info.Njack = Njack;
    fit_info.Nxen = { { B72_64, B72_96},
                        { C06}, {D54},
                        {E112},  {B72_64},
                        {C06}, {D54} };
    fit_info.init_N_etot_form_Nxen();

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.Nxen.size();n++) {
        for (int e : fit_info.Nxen[n]) {
            for (int j = 0;j < Njack;j++) {
                double my_mu, my_M, my_fpi;
                if (n < 4) {
                    my_mu = jackall.en[e].jack[165][j];
                    my_M = jackall.en[e].jack[1][j];
                    my_fpi = jackall.en[e].jack[163][j];
                }
                else if (n >= 4) {
                    my_mu = jackall.en[e].jack[166][j];
                    my_M = jackall.en[e].jack[123][j];
                    my_fpi = jackall.en[e].jack[164][j];
                }

                fit_info.x[0][count][j] = my_mu; // 
                fit_info.x[1][count][j] = my_M;  // 
                fit_info.x[2][count][j] = my_fpi;  //
                fit_info.x[3][count][j] = jackall.en[e].header.L;

                double xi = my_M / (4 * M_PI * my_fpi);
                xi *= xi;
                double delta_FVE = FVE_GL_Mpi(jackall.en[e].header.L, xi, my_fpi);
                xi *= (1 + delta_FVE) * (1 + delta_FVE) / (1 - 0.25 * delta_FVE) * (1 - 0.25 * delta_FVE);
                fit_info.x[4][count][j] = xi;

                fit_info.x[5][count][j] = jack_Mpi_phys_MeV[j] / hbarc;
                fit_info.x[6][count][j] = jack_fpi_phys_MeV[j] / hbarc;


                fit_info.x[7][count][j] = jack_Mpi_phys_MeV[j] / (4 * M_PI * jack_fpi_phys_MeV[j]);
                fit_info.x[7][count][j] *= fit_info.x[7][count][j];

            }
            count++;
        }
    }

    fit_info.corr_id = { 163, 163, 163, 163, 164, 164, 164 }; //  fpi(mu1), fpi(mu2)
    // fit_info.corr_id = {  1, 1, 1, 1, 123, 123, 123 }; //  Mpi(mu1), .... , Mpi(mu2), ...
    fit_info.function = rhs_afpi_phys_point;
    fit_info.linear_fit = false;
    fit_info.covariancey = false;
    // fit_info.acc= 1e-6;
    // fit_info.chi2_gap_jackboot=0.1;
    // fit_info.guess_per_jack=5;
    // fit_info.repeat_start=5;
    // fit_info.verbosity = 0;
    // fit_info.compute_cov_fit(argv, jackall, lhs_afpi);
    // ie = 0, ie1 = 0;
    // for (int n = 0;n < fit_info.N;n++) {
    //     for (int e = 0;e < fit_info.myen.size();e++) {
    //         ie1 = 0;
    //         for (int n1 = 0;n1 < fit_info.N;n1++) {
    //             for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
    //                 if (e != e1)   fit_info.cov[ie][ie1] = 0;
    //                 ie1++;
    //             }
    //         }
    //         ie++;
    //     }
    // }
    // fit_info.compute_cov1_fit();
    fit_info.guess = { 0.07951 , 0.06816, 0.05688, 0.04891 ,1 ,1 };
    mysprintf(namefit, NAMESIZE, "afpi_phys_point_cov");
    fit_result fit_afpi_phys_point = fit_all_data(argv, jackall, lhs_afpi, fit_info, namefit);
    fit_info.band_range = { 0.005,0.01 };

    xcont = {};
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "xi", fit_afpi_phys_point, fit_afpi_phys_point, 4, 0, 0.0005, xcont);

    // myres->write_jack_in_file(fit_afpi.P[0], "../../g-2_new_stat/out/a_fm_A.txt");
    // myres->write_jack_in_file(fit_afpi.P[1], "../../g-2_new_stat/out/a_fm_B.txt");
    // myres->write_jack_in_file(fit_afpi.P[2], "../../g-2_new_stat/out/a_fm_C.txt");
    // myres->write_jack_in_file(fit_afpi.P[3], "../../g-2_new_stat/out/a_fm_D.txt");
    // myres->write_jack_in_file(fit_afpi.P[4], "../../g-2_new_stat/out/a_fm_E.txt");

    fit_afpi_phys_point.clear();

    ///////////////////////////////////////////////////////////////////////////////////////////////
    fit_info.Npar = 8;
    fit_info.function = rhs_afpi_FVEres;

    fit_info.guess = { 0.0908026, 0.07951 , 0.06816, 0.05688, 0.04891 ,1 ,1 };
    mysprintf(namefit, NAMESIZE, "afpi_resFVE_cov");
    fit_result fit_afpi_1 = fit_all_data(argv, jackall, lhs_afpi, fit_info, namefit);
    fit_info.band_range = { 0.005,0.035 };
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "xi", fit_afpi_1, fit_afpi_1, 4, 0, 0.0005, xcont);

    mysprintf(namefit, NAMESIZE, "afpi_resFVE_corr_cov");
    fit_info.function = rhs_afpi;
    print_data_fit_corrected(argv, jackall, lhs_afpi_remove_FVE, fit_info, namefit, fit_afpi_1);
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "xi", fit_afpi_1, fit_afpi_1, 4, 0, 0.0005, xcont);

    fit_afpi_1.clear();

    //////////////////////////////////////////////////////////////
    // exclude C20
    //////////////////////////////////////////////////////////////
    fit_info.restore_default();
    fit_info.N = 8;
    fit_info.Nvar = 8;
    fit_info.Npar = 7;
    fit_info.Njack = Njack;
    fit_info.Nxen = { {A53, A40, A30}, {B25_48, B14_64, B72_64, B72_96},
                        { C06, C112}, {D54},
                        {E112},  {B72_64},
                        {C06}, {D54} };
    fit_info.init_N_etot_form_Nxen();

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.Nxen.size();n++) {
        for (int e : fit_info.Nxen[n]) {
            for (int j = 0;j < Njack;j++) {
                double my_mu, my_M, my_fpi;
                if (n < 5) {
                    my_mu = jackall.en[e].jack[165][j];
                    my_M = jackall.en[e].jack[1][j];
                    my_fpi = jackall.en[e].jack[163][j];
                }
                else if (n >= 5) {
                    my_mu = jackall.en[e].jack[166][j];
                    my_M = jackall.en[e].jack[123][j];
                    my_fpi = jackall.en[e].jack[164][j];
                }

                fit_info.x[0][count][j] = my_mu; // 
                fit_info.x[1][count][j] = my_M;  // 
                fit_info.x[2][count][j] = my_fpi;  //
                fit_info.x[3][count][j] = jackall.en[e].header.L;

                double xi = my_M / (4 * M_PI * my_fpi);
                xi *= xi;
                double delta_FVE = FVE_GL_Mpi(jackall.en[e].header.L, xi, my_fpi);
                xi *= (1 + delta_FVE) * (1 + delta_FVE) / (1 - 0.25 * delta_FVE) * (1 - 0.25 * delta_FVE);
                fit_info.x[4][count][j] = xi;

                fit_info.x[5][count][j] = jack_Mpi_phys_MeV[j] / hbarc;
                fit_info.x[6][count][j] = jack_fpi_phys_MeV[j] / hbarc;


                fit_info.x[7][count][j] = jack_Mpi_phys_MeV[j] / (4 * M_PI * jack_fpi_phys_MeV[j]);
                fit_info.x[7][count][j] *= fit_info.x[7][count][j];

            }
            count++;
        }
    }

    fit_info.corr_id = { 163, 163, 163, 163, 163, 164, 164, 164 }; //  fpi(mu1), fpi(mu2)
    fit_info.function = rhs_afpi;
    fit_info.linear_fit = false;
    fit_info.covariancey = false;
    fit_info.guess = { 0.0908026, 0.07951 , 0.06816, 0.05688, 0.04891 ,1 ,1 };
    mysprintf(namefit, NAMESIZE, "afpi_cov_noC20");
    fit_afpi_1 = fit_all_data(argv, jackall, lhs_afpi, fit_info, namefit);
    fit_info.band_range = { 0.005,0.035 };
    xcont = {};
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "xi", fit_afpi_1, fit_afpi_1, 4, 0, 0.0005, xcont);

    fit_afpi_1.clear();

    //////////////////////////////////////////////////////////////
    // Mpi
    //////////////////////////////////////////////////////////////
    {

        fit_info.restore_default();
        fit_info.N = 7;
        fit_info.Nvar = 9;
        fit_info.Npar = 5;
        fit_info.Njack = Njack;
        fit_info.Nxen = { {B25_48, B14_64, B72_64, B72_96},
                            {C20, C06, C112}, {D54},
                            {E112},  {B72_64},
                            {C06}, {D54} };
        fit_info.init_N_etot_form_Nxen();

        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
        count = 0;
        for (int n = 0;n < fit_info.Nxen.size();n++) {
            for (int e : fit_info.Nxen[n]) {
                for (int j = 0;j < Njack;j++) {
                    double my_mu, my_M, my_fpi;
                    if (n < 4) {
                        my_mu = jackall.en[e].jack[165][j];
                        my_M = jackall.en[e].jack[1][j];
                        my_fpi = jackall.en[e].jack[163][j];
                    }
                    else if (n >= 4) {
                        my_mu = jackall.en[e].jack[166][j];
                        my_M = jackall.en[e].jack[123][j];
                        my_fpi = jackall.en[e].jack[164][j];
                    }
                    fit_info.x[0][count][j] = my_mu; // 
                    fit_info.x[1][count][j] = my_M;  // 
                    fit_info.x[2][count][j] = my_fpi;  //
                    fit_info.x[3][count][j] = jackall.en[e].header.L;

                    double xi = my_M / (4 * M_PI * my_fpi);
                    xi *= xi;
                    double delta_FVE = FVE_GL_Mpi(jackall.en[e].header.L, xi, my_fpi);
                    xi *= (1 + delta_FVE) * (1 + delta_FVE) / (1 - 0.25 * delta_FVE) * (1 - 0.25 * delta_FVE);
                    fit_info.x[4][count][j] = xi;

                    fit_info.x[5][count][j] = jack_Mpi_phys_MeV[j] / hbarc;
                    fit_info.x[6][count][j] = jack_fpi_phys_MeV[j] / hbarc;


                    fit_info.x[7][count][j] = jack_Mpi_phys_MeV[j] / (4 * M_PI * jack_fpi_phys_MeV[j]);
                    fit_info.x[7][count][j] *= fit_info.x[7][count][j];

                    if (e == A53 || e == A40 || e == A30)
                        fit_info.x[8][count][j] = fit_afpi.P[0][j];
                    else if (e == B25_48 || e == B14_64 || e == B72_64 || e == B72_96)
                        fit_info.x[8][count][j] = fit_afpi.P[1][j];
                    else if (e == C20 || e == C06 || e == C112)
                        fit_info.x[8][count][j] = fit_afpi.P[2][j];
                    else if (e == D54)
                        fit_info.x[8][count][j] = fit_afpi.P[3][j];
                    else if (e == E112)
                        fit_info.x[8][count][j] = fit_afpi.P[4][j];


                }
                count++;
            }
        }



        // fit_info.corr_id = { 1, 123, 163, 164 }; // Mpi(mu1), Mpi(mu2), fpi(mu1), fpi(mu2)
        fit_info.corr_id = { 1, 1, 1, 1, 123, 123, 123 }; //  Mpi(mu1), Mpi(mu2)
        fit_info.function = rhs_aMpi2_over_afpi2;
        fit_info.linear_fit = false;
        fit_info.covariancey = false;
        // fit_info.acc= 1e-8;
        // fit_info.h=1e-7;
        // fit_info.maxiter=500;
        // fit_info.NM = true;
        // // fit_info.chi2_gap_jackboot=0.1;
        // // fit_info.guess_per_jack=5;
        // fit_info.repeat_start=1000;
        // fit_info.verbosity = 4;
        // fit_info.compute_cov_fit(argv, jackall, lhs_Mpi2_over_afpi2);
        // int ide = 0, ide1 = 0;
        // for (int n = 0;n < fit_info.Nxen.size();n++) {
        //     for (int e : fit_info.Nxen[n]) {
        //         ide1 = 0;
        //         for (int n1 = 0;n1 < fit_info.Nxen.size();n1++) {
        //             for (int e1 : fit_info.Nxen[n1]) {
        //                 if (e != e1)   fit_info.cov[ide][ide1] = 0;
        //                 ide1++;
        //             }
        //         }
        //         ide++;
        //     }
        // }
        // fit_info.compute_cov1_fit();
        // fit_info.guess = { 0.0072, 0.0060 , 0.0050, 0.0012, -8.75 , -1 };
        // fit_info.guess = { 866, 990 , 1169, 1346, -5.46, 0  };
        fit_info.guess = { 1, 1 , 1, 1, -8.75 };
        // fit_info.guess = { 1, 1 , 1, 1  };
        // fit_info.guess = { 0.05 };
        mysprintf(namefit, NAMESIZE, "aMpi2_over_afpi2_cov");
        fit_result fit_aMpi2_over_afpi2 = fit_all_data(argv, jackall, lhs_Mpi2_over_afpi2, fit_info, namefit);
        fit_info.band_range = { 0.00001,0.003 };

        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "amu", fit_aMpi2_over_afpi2, fit_aMpi2_over_afpi2, 0 /*amu */, 0, 0.0001, xcont);

        compute_amul_print_res( argv,  namefit,  fit_info, Njack, fit_aMpi2_over_afpi2,  Mpi2_fpi2_phys, {"B","C","D","E"});

    }
    //////////////////////////////////////////////////////////////
     // Mpi with A
     //////////////////////////////////////////////////////////////

    {
        fit_info.restore_default();
        fit_info.N = 8;
        fit_info.Nvar = 9;
        fit_info.Npar = 6;
        fit_info.Njack = Njack;
        fit_info.Nxen = { {A53, A40, A30}, {B25_48, B14_64, B72_64, B72_96},
                            {C20, C06, C112}, {D54},
                            {E112},  {B72_64},
                            {C06}, {D54} };
        fit_info.init_N_etot_form_Nxen();

        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
        count = 0;
        for (int n = 0;n < fit_info.Nxen.size();n++) {
            for (int e : fit_info.Nxen[n]) {
                for (int j = 0;j < Njack;j++) {
                    double my_mu, my_M, my_fpi;
                    if (n < 5) {
                        my_mu = jackall.en[e].jack[165][j];
                        my_M = jackall.en[e].jack[1][j];
                        my_fpi = jackall.en[e].jack[163][j];
                    }
                    else if (n >= 5) {
                        my_mu = jackall.en[e].jack[166][j];
                        my_M = jackall.en[e].jack[123][j];
                        my_fpi = jackall.en[e].jack[164][j];
                    }
                    fit_info.x[0][count][j] = my_mu; // 
                    fit_info.x[1][count][j] = my_M;  // 
                    fit_info.x[2][count][j] = my_fpi;  //
                    fit_info.x[3][count][j] = jackall.en[e].header.L;

                    double xi = my_M / (4 * M_PI * my_fpi);
                    xi *= xi;
                    double delta_FVE = FVE_GL_Mpi(jackall.en[e].header.L, xi, my_fpi);
                    xi *= (1 + delta_FVE) * (1 + delta_FVE) / (1 - 0.25 * delta_FVE) * (1 - 0.25 * delta_FVE);
                    fit_info.x[4][count][j] = xi;

                    fit_info.x[5][count][j] = jack_Mpi_phys_MeV[j] / hbarc;
                    fit_info.x[6][count][j] = jack_fpi_phys_MeV[j] / hbarc;


                    fit_info.x[7][count][j] = jack_Mpi_phys_MeV[j] / (4 * M_PI * jack_fpi_phys_MeV[j]);
                    fit_info.x[7][count][j] *= fit_info.x[7][count][j];

                    if (e == A53 || e == A40 || e == A30)
                        fit_info.x[8][count][j] = fit_afpi.P[0][j];
                    else if (e == B25_48 || e == B14_64 || e == B72_64 || e == B72_96)
                        fit_info.x[8][count][j] = fit_afpi.P[1][j];
                    else if (e == C20 || e == C06 || e == C112)
                        fit_info.x[8][count][j] = fit_afpi.P[2][j];
                    else if (e == D54)
                        fit_info.x[8][count][j] = fit_afpi.P[3][j];
                    else if (e == E112)
                        fit_info.x[8][count][j] = fit_afpi.P[4][j];

                }
                count++;
            }
        }



        // fit_info.corr_id = { 1, 123, 163, 164 }; // Mpi(mu1), Mpi(mu2), fpi(mu1), fpi(mu2)
        fit_info.corr_id = { 1, 1, 1, 1, 1, 123, 123, 123 }; //  Mpi(mu1), Mpi(mu2)
        fit_info.function = rhs_aMpi2_over_afpi2_with_A;
        fit_info.linear_fit = false;
        fit_info.covariancey = false;
        // fit_info.acc= 1e-8;
        // fit_info.h=1e-7;
        // fit_info.maxiter=500;
        // fit_info.NM = true;
        // // fit_info.chi2_gap_jackboot=0.1;
        // // fit_info.guess_per_jack=5;
        // fit_info.repeat_start=1000;
        // fit_info.verbosity = 4;
        // fit_info.compute_cov_fit(argv, jackall, lhs_Mpi2_over_afpi2);
        // int ide = 0, ide1 = 0;
        // for (int n = 0;n < fit_info.Nxen.size();n++) {
        //     for (int e : fit_info.Nxen[n]) {
        //         ide1 = 0;
        //         for (int n1 = 0;n1 < fit_info.Nxen.size();n1++) {
        //             for (int e1 : fit_info.Nxen[n1]) {
        //                 if (e != e1)   fit_info.cov[ide][ide1] = 0;
        //                 ide1++;
        //             }
        //         }
        //         ide++;
        //     }
        // }
        // fit_info.compute_cov1_fit();
        // fit_info.guess = { 0.0072, 0.0060 , 0.0050, 0.0012, -8.75 , -1 };
        // fit_info.guess = { 866, 990 , 1169, 1346, -5.46, 0  };
        fit_info.guess = { 1, 1, 1 , 1, 1, -8.75 };
        // fit_info.guess = { 1, 1 , 1, 1  };
        // fit_info.guess = { 0.05 };
        mysprintf(namefit, NAMESIZE, "aMpi2_over_afpi2_A_cov");
        fit_result fit_aMpi2_over_afpi2 = fit_all_data(argv, jackall, lhs_Mpi2_over_afpi2, fit_info, namefit);
        fit_info.band_range = { 0.00001,0.0055 };

        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "amu", fit_aMpi2_over_afpi2, fit_aMpi2_over_afpi2, 0 /*amu */, 0, 0.0001, xcont);

        compute_amul_print_res( argv,  namefit,  fit_info, Njack, fit_aMpi2_over_afpi2,  Mpi2_fpi2_phys, {"A","B","C","D","E"});


    }
    //////////////////////////////////////////////////////////////
    // Mpi with A a2
    //////////////////////////////////////////////////////////////

    {
        fit_info.restore_default();
        fit_info.N = 8;
        fit_info.Nvar = 9;
        fit_info.Npar = 7;
        fit_info.Njack = Njack;
        fit_info.Nxen = { {A53, A40, A30}, {B25_48, B14_64, B72_64, B72_96},
                            {C20, C06, C112}, {D54},
                            {E112},  {B72_64},
                            {C06}, {D54} };
        fit_info.init_N_etot_form_Nxen();

        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
        count = 0;
        for (int n = 0;n < fit_info.Nxen.size();n++) {
            for (int e : fit_info.Nxen[n]) {
                for (int j = 0;j < Njack;j++) {
                    double my_mu, my_M, my_fpi;
                    if (n < 5) {
                        my_mu = jackall.en[e].jack[165][j];
                        my_M = jackall.en[e].jack[1][j];
                        my_fpi = jackall.en[e].jack[163][j];
                    }
                    else if (n >= 5) {
                        my_mu = jackall.en[e].jack[166][j];
                        my_M = jackall.en[e].jack[123][j];
                        my_fpi = jackall.en[e].jack[164][j];
                    }
                    fit_info.x[0][count][j] = my_mu; // 
                    fit_info.x[1][count][j] = my_M;  // 
                    fit_info.x[2][count][j] = my_fpi;  //
                    fit_info.x[3][count][j] = jackall.en[e].header.L;

                    double xi = my_M / (4 * M_PI * my_fpi);
                    xi *= xi;
                    double delta_FVE = FVE_GL_Mpi(jackall.en[e].header.L, xi, my_fpi);
                    xi *= (1 + delta_FVE) * (1 + delta_FVE) / (1 - 0.25 * delta_FVE) * (1 - 0.25 * delta_FVE);
                    fit_info.x[4][count][j] = xi;

                    fit_info.x[5][count][j] = jack_Mpi_phys_MeV[j] / hbarc;
                    fit_info.x[6][count][j] = jack_fpi_phys_MeV[j] / hbarc;


                    fit_info.x[7][count][j] = jack_Mpi_phys_MeV[j] / (4 * M_PI * jack_fpi_phys_MeV[j]);
                    fit_info.x[7][count][j] *= fit_info.x[7][count][j];

                    if (e == A53 || e == A40 || e == A30)
                        fit_info.x[8][count][j] = fit_afpi.P[0][j];
                    else if (e == B25_48 || e == B14_64 || e == B72_64 || e == B72_96)
                        fit_info.x[8][count][j] = fit_afpi.P[1][j];
                    else if (e == C20 || e == C06 || e == C112)
                        fit_info.x[8][count][j] = fit_afpi.P[2][j];
                    else if (e == D54)
                        fit_info.x[8][count][j] = fit_afpi.P[3][j];
                    else if (e == E112)
                        fit_info.x[8][count][j] = fit_afpi.P[4][j];

                }
                count++;
            }
        }



        // fit_info.corr_id = { 1, 123, 163, 164 }; // Mpi(mu1), Mpi(mu2), fpi(mu1), fpi(mu2)
        fit_info.corr_id = { 1, 1, 1, 1, 1, 123, 123, 123 }; //  Mpi(mu1), Mpi(mu2)
        fit_info.function = rhs_aMpi2_over_afpi2_with_A;
        fit_info.linear_fit = false;
        fit_info.covariancey = false;
        // fit_info.acc= 1e-8;
        // fit_info.h=1e-7;
        // fit_info.maxiter=500;
        // fit_info.NM = true;
        // // fit_info.chi2_gap_jackboot=0.1;
        // // fit_info.guess_per_jack=5;
        // fit_info.repeat_start=1000;
        fit_info.verbosity = 1;
        // fit_info.compute_cov_fit(argv, jackall, lhs_Mpi2_over_afpi2);
        // int ide = 0, ide1 = 0;
        // for (int n = 0;n < fit_info.Nxen.size();n++) {
        //     for (int e : fit_info.Nxen[n]) {
        //         ide1 = 0;
        //         for (int n1 = 0;n1 < fit_info.Nxen.size();n1++) {
        //             for (int e1 : fit_info.Nxen[n1]) {
        //                 if (e != e1)   fit_info.cov[ide][ide1] = 0;
        //                 ide1++;
        //             }
        //         }
        //         ide++;
        //     }
        // }
        // fit_info.compute_cov1_fit();
        // fit_info.guess = { 0.0072, 0.0060 , 0.0050, 0.0012, -8.75 , -1 };
        // fit_info.guess = { 866, 990 , 1169, 1346, -5.46, 0  };
        fit_info.guess = { 812, 889, 1015 , 1200, 1381, 1.4,0 };
        // fit_info.guess = { 1, 1 , 1, 1  };
        // fit_info.guess = { 0.05 };
        mysprintf(namefit, NAMESIZE, "aMpi2_over_afpi2_a2_A_cov");
        fit_result fit_aMpi2_over_afpi2 = fit_all_data(argv, jackall, lhs_Mpi2_over_afpi2, fit_info, namefit);
        fit_info.band_range = { 0.00001,0.0055 };

        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "amu", fit_aMpi2_over_afpi2, fit_aMpi2_over_afpi2, 0 /*amu */, 0, 0.0001, xcont);
        
        compute_amul_print_res( argv,  namefit,  fit_info, Njack, fit_aMpi2_over_afpi2,  Mpi2_fpi2_phys, {"A","B","C","D","E"});

    }
    ///////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////
    // mK
    //////////////////////////////////////////////////////////////
}