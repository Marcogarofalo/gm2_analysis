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





void   do_analysis(char** argv, std::vector<int> ids, std::vector<std::string> M, std::string basename, data_all jackall) {

    fit_type fit_info;
    int Njack = jackall.en[0].Njack;
    int count = 0;
    double* jack_Mpi_MeV_exp = fake_sampling(argv[1], Mpi_MeV, Mpi_MeV_err, Njack, 1003);
    char namefit[NAMESIZE];
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n/////////////////////////////////     amu_W_l Lref   //////////////////\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////
    data_all  syst_amu_W_lphys_Lref;
    syst_amu_W_lphys_Lref.resampling = argv[1];
    std::string save_basename(basename);

    for (int iM = 0; iM < M.size(); iM++) {
        basename = save_basename + "_" + M[iM];
        ////// separate fits
        std::vector<std::string> integrations = { "reinman" };
        for (auto integration : integrations) {
            int id0, id1;
            if (integration == "reinman") { id0 = ids[0 + iM * 2]; id1 = ids[1 + iM * 2]; }

            for (int l = 0;l < 4;l++) {
                for (int a : { 1}) {
                    for (int w = 0;w < 2;w++) {
                        for (int OSTM = 0; OSTM < 2;OSTM++) {
                            fit_info.restore_default();
                            fit_info.Npar = 1;
                            if (integration == "reinman" && OSTM == 0) { id0 = ids[0 + iM * 2]; }
                            if (integration == "reinman" && OSTM == 1) { id0 = ids[1 + iM * 2]; }

                            if (a > 0) fit_info.Npar += 1;
                            if (a == 0 && l >= 1) continue;
                            if (a == 0 && w >= 1) continue;


                            fit_info.N = 1;
                            fit_info.Nvar = 8;
                            fit_info.Njack = Njack;
                            fit_info.myen = { B72_96, C06, D54 ,E112 };
                            if (fit_info.Npar >= fit_info.myen.size() * fit_info.N) { continue; }

                            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
                            count = 0;
                            for (int n = 0;n < fit_info.N;n++) {
                                for (int e : fit_info.myen) {
                                    for (int j = 0;j < Njack;j++) {
                                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
                                        fit_info.x[1][count][j] = jackall.en[e].jack[58][j];  // Delta_FV_GS
                                        fit_info.x[2][count][j] = jackall.en[e].jack[1][j];  //Mpi
                                        fit_info.x[3][count][j] = jack_Mpi_MeV_exp[j];
                                        fit_info.x[4][count][j] = l + 1e-6;
                                        fit_info.x[5][count][j] = a + 1e-6;
                                        fit_info.x[6][count][j] = 0 + 1e-6;
                                        fit_info.x[7][count][j] = w + 1e-6;
                                    }
                                    count++;
                                }
                            }
                            fit_info.corr_id = { id0 };
                            fit_info.function = rhs_amu_separate;
                            fit_info.linear_fit = true;
                            fit_info.covariancey = true;
                            // fit_info.acc= 1e-6;
                            // fit_info.chi2_gap_jackboot=0.1;
                            // fit_info.guess_per_jack=5;
                            // fit_info.repeat_start=5;
                            fit_info.verbosity = 0;
                            fit_info.compute_cov_fit(argv, jackall, lhs_amu_separate);
                            int ie = 0, ie1 = 0;
                            for (int n = 0;n < fit_info.N;n++) {
                                for (int e = 0;e < fit_info.myen.size();e++) {
                                    ie1 = 0;
                                    for (int n1 = 0;n1 < fit_info.N;n1++) {
                                        for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                                            if (e != e1)   fit_info.cov[ie][ie1] = 0;
                                            ie1++;
                                        }
                                    }
                                    ie++;
                                }
                            }
                            fit_info.compute_cov1_fit();

                            std::string logname;
                            if (l == 0) { logname = ""; }
                            if (l > 0) { logname = "log" + std::to_string(l); }


                            if (l == 0 && w > 0) continue;
                            std::string wname;
                            if (w == 0) { wname = "w1"; }
                            if (w == 1) { wname = "w3"; }



                            std::string aname;
                            if (a == 0) { aname = ""; }
                            if (a == 1) { aname = "a2"; }


                            std::string regname;
                            if (OSTM == 0) { regname = "OS"; }
                            if (OSTM == 1) { regname = "TM"; }

                            mysprintf(namefit, NAMESIZE, "amu_%s_separate_%s_%s_%s_%s_cov", basename.c_str(), regname.c_str(), logname.c_str(), wname.c_str(), aname.c_str());
                            fit_result amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu_separate, fit_info, namefit);
                            fit_info.band_range = { 0,0.0081 };
                            std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
                                 fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

                            // TODO: in order to print the band you need to subtract the
                            //    FVE ok
                            //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
                            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
                            syst_amu_W_lphys_Lref.add_fit(amu_SD_l_common_a4);

                            // if(iM==3 && a ==4 ) exit(1);
                            free_fit_result(fit_info, amu_SD_l_common_a4);
                            fit_info.restore_default();
                        }


                    }
                }
            }
        }
        fit_info.restore_default();

        ////// cuts fits
        integrations = { "reinman" };
        for (auto integration : integrations) {
            int id0, id1;
            if (integration == "reinman") { id0 = ids[0 + iM * 2]; id1 = ids[1 + iM * 2]; }

            for (int l : {0, 1, 2, 3, 4, 5, 8, 10, 12, 15}) {
                // for (int l=0; l<16 ;l++) {
                for (int a : { 1}) {
                    for (int w = 0;w < 2;w++) {
                        for (int icut = 0; icut < 4;icut++) {
                            fit_info.restore_default();
                            fit_info.Npar = 1;
                            if (integration == "reinman") { id0 = ids[0 + iM * 2]; id1 = ids[1 + iM * 2]; }


                            if (a > 0) fit_info.Npar += 2;
                            if (a == 0 && l >= 1) { fit_info.restore_default(); continue; }
                            if (a == 0 && w >= 1) { fit_info.restore_default();continue; }


                            fit_info.N = 2;
                            fit_info.Nvar = 8;
                            fit_info.Njack = Njack;
                            // fit_info.myen = { B72_64, C06, D54 };
                            if (icut == 0)
                                fit_info.Nxen = { {B72_96, C06 ,D54, E112},
                                                   {B72_96, C06, D54 ,E112} };
                            if (icut == 1)
                                fit_info.Nxen = { { C06 ,D54, E112},
                                                   {B72_96, C06, D54 ,E112} };
                            if (icut == 2)
                                fit_info.Nxen = { {B72_96, C06 ,D54 ,E112},
                                                   { C06, D54, E112} };
                            if (icut == 3)
                                fit_info.Nxen = { { C06 ,D54, E112},
                                                   { C06, D54, E112} };


                            fit_info.init_N_etot_form_Nxen();
                            if (fit_info.Npar >= fit_info.entot) { fit_info.restore_default(); continue; }

                            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
                            count = 0;
                            for (int n = 0;n < fit_info.Nxen.size();n++) {
                                for (int e : fit_info.Nxen[n]) {
                                    for (int j = 0;j < Njack;j++) {
                                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
                                        fit_info.x[1][count][j] = jackall.en[e].jack[58][j];  // Delta_FV_GS
                                        fit_info.x[2][count][j] = jackall.en[e].jack[1][j];  //Mpi
                                        fit_info.x[3][count][j] = jack_Mpi_MeV_exp[j];
                                        fit_info.x[4][count][j] = l + 1e-6;
                                        fit_info.x[5][count][j] = a + 1e-6;
                                        fit_info.x[6][count][j] = 0 + 1e-6;
                                        fit_info.x[7][count][j] = w + 1.0;
                                    }
                                    count++;
                                }
                            }
                            fit_info.corr_id = { id0, id1 };
                            fit_info.function = rhs_amu_cut;
                            fit_info.linear_fit = true;
                            fit_info.covariancey = true;
                            // fit_info.acc= 1e-6;
                            // fit_info.chi2_gap_jackboot=0.1;
                            // fit_info.guess_per_jack=5;
                            // fit_info.repeat_start=5;
                            fit_info.verbosity = 0;
                            fit_info.compute_cov_fit(argv, jackall, lhs_amu);
                            int ie = 0, ie1 = 0;
                            for (int n = 0;n < fit_info.N;n++) {
                                for (int e = 0;e < fit_info.Nxen[n].size();e++) {
                                    ie1 = 0;
                                    for (int n1 = 0;n1 < fit_info.N;n1++) {
                                        for (int e1 = 0;e1 < fit_info.Nxen[n1].size();e1++) {
                                            if (e != e1)   fit_info.cov[ie][ie1] = 0;
                                            ie1++;
                                        }
                                    }
                                    ie++;
                                }
                            }
                            fit_info.compute_cov1_fit();

                            std::string logname;
                            if (l == 0) { logname = ""; }
                            if (l > 0) { logname = "log" + std::to_string(l); }

                            if (l == 0 && w > 0) { fit_info.restore_default();continue; }
                            std::string wname;
                            if (w == 0) { wname = "w1"; }
                            if (w == 1) { wname = "w3"; }



                            std::string aname;
                            if (a == 0) { aname = ""; }
                            if (a == 1) { aname = "a2"; }


                            std::string regname;
                            if (icut == 0) { regname = ""; }
                            if (icut == 1) { regname = "OS2"; }
                            if (icut == 2) { regname = "TM2"; }
                            if (icut == 3) { regname = "OS2_TM2"; }

                            mysprintf(namefit, NAMESIZE, "amu_%s_%s_%s_%s_%s_cov", basename.c_str(), regname.c_str(), logname.c_str(), wname.c_str(), aname.c_str());
                            fit_result amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
                            fit_info.band_range = { 0,0.0081 };
                            std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
                                 fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

                            // TODO: in order to print the band you need to subtract the
                            //    FVE ok
                            //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
                            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
                            syst_amu_W_lphys_Lref.add_fit(amu_SD_l_common_a4);

                            // if(iM==3 && a ==4 ) exit(1);
                            free_fit_result(fit_info, amu_SD_l_common_a4);
                            fit_info.restore_default();
                        }


                    }
                }
            }
        }

        ////////////////////////////// Padè fit
        // std::map<int, std::string> stringOSTM = { {0,"OS"}, {1,"TM"} };

        // for (int OSTM = 0;OSTM < 2;OSTM++) {
        // for (int l = 0;l < 4;l++) {
        //     for (int w = 0;w < 2;w++) {
        //         fit_info.restore_default();
        //         fit_info.Npar = 4;
        //         int id0 = ids[0 + iM * 2], id1 = ids[1 + iM * 2];
        //         fit_info.N = 2;
        //         fit_info.Nvar = 8;
        //         fit_info.Njack = Njack;
        //         fit_info.myen = { B72_64, C06, D54 };
        //         fit_info.init_Nxen_from_N_myen();
        //         if (fit_info.Npar >= fit_info.myen.size() * fit_info.N) { exit(1); }

        //         fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
        //         count = 0;
        //         for (int n = 0;n < fit_info.Nxen.size();n++) {
        //             for (int e : fit_info.Nxen[n]) {
        //                 for (int j = 0;j < Njack;j++) {
        //                     fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
        //                     fit_info.x[1][count][j] = jackall.en[e].jack[58][j];  // Delta_FV_GS
        //                     fit_info.x[2][count][j] = jackall.en[e].jack[1][j];  //Mpi
        //                     fit_info.x[3][count][j] = jack_Mpi_MeV_exp[j];
        //                     fit_info.x[4][count][j] = l + 1e-6; //log
        //                     fit_info.x[5][count][j] = 0 + 1e-6;
        //                     fit_info.x[6][count][j] = OSTM + 1e-6;
        //                     fit_info.x[7][count][j] = w + 1.0; //w0
        //                 }
        //                 count++;
        //             }
        //         }
        //         fit_info.corr_id = { id0, id1 };
        //         fit_info.function = rhs_amu_pade;
        //         fit_info.repeat_start = 1;
        //         fit_info.linear_fit = false;
        //         if (l == 0 && stringOSTM[OSTM] == "TM") {
        //             fit_info.repeat_start = 0;
        //             fit_info.guess = { 2.03448e-08, 2.78364e-09, -2.01282e-05, -900 };
        //             // fit_info.guess = { 2.03448e-08, 2.78364e-09, 2.01282e-05, 10 };
        //             // fit_info.chi2_gap_jackboot=1e-2; // not implemented in this fit
        //             // fit_info.guess_per_jack=1;
        //             fit_info.h = { 2.05542e-10      ,3.70708e-11         , 1e-8     ,0.1 };
        //             // fit_info.NM = true;
        //             // fit_info.noderiv=true;
        //             fit_info.acc = 1e-2;
        //             fit_info.verbosity = 0;
        //             // fit_info.manual=true;

        //         }
        //         // fit_info.verbosity = 0;

        //         fit_info.covariancey = true;
        //         fit_info.compute_cov_fit(argv, jackall, lhs_amu);
        //         int ie = 0, ie1 = 0;
        //         for (int n = 0;n < fit_info.N;n++) {
        //             for (int e = 0;e < fit_info.Nxen[n].size();e++) {
        //                 ie1 = 0;
        //                 for (int n1 = 0;n1 < fit_info.N;n1++) {
        //                     for (int e1 = 0;e1 < fit_info.Nxen[n1].size();e1++) {
        //                         if (e != e1)   fit_info.cov[ie][ie1] = 0;
        //                         ie1++;
        //                     }
        //                 }
        //                 ie++;
        //             }
        //         }
        //         fit_info.compute_cov1_fit();

        //         std::string logname;
        //         if (l == 0) { logname = ""; }
        //         if (l == 1) { logname = "log1"; }
        //         if (l == 2) { logname = "log2"; }
        //         if (l == 3) { logname = "log3"; }

        //         if (l == 0 && w > 0) { fit_info.restore_default(); continue; }
        //         std::string wname;
        //         if (w == 0) { wname = "w1"; }
        //         if (w == 1) { wname = "w3"; }

        //         mysprintf(namefit, NAMESIZE, "amu_%s_pade_%s_%s_%s_cov", basename.c_str(), stringOSTM[OSTM].c_str(), logname.c_str(), wname.c_str());
        //         fit_result amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
        //         fit_info.band_range = { 0,0.0081 };
        //         std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
        //              fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };


        //         print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0001, xcont);
        //         syst_amu_W_lphys_Lref.add_fit(amu_SD_l_common_a4);

        //         // for (int j = 0;j < Njack;j++) {
        //         //     printf("jack %2d: chi2=%-15.2g", j, amu_SD_l_common_a4.chi2[j]);
        //         //     for (int p = 0;p < fit_info.Npar;p++) {
        //         //         printf("%-25.6g", amu_SD_l_common_a4.P[p][j]);
        //         //     }
        //         //     printf("\n");
        //         // }
        //         free_fit_result(fit_info, amu_SD_l_common_a4);
        //         fit_info.restore_default();

        //         // if (l == 0 && stringOSTM[OSTM] == "TM") { exit(1); }//stringOSTM[OSTM].c_str() == "OS" &&
        //     }
        // }


        ////////////////////////////// Ratio fit
        integrations = { "reinman" };
        for (auto integration : integrations) {
            int id0, id1;
            if (integration == "reinman") { id0 = ids[0 + iM * 2]; id1 = ids[1 + iM * 2]; }

            for (int a : { 1, 2}) {
                std::vector<int> vecl;
                if (a == 0) vecl = std::vector<int>({ 0, 5, 10, 15 });
                if (a == 1) vecl = std::vector<int>({ 0, 5, 10, 15 });
                // if (a == 2) vecl = std::vector<int>({  1, 2, 3, 4,  8, 12 });
                if (a == 2) vecl = std::vector<int>({ 1, 2, 3, 4, 5,  8, 10, 12, 15 });
                for (int l : vecl) {
                    for (int w = 0;w < 2;w++) {
                        for (int iR = 0; iR < 2;iR++) {
                            fit_info.restore_default();
                            fit_info.Npar = 2;
                            if (a >= 1) fit_info.Npar += 2;
                            if (l == 1 || l == 2 || l == 3 || l == 4 || l == 8 || l == 12) fit_info.Npar--;

                            if (integration == "reinman" && iR == 0) { id0 = ids[0 + iM * 2]; id1 = ids[1 + iM * 2]; }
                            if (integration == "reinman" && iR == 1) { id0 = ids[1 + iM * 2]; id1 = ids[0 + iM * 2]; }



                            fit_info.N = 2;
                            fit_info.Nvar = 8;
                            fit_info.Njack = Njack;
                            fit_info.myen = { B72_96, C06, D54 , E112 };
                            if (fit_info.Npar >= fit_info.myen.size() * fit_info.N) { fit_info.restore_default(); continue; }

                            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
                            count = 0;
                            for (int n = 0;n < fit_info.N;n++) {
                                for (int e : fit_info.myen) {
                                    for (int j = 0;j < Njack;j++) {
                                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
                                        fit_info.x[1][count][j] = jackall.en[e].jack[58][j];  // Delta_FV_GS
                                        fit_info.x[2][count][j] = jackall.en[e].jack[1][j];  //Mpi
                                        fit_info.x[3][count][j] = jack_Mpi_MeV_exp[j];
                                        fit_info.x[4][count][j] = l + 1e-6;
                                        fit_info.x[5][count][j] = a + 1e-6;
                                        fit_info.x[6][count][j] = 0 + 1e-6;
                                        fit_info.x[7][count][j] = w + 1e-6;
                                    }
                                    count++;
                                }
                            }
                            fit_info.corr_id = { id0 , id1 };
                            fit_info.function = rhs_amu_diff_ratio_charm;
                            fit_info.linear_fit = true;
                            fit_info.covariancey = true;
                            // fit_info.repeat_start=10;
                            // fit_info.acc= 1e-6;
                            // fit_info.chi2_gap_jackboot=0.1;
                            // fit_info.guess_per_jack=5;
                            // fit_info.repeat_start=5;
                            fit_info.verbosity = 0;
                            fit_info.compute_cov_fit(argv, jackall, lhs_amu_diff_ratio);
                            int ie = 0, ie1 = 0;
                            for (int n = 0;n < fit_info.N;n++) {
                                for (int e = 0;e < fit_info.myen.size();e++) {
                                    ie1 = 0;
                                    for (int n1 = 0;n1 < fit_info.N;n1++) {
                                        for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                                            if (e != e1)   fit_info.cov[ie][ie1] = 0;
                                            ie1++;
                                        }
                                    }
                                    ie++;
                                }
                            }
                            fit_info.compute_cov1_fit();

                            std::string logname;
                            if (l == 0) { logname = ""; }
                            if (l > 0) { logname = "log" + std::to_string(l); }
                            // if (l == 2) { logname = "log2"; }
                            // if (l == 3) { logname = "log3"; }

                            if (l == 0 && w > 0) continue;
                            std::string wname;
                            if (w == 0) { wname = "w1"; }
                            if (w == 1) { wname = "w3"; }



                            std::string aname;
                            if (a == 0) { aname = ""; }
                            if (a == 1) { aname = "a4"; }
                            if (a == 2) { aname = "+log"; }
                            if (a == 2 && l == 0) { fit_info.restore_default(); continue; }


                            std::string regname;
                            if (iR == 0) { regname = "R"; }
                            if (iR == 1) { regname = "R1"; }

                            mysprintf(namefit, NAMESIZE, "amu_%s_DR_%s_%s_%s_%s_cov", basename.c_str(), regname.c_str(), logname.c_str(), wname.c_str(), aname.c_str());
                            fit_result fit_DR = fit_all_data(argv, jackall, lhs_amu_diff_ratio, fit_info, namefit);
                            fit_info.band_range = { 0,0.0081 };
                            std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
                                 fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

                            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", fit_DR, fit_DR, 0, fit_info.myen.size() - 1, 0.0005, xcont);

                            if (fabs(fit_DR.P[0][Njack - 1]) - myres->comp_error(fit_DR.P[0]) <= 0 ||
                                fabs(fit_DR.P[1][Njack - 1]) - myres->comp_error(fit_DR.P[1]) <= 0 ||
                                std::isnan(fit_DR.P[0][Njack - 1]) ||
                                std::isnan(fit_DR.P[1][Njack - 1])
                                ) {
                                printf("fit DR produces poles excluding \n ");
                                fit_info.restore_default();
                                continue;

                            }

                            for (int j = 00;j < Njack;j++) {
                                fit_DR.P[0][j] = fit_DR.P[0][j] / fit_DR.P[1][j];
                            }
                            syst_amu_W_lphys_Lref.add_fit(fit_DR);

                            // if(l==4 && a==2){
                            //     exit(1);
                            // }

                            // for (int j = 0;j < Njack;j++) {
                            //     printf("jack %2d: chi2=%-15.2g", j, fit_DR.chi2[j]);
                            //     for (int p = 0;p < fit_info.Npar;p++) {
                            //         printf("%-25.6g", fit_DR.P[p][j]);
                            //     }
                            //     printf("\n");
                            // } 
                            // if ( strcmp(namefit,"amu_W_sphys_Mphi_DR_R1_log3_w1_+log_cov")==0 ) { exit(1); }

                            // if (l == 3 && w == 0 && a==2 && basename =="W_sphys_Mphi" && iR==1) { exit(1); }
                                                // if(iM==3 && a ==4 ) exit(1);
                            free_fit_result(fit_info, fit_DR);
                            fit_info.restore_default();
                        }


                    }
                }
            }
        }

        ////////////////////////////// a4 fit
        integrations = { "reinman" };
        for (auto integration : integrations) {
            int id0, id1;
            if (integration == "reinman") { id0 = ids[0 + iM * 2]; id1 = ids[1 + iM * 2]; }


            for (int a : { 0, 1, 2}) {
                std::vector<int> vec;
                for (int l : {0, 1, 2, 3, 4, 5, 8, 10, 12, 15}) {
                    if (a == 0) vec = std::vector<int>({ 0,  4, 8, 12 });
                    if (a == 1) vec = std::vector<int>({ 0,  1, 2, 3 });
                    if (a == 2) vec = std::vector<int>({ 0 });//vec = std::vector<int>({ 0,  5, 6, 7, 9, 10, 11,  13, 14, 15 });
                    for (int al : vec) {

                        for (int w = 0;w < 2;w++) {

                            if (al > 0 && l > 0) { fit_info.restore_default(); continue; }
                            fit_info.restore_default();
                            fit_info.Npar = 4;
                            if (a == 2)fit_info.Npar++;
                            if (integration == "reinman") { id0 = ids[0 + iM * 2]; id1 = ids[1 + iM * 2]; }



                            fit_info.N = 2;
                            fit_info.Nvar = 8;
                            fit_info.Njack = Njack;
                            fit_info.myen = { B72_96, C06, D54, E112 };
                            if (fit_info.Npar >= fit_info.myen.size() * fit_info.N) { fit_info.restore_default(); continue; }

                            fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
                            count = 0;
                            for (int n = 0;n < fit_info.N;n++) {
                                for (int e : fit_info.myen) {
                                    for (int j = 0;j < Njack;j++) {
                                        fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
                                        fit_info.x[1][count][j] = jackall.en[e].jack[58][j];  // Delta_FV_GS
                                        fit_info.x[2][count][j] = jackall.en[e].jack[1][j];  //Mpi
                                        fit_info.x[3][count][j] = jack_Mpi_MeV_exp[j];
                                        fit_info.x[4][count][j] = l + 1e-6;
                                        fit_info.x[5][count][j] = a + 1e-6;
                                        fit_info.x[6][count][j] = al + 1e-6;
                                        fit_info.x[7][count][j] = w + 1e-6;
                                    }
                                    count++;
                                }
                            }
                            fit_info.corr_id = { id0 , id1 };
                            fit_info.function = rhs_amu_a4;
                            fit_info.linear_fit = true;
                            fit_info.covariancey = true;
                            // fit_info.acc= 1e-6;
                            // fit_info.chi2_gap_jackboot=0.1;
                            // fit_info.guess_per_jack=5;
                            // fit_info.repeat_start=5;
                            fit_info.verbosity = 0;
                            fit_info.compute_cov_fit(argv, jackall, lhs_amu);
                            int ie = 0, ie1 = 0;
                            for (int n = 0;n < fit_info.N;n++) {
                                for (int e = 0;e < fit_info.myen.size();e++) {
                                    ie1 = 0;
                                    for (int n1 = 0;n1 < fit_info.N;n1++) {
                                        for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                                            if (e != e1)   fit_info.cov[ie][ie1] = 0;
                                            ie1++;
                                        }
                                    }
                                    ie++;
                                }
                            }
                            fit_info.compute_cov1_fit();

                            std::string logname;
                            if (l == 0) { logname = ""; }
                            if (l > 0) { logname = "log" + std::to_string(l); }


                            if (l == 0 && w > 0 && al == 0) continue;

                            std::string wname;
                            if (w == 0) { wname = "w1"; }
                            if (w == 1) { wname = "w3"; }

                            std::string logPname;
                            if (al == 0) { logPname = ""; }
                            if (al > 0) { logPname = "log" + std::to_string(al); }

                            std::string aname;
                            if (a == 0 && al == 0) { aname = "a4OS"; }
                            if (a == 0 && al > 0) { aname = "+" + logPname + "OS"; }
                            if (a == 1 && al == 0) { aname = "a4TM"; }
                            if (a == 1 && al > 0) { aname = "+" + logPname + "TM"; }
                            if (a == 2) { aname = "a4OS_a4TM"; }




                            mysprintf(namefit, NAMESIZE, "amu_%s_poly_%s_%s_%s_cov", basename.c_str(), logname.c_str(), wname.c_str(), aname.c_str());
                            fit_result amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
                            fit_info.band_range = { 0,0.0081 };
                            std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
                                 fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

                            print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);

                            syst_amu_W_lphys_Lref.add_fit(amu_SD_l_common_a4);

                            // if(iM==3 && a ==4 ) exit(1);
                            free_fit_result(fit_info, amu_SD_l_common_a4);
                            fit_info.restore_default();

                            // if (strcmp("amu_W_lphys_Lref_poly_log10_w3_a4TM_cov", namefit) == 0) { exit(1); }

                        }
                    }
                }
            }
        }


    }
    if (M.size() == 1)
        mysprintf(namefit, NAMESIZE, "Systematics_amu_%s_%s.txt", save_basename.c_str(), M[0].c_str());
    else
        mysprintf(namefit, NAMESIZE, "Systematics_amu_%s.txt", save_basename.c_str());
    compute_syst_eq28(syst_amu_W_lphys_Lref, argv[3], namefit);
    free(jack_Mpi_MeV_exp);
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
    fit_info.Npar = 3;
    fit_info.N = 2;
    fit_info.Nvar = 8;
    fit_info.Njack = Njack;
    fit_info.myen = { B72_64, C06, D54 , E112 };

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e : fit_info.myen) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
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
    fit_info.corr_id = { 146, 147 };
    fit_info.function = rhs_amu_common;
    fit_info.linear_fit = false;
    fit_info.covariancey = true;
    // fit_info.acc= 1e-6;
    // fit_info.chi2_gap_jackboot=0.1;
    // fit_info.guess_per_jack=5;
    // fit_info.repeat_start=5;
    fit_info.verbosity = 0;
    fit_info.compute_cov_fit(argv, jackall, lhs_amu);
    int ie = 0, ie1 = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            ie1 = 0;
            for (int n1 = 0;n1 < fit_info.N;n1++) {
                for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                    if (e != e1)   fit_info.cov[ie][ie1] = 0;
                    ie1++;
                }
            }
            ie++;
        }
    }
    fit_info.compute_cov1_fit();
    mysprintf(namefit, NAMESIZE, "amu_full_a2_MK_cov");
    fit_result amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
    fit_info.band_range = { 0,0.0081 };
    std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
         fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };


    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);

    free_fit_result(fit_info, amu_SD_l_common_a4);
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.restore_default();
    fit_info.Npar = 3;
    fit_info.N = 2;
    fit_info.Nvar = 8;
    fit_info.Njack = Njack;
    fit_info.myen = { B72_64, C06, D54 };

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e : fit_info.myen) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
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
    fit_info.corr_id = { 146, 147 };
    fit_info.function = rhs_amu_common;
    fit_info.linear_fit = false;
    fit_info.covariancey = true;
    // fit_info.acc= 1e-6;
    // fit_info.chi2_gap_jackboot=0.1;
    // fit_info.guess_per_jack=5;
    // fit_info.repeat_start=5;
    fit_info.verbosity = 0;
    fit_info.compute_cov_fit(argv, jackall, lhs_amu);
    ie = 0; ie1 = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            ie1 = 0;
            for (int n1 = 0;n1 < fit_info.N;n1++) {
                for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                    if (e != e1)   fit_info.cov[ie][ie1] = 0;
                    ie1++;
                }
            }
            ie++;
        }
    }
    fit_info.compute_cov1_fit();
    mysprintf(namefit, NAMESIZE, "amu_full_a2_MK_3beta_cov");
    amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
    fit_info.band_range = { 0,0.0081 };
    xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
         fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
    free_fit_result(fit_info, amu_SD_l_common_a4);

    ////////////////////////////////     a4  ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.restore_default();
    fit_info.Npar = 5;
    fit_info.N = 2;
    fit_info.Nvar = 8;
    fit_info.Njack = Njack;
    fit_info.myen = { B72_64, C06, D54 , E112 };

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e : fit_info.myen) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
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
    fit_info.corr_id = { 146, 147 };
    fit_info.function = rhs_amu_a4_common;
    fit_info.linear_fit = false;
    fit_info.covariancey = true;
    // fit_info.acc= 1e-6;
    // fit_info.chi2_gap_jackboot=0.1;
    // fit_info.guess_per_jack=5;
    // fit_info.repeat_start=5;
    fit_info.verbosity = 0;
    fit_info.compute_cov_fit(argv, jackall, lhs_amu);
    ie = 0; ie1 = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            ie1 = 0;
            for (int n1 = 0;n1 < fit_info.N;n1++) {
                for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                    if (e != e1)   fit_info.cov[ie][ie1] = 0;
                    ie1++;
                }
            }
            ie++;
        }
    }
    fit_info.compute_cov1_fit();
    mysprintf(namefit, NAMESIZE, "amu_full_a2_a4_MK_cov");
    amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
    fit_info.band_range = { 0,0.0081 };
    xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
         fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
    free_fit_result(fit_info, amu_SD_l_common_a4);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.restore_default();
    fit_info.Npar = 5;
    fit_info.N = 2;
    fit_info.Nvar = 8;
    fit_info.Njack = Njack;
    fit_info.myen = { B72_64, C06, D54 };

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e : fit_info.myen) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
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
    fit_info.corr_id = { 146, 147 };
    fit_info.function = rhs_amu_a4_common;
    fit_info.linear_fit = false;
    fit_info.covariancey = true;
    // fit_info.acc= 1e-6;
    // fit_info.chi2_gap_jackboot=0.1;
    // fit_info.guess_per_jack=5;
    // fit_info.repeat_start=5;
    fit_info.verbosity = 0;
    fit_info.compute_cov_fit(argv, jackall, lhs_amu);
    ie = 0; ie1 = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            ie1 = 0;
            for (int n1 = 0;n1 < fit_info.N;n1++) {
                for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                    if (e != e1)   fit_info.cov[ie][ie1] = 0;
                    ie1++;
                }
            }
            ie++;
        }
    }
    fit_info.compute_cov1_fit();
    mysprintf(namefit, NAMESIZE, "amu_full_a2_a4_MK_3beta_cov");
    amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
    fit_info.band_range = { 0,0.0081 };
    xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
         fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
    free_fit_result(fit_info, amu_SD_l_common_a4);



    ////////////////////////////////    s W  ///////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.restore_default();
    fit_info.Npar = 3;
    fit_info.N = 2;
    fit_info.Nvar = 8;
    fit_info.Njack = Njack;
    fit_info.myen = { B72_64, C06, D54 , E112 };

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e : fit_info.myen) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
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
    fit_info.corr_id = { 151, 152 };
    fit_info.function = rhs_amu_common;
    fit_info.linear_fit = false;
    fit_info.covariancey = true;
    // fit_info.acc= 1e-6;
    // fit_info.chi2_gap_jackboot=0.1;
    // fit_info.guess_per_jack=5;
    // fit_info.repeat_start=5;
    fit_info.verbosity = 0;
    fit_info.compute_cov_fit(argv, jackall, lhs_amu);
    ie = 0; ie1 = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            ie1 = 0;
            for (int n1 = 0;n1 < fit_info.N;n1++) {
                for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                    if (e != e1)   fit_info.cov[ie][ie1] = 0;
                    ie1++;
                }
            }
            ie++;
        }
    }
    fit_info.compute_cov1_fit();
    mysprintf(namefit, NAMESIZE, "amu_W_a2_MK_cov");
    amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
    fit_info.band_range = { 0,0.0081 };
    xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
         fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
    free_fit_result(fit_info, amu_SD_l_common_a4);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.restore_default();
    fit_info.Npar = 3;
    fit_info.N = 2;
    fit_info.Nvar = 8;
    fit_info.Njack = Njack;
    fit_info.myen = { B72_64, C06, D54 };

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e : fit_info.myen) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
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
    fit_info.corr_id = { 151, 152 };
    fit_info.function = rhs_amu_common;
    fit_info.linear_fit = false;
    fit_info.covariancey = true;
    // fit_info.acc= 1e-6;
    // fit_info.chi2_gap_jackboot=0.1;
    // fit_info.guess_per_jack=5;
    // fit_info.repeat_start=5;
    fit_info.verbosity = 0;
    fit_info.compute_cov_fit(argv, jackall, lhs_amu);
    ie = 0; ie1 = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            ie1 = 0;
            for (int n1 = 0;n1 < fit_info.N;n1++) {
                for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                    if (e != e1)   fit_info.cov[ie][ie1] = 0;
                    ie1++;
                }
            }
            ie++;
        }
    }
    fit_info.compute_cov1_fit();
    mysprintf(namefit, NAMESIZE, "amu_W_a2_MK_3beta_cov");
    amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
    fit_info.band_range = { 0,0.0081 };
    xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
         fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
    free_fit_result(fit_info, amu_SD_l_common_a4);

    ////////////////////////////////    s W  ///////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.restore_default();
    fit_info.Npar = 5;
    fit_info.N = 2;
    fit_info.Nvar = 8;
    fit_info.Njack = Njack;
    fit_info.myen = { B72_96, C06, D54 , E112 };

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e : fit_info.myen) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
                fit_info.x[1][count][j] = jackall.en[e].jack[58][j];  // Delta_FV_GS
                fit_info.x[2][count][j] = jackall.en[e].jack[1][j];  //Mpi
                fit_info.x[3][count][j] = jack_Mpi_MeV_exp[j];
                fit_info.x[4][count][j] = 0/* l */ + 1e-6;
                fit_info.x[5][count][j] = 2/* a */ + 1e-6;
                fit_info.x[6][count][j] = 0 + 1e-6;
                fit_info.x[7][count][j] = 0/* w */ + 1e-6;
            }
            count++;
        }
    }
    fit_info.corr_id = { 151, 152 };
    fit_info.function = rhs_amu_a4_common;
    fit_info.linear_fit = false;
    fit_info.covariancey = true;
    // fit_info.acc= 1e-6;
    // fit_info.chi2_gap_jackboot=0.1;
    // fit_info.guess_per_jack=5;
    // fit_info.repeat_start=5;
    fit_info.verbosity = 0;
    fit_info.compute_cov_fit(argv, jackall, lhs_amu);
    ie = 0; ie1 = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            ie1 = 0;
            for (int n1 = 0;n1 < fit_info.N;n1++) {
                for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                    if (e != e1)   fit_info.cov[ie][ie1] = 0;
                    ie1++;
                }
            }
            ie++;
        }
    }
    fit_info.compute_cov1_fit();
    mysprintf(namefit, NAMESIZE, "amu_W_a2_a4_MK_cov");
    amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
    fit_info.band_range = { 0,0.0081 };
    xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
         fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
    free_fit_result(fit_info, amu_SD_l_common_a4);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    fit_info.restore_default();
    fit_info.Npar = 5;
    fit_info.N = 2;
    fit_info.Nvar = 8;
    fit_info.Njack = Njack;
    fit_info.myen = { B72_96, C06, D54 };

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e : fit_info.myen) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
                fit_info.x[1][count][j] = jackall.en[e].jack[58][j];  // Delta_FV_GS
                fit_info.x[2][count][j] = jackall.en[e].jack[1][j];  //Mpi
                fit_info.x[3][count][j] = jack_Mpi_MeV_exp[j];
                fit_info.x[4][count][j] = 0/* l */ + 1e-6;
                fit_info.x[5][count][j] = 2/* a */ + 1e-6;
                fit_info.x[6][count][j] = 0 + 1e-6;
                fit_info.x[7][count][j] = 0/* w */ + 1e-6;
            }
            count++;
        }
    }
    fit_info.corr_id = { 151, 152 };
    fit_info.function = rhs_amu_a4_common;
    fit_info.linear_fit = false;
    fit_info.covariancey = true;
    // fit_info.acc= 1e-6;
    // fit_info.chi2_gap_jackboot=0.1;
    // fit_info.guess_per_jack=5;
    // fit_info.repeat_start=5;
    fit_info.verbosity = 0;
    fit_info.compute_cov_fit(argv, jackall, lhs_amu);
    ie = 0; ie1 = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e = 0;e < fit_info.myen.size();e++) {
            ie1 = 0;
            for (int n1 = 0;n1 < fit_info.N;n1++) {
                for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                    if (e != e1)   fit_info.cov[ie][ie1] = 0;
                    ie1++;
                }
            }
            ie++;
        }
    }
    fit_info.compute_cov1_fit();
    mysprintf(namefit, NAMESIZE, "amu_W_a2_a4_MK_3beta_cov");
    amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
    fit_info.band_range = { 0,0.0081 };
    xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
         fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
    free_fit_result(fit_info, amu_SD_l_common_a4);


    //////////////////////////////////////////////////////////////
    // charm
    //////////////////////////////////////////////////////////////

    ////////////////////////// fit the charm A to a constant

    fit_info.restore_default();
    fit_info.corr_id = { 156, 157, 161, 162 };
    fit_info.N = fit_info.corr_id.size();
    fit_info.Nvar = 1;
    fit_info.Npar = fit_info.corr_id.size();
    fit_info.Njack = Njack;
    fit_info.myen = { A30, A40, A53 };
    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.N;n++) {
        for (int e : fit_info.myen) {
            for (int j = 0;j < Njack;j++) {
                fit_info.x[0][count][j] = jackall.en[e].header.L * jackall.en[e].jack[1][j];// ;//jackall.en[e].jack[41][j] * jack_Mpi_MeV_exp[j] / hbarc;//* jackall.en[e].jack[1][j] 
            }
            count++;
        }
    }
    fit_info.function = const_A;
    fit_info.linear_fit = true;

    mysprintf(namefit, NAMESIZE, "interpolate_A_charm");
    fit_result Acharm = fit_all_data(argv, jackall, lhs_Acharm, fit_info, namefit);
    fit_info.band_range = { 0, 8 };
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "Acharm", Acharm, Acharm, 0, fit_info.myen.size() - 1, 0.1);
    for (int j = 0; j < Njack;j++) {
        double a = jackall.en[B72_64].jack[41][j];
        double Mpi_MeV = jack_Mpi_MeV_exp[j];
        double L = jackall.en[B72_64].header.L;
        for (int i = 0; i < fit_info.N; i++)
            jackall.en[A30].jack[fit_info.corr_id[i]][j] = Acharm.P[i][j];

    }
    free_fit_result(fit_info, Acharm);
    fit_info.restore_default();


    ///////////////////////////////////////////////////////////////////////

    ////////////////////////////////     a4  ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    {
        fit_info.restore_default();
        fit_info.Npar = 5;
        fit_info.N = 2;
        fit_info.Nvar = 8;
        fit_info.Njack = Njack;
        fit_info.myen = { A30,  B72_64, C06, D54 , E112 };

        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
        count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e : fit_info.myen) {
                for (int j = 0;j < Njack;j++) {
                    fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
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
        fit_info.corr_id = { 156, 161 };
        fit_info.function = rhs_amu_a4_common;
        fit_info.linear_fit = false;
        fit_info.covariancey = true;
        // fit_info.acc= 1e-6;
        // fit_info.chi2_gap_jackboot=0.1;
        // fit_info.guess_per_jack=5;
        // fit_info.repeat_start=5;
        fit_info.verbosity = 0;
        fit_info.compute_cov_fit(argv, jackall, lhs_amu);
        ie = 0; ie1 = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < fit_info.myen.size();e++) {
                ie1 = 0;
                for (int n1 = 0;n1 < fit_info.N;n1++) {
                    for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                        if (e != e1)   fit_info.cov[ie][ie1] = 0;
                        ie1++;
                    }
                }
                ie++;
            }
        }
        fit_info.compute_cov1_fit();
        mysprintf(namefit, NAMESIZE, "amu_full_a2_a4_MJPsi_cov");
        amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
        fit_info.band_range = { 0,0.0085 };
        xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
             fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

        //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
        free_fit_result(fit_info, amu_SD_l_common_a4);


        ///////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////

        fit_info.restore_default();
        fit_info.Npar = 5;
        fit_info.N = 2;
        fit_info.Nvar = 8;
        fit_info.Njack = Njack;
        fit_info.myen = {  A30, B72_64, C06, D54 };

        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.myen.size() * fit_info.N, fit_info.Njack);
        count = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e : fit_info.myen) {
                for (int j = 0;j < Njack;j++) {
                    fit_info.x[0][count][j] = pow(jackall.en[e].jack[41][j], 2); // a^2
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
        fit_info.corr_id = { 156, 161 };
        fit_info.function = rhs_amu_a4_common;
        fit_info.linear_fit = false;
        fit_info.covariancey = true;
        // fit_info.acc= 1e-6;
        // fit_info.chi2_gap_jackboot=0.1;
        // fit_info.guess_per_jack=5;
        // fit_info.repeat_start=5;
        fit_info.verbosity = 0;
        fit_info.compute_cov_fit(argv, jackall, lhs_amu);
        ie = 0; ie1 = 0;
        for (int n = 0;n < fit_info.N;n++) {
            for (int e = 0;e < fit_info.myen.size();e++) {
                ie1 = 0;
                for (int n1 = 0;n1 < fit_info.N;n1++) {
                    for (int e1 = 0;e1 < fit_info.myen.size();e1++) {
                        if (e != e1)   fit_info.cov[ie][ie1] = 0;
                        ie1++;
                    }
                }
                ie++;
            }
        }
        fit_info.compute_cov1_fit();
        mysprintf(namefit, NAMESIZE, "amu_full_a2_a4_MJPsi_4beta_cov");
        amu_SD_l_common_a4 = fit_all_data(argv, jackall, lhs_amu, fit_info, namefit);
        fit_info.band_range = { 0,0.0085 };
        xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
             fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };

        //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "afm", amu_SD_l_common_a4, amu_SD_l_common_a4, 0, fit_info.myen.size() - 1, 0.0005, xcont);
        free_fit_result(fit_info, amu_SD_l_common_a4);
    }

}