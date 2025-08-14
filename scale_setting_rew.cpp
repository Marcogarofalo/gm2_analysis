#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
// #include <mdspan>

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
    C20,
    A12
};


constexpr double fpi_MeV = 130.5;
constexpr double fpi_MeV_err = 0.04;

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
        printf("amu(%s)=%-20.12g %-20.12g\n", lattices[l].c_str(), amu[Njack - 1], myres->comp_error(amu.data()));
        fprintf(famul, "%s %-20.12g %-20.12g\n", lattices[l].c_str(), amu[Njack - 1], myres->comp_error(amu.data()));

        char file_amul_jack[NAMESIZE];
        mysprintf(file_amul_jack, NAMESIZE, "%s/%s_amul_jack_%s.txt", argv[3], namefit, lattices[l].c_str());
        myres->write_jack_in_file(amu.data(), file_amul_jack);

        count += fit_info.Nxen[l].size();
    }
}

fit_result read_file_P(std::string file_der) {

    std::string line;
    std::ifstream f_der(file_der);

    fit_result fit_out;
    std::vector<std::string>  what = { "npar", "chi2dof" };
    for (int i = 0;i < 2;i++) {
        if (std::getline(f_der, line)) {
            std::vector<std::string> array_l = split(line, ' ');
            error(array_l.size() != 2, 1, "main", "error reading %s", file_der.c_str());
            if (i == 0)        fit_out.Npar = stoi(array_l[1]);
            else fit_out.chi2[myres->Njack - 1] = stod(array_l[1]);
            error(array_l[0].compare(what[i]) != 0, 1, "main", "error reading %s", file_der.c_str());
        }
    }
    std::vector<double> par(fit_out.Npar);
    std::vector<double> err_p(fit_out.Npar);
    for (int i = 0;i < fit_out.Npar;i++) {
        if (std::getline(f_der, line)) {
            std::vector<std::string> array_l = split(line, ' ');
            // printf("read %s    %d\n", line.c_str(), (int)array_l.size());
            error(array_l.size() != 3, 1, "main", "error reading %s", file_der.c_str());
            par[i] = stod(array_l[1]);
            err_p[i] = stod(array_l[2]);
        }
    }
    std::vector<std::vector<double>> cov(fit_out.Npar, std::vector<double>(fit_out.Npar));
    for (int i = 0;i < fit_out.Npar;i++) {
        if (std::getline(f_der, line)) {
            std::vector<std::string> array_l = split(line, ' ');
            for (int j = 0;j < fit_out.Npar;j++) {
                error(array_l.size() != fit_out.Npar, 1, "main", "error reading %s", file_der.c_str());
                cov[i][j] = stod(array_l[j]) * err_p[i] * err_p[j];
            }
        }
    }
    // Convert cov to double**
    double** cov_ptr = new double* [fit_out.Npar];
    for (int i = 0; i < fit_out.Npar; ++i) {
        cov_ptr[i] = cov[i].data();
    }
    fit_out.P = myres->create_fake_covariance(par.data(), fit_out.Npar, cov_ptr, 1);
    delete[] cov_ptr;
    f_der.close();
    return fit_out;
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

    mysprintf(namefile, NAMESIZE, "%s/%s_cA.12.48_mu.0.001200", argv[2], argv[1]);
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



    //////////////////////////////////////////////////////////////
    // read mpcac
    //////////////////////////////////////////////////////////////
    FILE* f_mpcac = open_file("/home/garofalo/analysis/g-2_new_stat/mpcac/mpcac_over_tau.txt", "r+");
    std::vector<double*> vev_mpcac(files.size());
    std::vector<char[NAMESIZE]> ens_name(files.size());
    for (int e = 0; e < files.size();e++) {
        double mean, err;
        fscanf(f_mpcac, "%s  %lf  %lf\n", ens_name[e], &mean, &err);
        vev_mpcac[e] = myres->create_fake(mean, err, e + 50);
        printf("%s  %lf \n", ens_name[e], (vev_mpcac[e][Njack - 1]));
    }
    //////////////////////////////////////////////////////////////
    // fix ZA
    //////////////////////////////////////////////////////////////
    int seed = 300;
    jackall.en[A53].jack[23] = myres->create_fake(0.728393, 0.0018, seed);
    jackall.en[A40].jack[23] = myres->create_fake(0.728393, 0.0018, seed);
    jackall.en[A30].jack[23] = myres->create_fake(0.728393, 0.0018, seed);
    jackall.en[A12].jack[23] = myres->create_fake(0.728393, 0.0018, seed);
    jackall.en[B14_64].jack[23] = jackall.en[B72_64].jack[23];
    jackall.en[B25_48].jack[23] = jackall.en[B72_64].jack[23];
    jackall.en[C20].jack[23] = jackall.en[C06].jack[23];

    //////////////////////////////////////////////////////////////
    // extra jack
    //////////////////////////////////////////////////////////////
    data_all jackextra;
    jackextra.resampling = jackall.resampling;

    jackextra.en = new data_single[jackall.ens];
    jackextra.ens = jackall.ens;
    for (int count = 0; count < jackextra.ens; count++) {

        data_single dj;
        dj.header = jackall.en[count].header;
        dj.Nobs = jackall.en[0].Nobs + 10 + 3 * 2 + 1;
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
    int id_fpi_cor = jackall.en[0].Nobs + 0;
    int id_fpi_cor_mu1 = jackall.en[0].Nobs + 1;
    int id_Mpi_cor = jackall.en[0].Nobs + 2;
    int id_Mpi_cor_mu1 = jackall.en[0].Nobs + 3;

    int id_fpi_rew = jackall.en[0].Nobs + 4;
    int id_Mpi_rew = jackall.en[0].Nobs + 5;

    double* dsfpi;
    double* dsMpi;
    double* dms;

    //////////////////////////////////////////////////////////////
    // read iso and sim params
    //////////////////////////////////////////////////////////////

    char namefile_plateaux[NAMESIZE];
    mysprintf(namefile_plateaux, NAMESIZE, "plateaux.txt");

    char** option;
    option = (char**)malloc(sizeof(char*) * 7);
    option[0] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[1] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[2] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[3] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[4] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[5] = (char*)malloc(sizeof(char) * NAMESIZE);
    option[6] = (char*)malloc(sizeof(char) * NAMESIZE);

    mysprintf(option[1], NAMESIZE, "read_plateaux"); // blind/see/read_plateaux
    mysprintf(option[2], NAMESIZE, "-p");            // -p
    mysprintf(option[3], NAMESIZE, "/home/garofalo/analysis/flow/w0_analysis/");         // path
    mysprintf(option[4], NAMESIZE, argv[1]);         // resampling
    mysprintf(option[5], NAMESIZE, "no");            // pdf
    mysprintf(option[6], NAMESIZE, "file_ensemble");         // infile

    std::vector<double*> amuiso(3);
    std::vector<std::vector<int>> id_amuiso(files.size(), std::vector<int>(5));
    std::vector<std::vector<int>> id_amusim(files.size(), std::vector<int>(5));
    std::vector<int> id_a_fm(5);

    // auto span_iso = std::mdspan(id_amuiso.data(), files.size(), 3);
    // auto span_sim = std::mdspan(id_amuiso.data(), files.size(), 3);
    // read the fit result
    std::string file_der = "/home/garofalo/analysis/flow/data/fit_all_beta/der_fpi_diffplat_full_mu_mu2_mu3_chi2d1_fit_P.dat";
    fit_result fit_der_fpi = read_file_P(file_der);
    std::string file_der_M = "/home/garofalo/analysis/flow/data/fit_all_beta/der_Mpi_diffplat_full_mu_mu2_mu3_chi2d1_fit_P.dat";
    fit_result fit_der_Mpi = read_file_P(file_der_M);

    // for (int i = 0;i < npar[0];i++) {
    //     printf("par[%d]=%g +- %g\n", i, Pj[i][Njack - 1], myres->comp_error(Pj[i]));
    // }

    double mean, err;
    double* tif = (double*)malloc(sizeof(double) * fit_der_fpi.Npar);
    double* tmp_x = (double*)malloc(sizeof(double) * 1);
    double* muder_s_MeV = (double*)malloc(sizeof(double) * myres->Njack);

    for (int e = 0; e < files.size(); e++) {
        printf("ensemble %s", files[e].c_str());
        if (e == B72_64 || e == B72_96 || e == B14_64 || e == B25_48) {
            mysprintf(option[6], NAMESIZE, "onlinemeas_B64.dat");
        }
        else if (e == C06 || e == C112 || e == C20) {
            mysprintf(option[6], NAMESIZE, "onlinemeas_C80_LMA.dat");
        }
        else if (e == D54) {
            mysprintf(option[6], NAMESIZE, "onlinemeas_D96.dat");
        }
        else if (e == E112) {
            mysprintf(option[6], NAMESIZE, "onlinemeas_E112_LMA.dat");
        }
        else {
            // for ensemble A we took the maximum
            dms = myres->create_fake(0.0006, 1e-16, 1);
            // mysprintf(option[6], NAMESIZE, "onlinemeas_A12.dat");
            continue; // skip A ensembles
        }
        id_amusim[e][0] = jackall.en[e].Nobs + 6;
        id_amusim[e][1] = jackall.en[e].Nobs + 7;
        id_amusim[e][2] = jackall.en[e].Nobs + 8;

        id_amuiso[e][0] = jackall.en[e].Nobs + 9;
        id_amuiso[e][1] = jackall.en[e].Nobs + 10;
        id_amuiso[e][2] = jackall.en[e].Nobs + 11;

        id_a_fm[e] = jackall.en[e].Nobs + 12;

        line_read_param(option, "muliso", mean, err, seed, namefile_plateaux);
        free(jackextra.en[e].jack[id_amuiso[e][0]]);
        jackextra.en[e].jack[id_amuiso[e][0]] = myres->create_fake(mean, err, seed);

        line_read_param(option, "musiso", mean, err, seed, namefile_plateaux);
        free(jackextra.en[e].jack[id_amuiso[e][1]]);
        jackextra.en[e].jack[id_amuiso[e][1]] = myres->create_fake(mean, err, seed);

        line_read_param(option, "muciso", mean, err, seed, namefile_plateaux);
        free(jackextra.en[e].jack[id_amuiso[e][2]]);
        jackextra.en[e].jack[id_amuiso[e][2]] = myres->create_fake(mean, err, seed);



        line_read_param(option, "mulsim", mean, err, seed, namefile_plateaux);
        free(jackextra.en[e].jack[id_amusim[e][0]]);
        jackextra.en[e].jack[id_amusim[e][0]] = myres->create_fake(mean, err, seed);

        line_read_param(option, "mussim", mean, err, seed, namefile_plateaux);
        free(jackextra.en[e].jack[id_amusim[e][1]]);
        jackextra.en[e].jack[id_amusim[e][1]] = myres->create_fake(mean, err, seed);

        line_read_param(option, "mucsim", mean, err, seed, namefile_plateaux);
        free(jackextra.en[e].jack[id_amusim[e][2]]);
        jackextra.en[e].jack[id_amusim[e][2]] = myres->create_fake(mean, err, seed);

        line_read_param(option, "mucsim", mean, err, seed, namefile_plateaux);
        free(jackextra.en[e].jack[id_amusim[e][2]]);
        jackextra.en[e].jack[id_amusim[e][2]] = myres->create_fake(mean, err, seed);

        line_read_param(option, "a", mean, err, seed, namefile_plateaux);
        free(jackextra.en[e].jack[id_a_fm[e]]);
        jackextra.en[e].jack[id_a_fm[e]] = myres->create_fake(mean, err, seed);


        // printf("ensemble %s: iso: %g %g %g  sim: %g %g %g\n",
        //     files[e].c_str(),
        //     jackextra.en[e].jack[id_amuiso[e][0]][Njack - 1],
        //     jackextra.en[e].jack[id_amuiso[e][1]][Njack - 1],
        //     jackextra.en[e].jack[id_amuiso[e][2]][Njack - 1],
        //     jackextra.en[e].jack[id_amusim[e][0]][Njack - 1],
        //     jackextra.en[e].jack[id_amusim[e][1]][Njack - 1],
        //     jackextra.en[e].jack[id_amusim[e][2]][Njack - 1]);

        auto fit_fun = [](int n, int Nvar, double* x, int Npar, double* P) {
            double amu = x[0];
            double r = P[0] / amu + P[1] / (amu * amu) + P[2] / (amu * amu * amu);
            return r;
            };

        for (int j = 0; j < Njack; j++) {
            double amul_iso = jackextra.en[e].jack[id_amuiso[e][0]][j];
            double amus_iso = jackextra.en[e].jack[id_amuiso[e][1]][j];
            double amuc_iso = jackextra.en[e].jack[id_amuiso[e][2]][j];
            double amul_sim = jackextra.en[e].jack[id_amusim[e][0]][j];
            double amus_sim = jackextra.en[e].jack[id_amusim[e][1]][j];
            double amuc_sim = jackextra.en[e].jack[id_amusim[e][2]][j];
            double a_fm = jackextra.en[e].jack[id_a_fm[e]][j];

            tmp_x[0] = amus_sim / amul_iso;
            double dmu =  (amus_iso - amus_sim);
            for (int p = 0; p < fit_der_fpi.Npar; p++)
                tif[p] = fit_der_fpi.P[p][j];
            muder_s_MeV[j] = fit_fun(0, 0, tmp_x, fit_der_fpi.Npar, tif) * (a_fm / hbarc) / amul_iso;

            jackextra.en[e].jack[id_fpi_rew][j] = jackextra.en[e].jack[163][j] + dmu * muder_s_MeV[j];

            for (int p = 0; p < fit_der_Mpi.Npar; p++)
                tif[p] = fit_der_Mpi.P[p][j];
            muder_s_MeV[j] = fit_fun(0, 0, tmp_x, fit_der_Mpi.Npar, tif) * (a_fm / hbarc) / amul_iso;

            jackextra.en[e].jack[id_Mpi_rew][j] = jackextra.en[e].jack[1][j] + dmu * muder_s_MeV[j];

            tmp_x[0] = amuc_sim / amul_iso;
            dmu =  (amuc_iso - amuc_sim);
            
            for (int p = 0; p < fit_der_fpi.Npar; p++)
                tif[p] = fit_der_fpi.P[p][j];
            muder_s_MeV[j] = fit_fun(0, 0, tmp_x, fit_der_fpi.Npar, tif) * (a_fm / hbarc) / amul_iso;

            jackextra.en[e].jack[id_fpi_rew][j] +=  dmu * muder_s_MeV[j];

            for (int p = 0; p < fit_der_Mpi.Npar; p++)
                tif[p] = fit_der_Mpi.P[p][j];
            muder_s_MeV[j] = fit_fun(0, 0, tmp_x, fit_der_Mpi.Npar, tif) * (a_fm / hbarc) / amul_iso;

            jackextra.en[e].jack[id_Mpi_rew][j] += dmu * muder_s_MeV[j];


        }


    }
    free(tif);
    free(tmp_x);
    free(muder_s_MeV);

    //////////////////////////////////////////////////////////////
    // gm2 way of correcting fpi and Mpi
    //////////////////////////////////////////////////////////////

    for (int e = 0; e < files.size(); e++) {
        if (e == B72_64 || e == B72_96) {
            // dsfpi = myres->create_fake(0.0655, 0.11, 2000 + B72_64);
            // dsMpi = myres->create_fake(0.0825, 0.0682, 2000 + 1 + B72_64);
            dms = myres->create_fake(0.00042039860, 1e-16, 1);
        }
        else if (e == C06 || e == C112) {
            // dsfpi = myres->create_fake(0.31, 0.0588, 2000 + C06);
            // dsMpi = myres->create_fake(-0.0732, 0.0335, 2000 + 1 + C06);
            dms = myres->create_fake(-0.00059425330, 1e-16, 1);

        }
        else if (e == D54) {
            // dsfpi = myres->create_fake(0.171, 0.0635, 2000 + D54);
            // dsMpi = myres->create_fake(0.111, 0.0568, 2000 + 1 + D54);
            dms = myres->create_fake(0.00016284, 1e-16, 1);

        }
        else if (e == E112) {
            // dsfpi = myres->create_fake(0.127, 0.04, 2000 + E112);
            // dsMpi = myres->create_fake(0.0204, 0.03, 2000 + 1 + E112);
            dms = myres->create_fake(0.00029011518, 1e-16, 1);

        }
        else {
            // dsfpi = myres->create_fake(0, 1e-12, 2000 + B72_64);
            // dsMpi = myres->create_fake(0, 1e-12, 2000 + 1 + B72_64);
            // for ensemble A we took the maximum
            dms = myres->create_fake(0.0006, 1e-16, 1);
        }
        dsMpi = myres->create_fake(0.035175, 1e-12, 2000 + 1 + B72_64);
        dsfpi = myres->create_fake(0.174, 0.032, 2000 + B72_64);
        myres->mult(dsfpi, dsfpi, dms);
        myres->mult(dsMpi, dsMpi, dms);

        myres->add_error_quadrature(jackextra.en[e].jack[id_fpi_cor], jackall.en[e].jack[163], dsfpi[Njack - 1]);
        myres->add_error_quadrature(jackextra.en[e].jack[id_fpi_cor_mu1], jackall.en[e].jack[164], dsfpi[Njack - 1]);
        myres->add_error_quadrature(jackextra.en[e].jack[id_Mpi_cor], jackall.en[e].jack[1], dsMpi[Njack - 1]);
        myres->add_error_quadrature(jackextra.en[e].jack[id_Mpi_cor_mu1], jackall.en[e].jack[123], dsMpi[Njack - 1]);


        // iter 0 is done with jackall
        myres->add(jackextra.en[e].jack[163], jackall.en[e].jack[163], dsfpi);
        myres->add(jackextra.en[e].jack[164], jackall.en[e].jack[164], dsfpi);
        myres->add(jackextra.en[e].jack[1], jackall.en[e].jack[1], dsMpi);
        myres->add(jackextra.en[e].jack[123], jackall.en[e].jack[123], dsMpi);
        free(dsfpi);
        free(dsMpi);
        free(dms);
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////


    std::vector<double> xcont = {};

    //////////////////////////////////////////////////////////////
    // A12 no C20 , unitary only, iter 0
    //////////////////////////////////////////////////////////////

    fit_info.restore_default();
    // fit_info.N = 8;
    fit_info.Nvar = 10;
    fit_info.Npar = 7;
    fit_info.Njack = Njack;
    fit_info.Nxen = { {A53, A40, A30, A12}, {B25_48, B14_64, B72_64, B72_96},
                        { C06, C112}, {D54},
                        {E112} };
    fit_info.init_N_etot_form_Nxen();

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.Nxen.size();n++) {
        for (int e : fit_info.Nxen[n]) {
            for (int j = 0;j < Njack;j++) {
                double my_mu, my_M, my_fpi;

                my_mu = jackall.en[e].jack[165][j];
                my_M = jackall.en[e].jack[1][j];
                my_fpi = jackall.en[e].jack[163][j];

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

                fit_info.x[8][count][j] = vev_mpcac[e][j];// mpcac/mu
                fit_info.x[9][count][j] = jackall.en[e].jack[23][j];// Z_A

                double mpcac_mu = fit_info.x[8][count][j];
                double ZA = fit_info.x[9][count][j];

                double mr = 0;
                mr = ZA * mpcac_mu;
                double cl = sqrt(1 + mr * mr);
                xi /= cl * cl * cl;
                fit_info.x[4][count][j] = xi;
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
    fit_info.corr_id = { 163, 163, 163, 163, 163 }; //  fpi(mu1), fpi(mu2)
    fit_info.function = rhs_afpi;
    fit_info.linear_fit = false;
    // fit_info.covariancey = true;
    // fit_info.acc = 1e-10;
    // fit_info.h = 1e-7;
    // fit_info.maxiter = 500;
    // fit_info.NM = true;
    // fit_info.chi2_gap_jackboot = 0.01;
    // fit_info.guess_per_jack = 5;
    // fit_info.repeat_start = 100;
    // // fit_info.verbosity = 1;
    // fit_info.compute_cov_fit(argv, jackall, lhs_afpi_max_twist);
    // int ide = 0, ide1 = 0;
    // for (int n = 0;n < fit_info.Nxen.size();n++) {
    //     for (int e : fit_info.Nxen[n]) {
    //         ide1 = 0;
    //         for (int n1 = 0;n1 < fit_info.Nxen.size();n1++) {
    //             for (int e1 : fit_info.Nxen[n1]) {
    //                 if (e != e1)   fit_info.cov[ide][ide1] = 0;
    //                 printf("%-12.5g ", fit_info.cov[ide][ide1] / sqrt(fit_info.cov[ide][ide] * fit_info.cov[ide1][ide1]));
    //                 ide1++;
    //             }
    //         }
    //         printf("\n");
    //         ide++;
    //     }
    // }
    // fit_info.compute_cov1_fit();
    fit_info.guess = { 0.0908026, 0.07951 , 0.06816, 0.05688, 0.04891 ,-8.75 ,-500 };
    mysprintf(namefit, NAMESIZE, "afpi_max_twist_A12_noC20_cov_unitary_iter0");
    fit_result fit_afpi_max_twist_A12_noC20_corr_unitary_iter0 = fit_all_data(argv, jackall, lhs_afpi_max_twist, fit_info, namefit);
    fit_info.band_range = { 0.005,0.035 };
    // std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
    //      fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };


    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    // std::vector<double> xcont = {};
    print_fit_band(argv, jackall, fit_info, fit_info, namefit, "xi", fit_afpi_max_twist_A12_noC20_corr_unitary_iter0, fit_afpi_max_twist_A12_noC20_corr_unitary_iter0, 4, 0, 0.0005, xcont);

    myres->write_jack_in_file(fit_afpi_max_twist_A12_noC20_corr_unitary_iter0.P[0], "../../g-2_new_stat/out/a_fm_A_A12_noC20_unitary_iter0.txt");
    myres->write_jack_in_file(fit_afpi_max_twist_A12_noC20_corr_unitary_iter0.P[1], "../../g-2_new_stat/out/a_fm_B_A12_noC20_unitary_iter0.txt");
    myres->write_jack_in_file(fit_afpi_max_twist_A12_noC20_corr_unitary_iter0.P[2], "../../g-2_new_stat/out/a_fm_C_A12_noC20_unitary_iter0.txt");
    myres->write_jack_in_file(fit_afpi_max_twist_A12_noC20_corr_unitary_iter0.P[3], "../../g-2_new_stat/out/a_fm_D_A12_noC20_unitary_iter0.txt");
    myres->write_jack_in_file(fit_afpi_max_twist_A12_noC20_corr_unitary_iter0.P[4], "../../g-2_new_stat/out/a_fm_E_A12_noC20_unitary_iter0.txt");
    //////////////////////////////////////////////////////////////
    // A12 no C20 + strange misstuning unitary only  (main_analysis)
    //////////////////////////////////////////////////////////////

    fit_info.restore_default();
    // fit_info.N = 8;
    fit_info.Nvar = 10;
    fit_info.Npar = 7;
    fit_info.Njack = Njack;
    fit_info.Nxen = { {A53, A40, A30, A12}, {B25_48, B14_64, B72_64, B72_96},
                        { C06, C112}, {D54},
                        {E112} };
    fit_info.init_N_etot_form_Nxen();

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
    count = 0;
    for (int n = 0;n < fit_info.Nxen.size();n++) {
        for (int e : fit_info.Nxen[n]) {
            for (int j = 0;j < Njack;j++) {
                double my_mu, my_M, my_fpi;

                my_mu = jackextra.en[e].jack[165][j];
                my_M = jackextra.en[e].jack[id_Mpi_cor][j];
                my_fpi = jackextra.en[e].jack[id_fpi_cor][j];


                fit_info.x[0][count][j] = my_mu; // 
                fit_info.x[1][count][j] = my_M;  // 
                fit_info.x[2][count][j] = my_fpi;  //
                fit_info.x[3][count][j] = jackextra.en[e].header.L;

                double xi = my_M / (4 * M_PI * my_fpi);
                xi *= xi;
                double delta_FVE = FVE_GL_Mpi(jackextra.en[e].header.L, xi, my_fpi);
                xi *= (1 + delta_FVE) * (1 + delta_FVE) / (1 - 0.25 * delta_FVE) * (1 - 0.25 * delta_FVE);
                fit_info.x[4][count][j] = xi;

                fit_info.x[5][count][j] = jack_Mpi_phys_MeV[j] / hbarc;
                fit_info.x[6][count][j] = jack_fpi_phys_MeV[j] / hbarc;


                fit_info.x[7][count][j] = jack_Mpi_phys_MeV[j] / (4 * M_PI * jack_fpi_phys_MeV[j]);
                fit_info.x[7][count][j] *= fit_info.x[7][count][j];

                fit_info.x[8][count][j] = vev_mpcac[e][j];// mpcac/mu
                fit_info.x[9][count][j] = jackextra.en[e].jack[23][j];// Z_A

                double mpcac_mu = fit_info.x[8][count][j];
                double ZA = fit_info.x[9][count][j];
                double mr = 0;
                mr = ZA * mpcac_mu;
                double cl = sqrt(1 + mr * mr);
                xi /= cl * cl * cl;
                fit_info.x[4][count][j] = xi;
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
    fit_info.corr_id = { 163, 163, 163, 163, 163 }; //  fpi(mu1), fpi(mu2)
    fit_info.function = rhs_afpi;
    fit_info.linear_fit = false;
    // fit_info.covariancey = true;
    // fit_info.acc = 1e-10;
    // fit_info.h = 1e-7;
    // fit_info.maxiter = 500;
    // fit_info.NM = true;
    // fit_info.chi2_gap_jackboot = 0.01;
    // fit_info.guess_per_jack = 5;
    // fit_info.repeat_start = 100;
    // // fit_info.verbosity = 1;
    // fit_info.compute_cov_fit(argv, jackextra, lhs_afpi_max_twist);
    // int ide = 0, ide1 = 0;
    // for (int n = 0;n < fit_info.Nxen.size();n++) {
    //     for (int e : fit_info.Nxen[n]) {
    //         ide1 = 0;
    //         for (int n1 = 0;n1 < fit_info.Nxen.size();n1++) {
    //             for (int e1 : fit_info.Nxen[n1]) {
    //                 if (e != e1)   fit_info.cov[ide][ide1] = 0;
    //                 printf("%-12.5g ", fit_info.cov[ide][ide1] / sqrt(fit_info.cov[ide][ide] * fit_info.cov[ide1][ide1]));
    //                 ide1++;
    //             }
    //         }
    //         printf("\n");
    //         ide++;
    //     }
    // }
    // fit_info.compute_cov1_fit();
    fit_info.guess = { 0.0908026, 0.07951 , 0.06816, 0.05688, 0.04891 ,-8.75 ,-500 };
    mysprintf(namefit, NAMESIZE, "afpi_max_twist_A12_noC20_cor_cov_unitary");
    fit_result fit_afpi_max_twist_A12_noC20_cor_unitary = fit_all_data(argv, jackextra, lhs_afpi_max_twist, fit_info, namefit);
    fit_info.band_range = { 0.005,0.035 };
    // std::vector<double> xcont = { 0, 0 /*Delta*/, 0, 0,/*l, a,m*/ fit_info.x[4][0][Njack - 1],
    //      fit_info.x[5][0][Njack - 1] , fit_info.x[6][0][Njack - 1], fit_info.x[7][0][Njack - 1] };


    //    Mpi:   the index of the parameter do not match!   P[i]*(M_pi- M_pi_phys ) 
    // std::vector<double> xcont = {};
    print_fit_band(argv, jackextra, fit_info, fit_info, namefit, "xi", fit_afpi_max_twist_A12_noC20_cor_unitary, fit_afpi_max_twist_A12_noC20_cor_unitary, 4, 0, 0.0005, xcont);

    myres->write_jack_in_file(fit_afpi_max_twist_A12_noC20_cor_unitary.P[0], "../../g-2_new_stat/out/a_fm_A_A12_noC20_cor_unitary.txt");
    myres->write_jack_in_file(fit_afpi_max_twist_A12_noC20_cor_unitary.P[1], "../../g-2_new_stat/out/a_fm_B_A12_noC20_cor_unitary.txt");
    myres->write_jack_in_file(fit_afpi_max_twist_A12_noC20_cor_unitary.P[2], "../../g-2_new_stat/out/a_fm_C_A12_noC20_cor_unitary.txt");
    myres->write_jack_in_file(fit_afpi_max_twist_A12_noC20_cor_unitary.P[3], "../../g-2_new_stat/out/a_fm_D_A12_noC20_cor_unitary.txt");
    myres->write_jack_in_file(fit_afpi_max_twist_A12_noC20_cor_unitary.P[4], "../../g-2_new_stat/out/a_fm_E_A12_noC20_cor_unitary.txt");


    //////////////////////////////////////////////////////////////
    // Mpi
    //////////////////////////////////////////////////////////////



    //////////////////////////////////////////////////////////////
    // Mpi with A12,  noC20,  unitary, iter0
    //////////////////////////////////////////////////////////////

    {
        fit_info.restore_default();
        fit_info.N = 5;
        fit_info.Nvar = 11;
        fit_info.Npar = 7;
        fit_info.Njack = Njack;
        fit_info.Nxen = { {A53, A40, A30, A12}, {B25_48, B14_64, B72_64, B72_96},
                            { C06, C112}, {D54},
                            {E112} };
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

                    if (e == A53 || e == A40 || e == A30 || e == A12)
                        fit_info.x[8][count][j] = fit_afpi_max_twist_A12_noC20_corr_unitary_iter0.P[0][j];
                    else if (e == B25_48 || e == B14_64 || e == B72_64 || e == B72_96)
                        fit_info.x[8][count][j] = fit_afpi_max_twist_A12_noC20_corr_unitary_iter0.P[1][j];
                    else if (e == C20 || e == C06 || e == C112)
                        fit_info.x[8][count][j] = fit_afpi_max_twist_A12_noC20_corr_unitary_iter0.P[2][j];
                    else if (e == D54)
                        fit_info.x[8][count][j] = fit_afpi_max_twist_A12_noC20_corr_unitary_iter0.P[3][j];
                    else if (e == E112)
                        fit_info.x[8][count][j] = fit_afpi_max_twist_A12_noC20_corr_unitary_iter0.P[4][j];
                    else { printf("error missing ensemble\n"); exit(1); }

                    fit_info.x[9][count][j] = 0;
                    fit_info.x[9][count][j] = vev_mpcac[e][j];// mpcac/mu
                    fit_info.x[10][count][j] = jackall.en[e].jack[23][j];// Z_A

                }
                count++;
            }
        }



        // fit_info.corr_id = { 1, 123, 163, 164 }; // Mpi(mu1), Mpi(mu2), fpi(mu1), fpi(mu2)
        fit_info.corr_id = { 1, 1, 1, 1, 1, 123, 123, 123 }; //  Mpi(mu1), Mpi(mu2)
        fit_info.function = rhs_aMpi2_over_afpi2_with_A;
        fit_info.linear_fit = false;
        fit_info.covariancey = false;
        fit_info.acc = 1e-12;
        fit_info.h = 1e-7;
        fit_info.maxiter = 500;
        fit_info.NM = true;
        fit_info.chi2_gap_jackboot = 0.01;
        fit_info.guess_per_jack = 5;
        fit_info.repeat_start = 100;
        // fit_info.verbosity = 1;
        fit_info.compute_cov_fit(argv, jackall, lhs_Mpi2_over_afpi2_max_twist);
        int ide = 0, ide1 = 0;
        for (int n = 0;n < fit_info.Nxen.size();n++) {
            for (int e : fit_info.Nxen[n]) {
                ide1 = 0;
                for (int n1 = 0;n1 < fit_info.Nxen.size();n1++) {
                    for (int e1 : fit_info.Nxen[n1]) {
                        if (e != e1)   fit_info.cov[ide][ide1] = 0;
                        printf("%-12.5g ", fit_info.cov[ide][ide1] / sqrt(fit_info.cov[ide][ide] * fit_info.cov[ide1][ide1]));
                        ide1++;
                    }
                }
                printf("\n");
                ide++;
            }
        }
        fit_info.compute_cov1_fit();
        // fit_info.guess = { 0.0072, 0.0060 , 0.0050, 0.0012, -8.75 , -1 };
        // fit_info.guess = { 866, 990 , 1169, 1346, -5.46, 0  };
        fit_info.guess = { 683, 804, 948 , 1164, 1363, 6.0,20 };
        // fit_info.guess = { 1, 1 , 1, 1  };
        // fit_info.guess = { 0.05 };
        mysprintf(namefit, NAMESIZE, "aMpi2_over_afpi2_A12_noC20_cov_unitary_iter0");
        fit_result fit_aMpi2_over_afpi2 = fit_all_data(argv, jackall, lhs_Mpi2_over_afpi2_max_twist, fit_info, namefit);
        fit_info.band_range = { 0.00001,0.0055 };

        print_fit_band(argv, jackall, fit_info, fit_info, namefit, "amu", fit_aMpi2_over_afpi2, fit_aMpi2_over_afpi2, 0 /*amu */, 0, 0.0001, xcont);

        compute_amul_print_res(argv, namefit, fit_info, Njack, fit_aMpi2_over_afpi2, Mpi2_fpi2_phys, { "A","B","C","D","E" });

    }
    //////////////////////////////////////////////////////////////
    // Mpi with A12, max twist corrections, noC20, stange misstuning, unitary (main analysis)
    //////////////////////////////////////////////////////////////

    {
        fit_info.restore_default();
        fit_info.N = 5;
        fit_info.Nvar = 11;
        fit_info.Npar = 7;
        fit_info.Njack = Njack;
        fit_info.Nxen = { {A53, A40, A30, A12}, {B25_48, B14_64, B72_64, B72_96},
                            { C06, C112}, {D54},
                            {E112} };
        fit_info.init_N_etot_form_Nxen();

        fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
        count = 0;
        for (int n = 0;n < fit_info.Nxen.size();n++) {
            for (int e : fit_info.Nxen[n]) {
                for (int j = 0;j < Njack;j++) {
                    double my_mu, my_M, my_fpi;
                    if (n < 5) {
                        my_mu = jackextra.en[e].jack[165][j];
                        my_M = jackextra.en[e].jack[id_Mpi_cor][j];
                        my_fpi = jackextra.en[e].jack[id_fpi_cor][j];
                    }
                    else if (n >= 5) {
                        my_mu = jackextra.en[e].jack[166][j];
                        my_M = jackextra.en[e].jack[id_Mpi_cor_mu1][j];
                        my_fpi = jackextra.en[e].jack[id_fpi_cor_mu1][j];
                    }
                    fit_info.x[0][count][j] = my_mu; // 
                    fit_info.x[1][count][j] = my_M;  // 
                    fit_info.x[2][count][j] = my_fpi;  //
                    fit_info.x[3][count][j] = jackextra.en[e].header.L;

                    double xi = my_M / (4 * M_PI * my_fpi);
                    xi *= xi;
                    double delta_FVE = FVE_GL_Mpi(jackextra.en[e].header.L, xi, my_fpi);
                    xi *= (1 + delta_FVE) * (1 + delta_FVE) / (1 - 0.25 * delta_FVE) * (1 - 0.25 * delta_FVE);
                    fit_info.x[4][count][j] = xi;

                    fit_info.x[5][count][j] = jack_Mpi_phys_MeV[j] / hbarc;
                    fit_info.x[6][count][j] = jack_fpi_phys_MeV[j] / hbarc;


                    fit_info.x[7][count][j] = jack_Mpi_phys_MeV[j] / (4 * M_PI * jack_fpi_phys_MeV[j]);
                    fit_info.x[7][count][j] *= fit_info.x[7][count][j];

                    if (e == A53 || e == A40 || e == A30 || e == A12)
                        fit_info.x[8][count][j] = fit_afpi_max_twist_A12_noC20_cor_unitary.P[0][j];
                    else if (e == B25_48 || e == B14_64 || e == B72_64 || e == B72_96)
                        fit_info.x[8][count][j] = fit_afpi_max_twist_A12_noC20_cor_unitary.P[1][j];
                    else if (e == C20 || e == C06 || e == C112)
                        fit_info.x[8][count][j] = fit_afpi_max_twist_A12_noC20_cor_unitary.P[2][j];
                    else if (e == D54)
                        fit_info.x[8][count][j] = fit_afpi_max_twist_A12_noC20_cor_unitary.P[3][j];
                    else if (e == E112)
                        fit_info.x[8][count][j] = fit_afpi_max_twist_A12_noC20_cor_unitary.P[4][j];
                    else { printf("error missing ensemble\n"); exit(1); }

                    fit_info.x[9][count][j] = vev_mpcac[e][j];// mpcac/mu
                    fit_info.x[10][count][j] = jackextra.en[e].jack[23][j];// Z_A

                }
                count++;
            }
        }



        // fit_info.corr_id = { 1, 123, 163, 164 }; // Mpi(mu1), Mpi(mu2), fpi(mu1), fpi(mu2)
        fit_info.corr_id = { 1, 1, 1, 1, 1, 123, 123, 123 }; //  Mpi(mu1), Mpi(mu2)
        fit_info.function = rhs_aMpi2_over_afpi2_with_A;
        fit_info.linear_fit = false;
        fit_info.covariancey = false;
        fit_info.acc = 1e-12;
        fit_info.h = 1e-7;
        fit_info.maxiter = 500;
        fit_info.NM = true;
        fit_info.chi2_gap_jackboot = 0.01;
        fit_info.guess_per_jack = 5;
        fit_info.repeat_start = 100;
        // fit_info.verbosity = 1;
        fit_info.compute_cov_fit(argv, jackextra, lhs_Mpi2_over_afpi2_max_twist);
        int ide = 0, ide1 = 0;
        for (int n = 0;n < fit_info.Nxen.size();n++) {
            for (int e : fit_info.Nxen[n]) {
                ide1 = 0;
                for (int n1 = 0;n1 < fit_info.Nxen.size();n1++) {
                    for (int e1 : fit_info.Nxen[n1]) {
                        if (e != e1)   fit_info.cov[ide][ide1] = 0;
                        printf("%-12.5g ", fit_info.cov[ide][ide1] / sqrt(fit_info.cov[ide][ide] * fit_info.cov[ide1][ide1]));
                        ide1++;
                    }
                }
                printf("\n");
                ide++;
            }
        }
        fit_info.compute_cov1_fit();
        // fit_info.guess = { 0.0072, 0.0060 , 0.0050, 0.0012, -8.75 , -1 };
        // fit_info.guess = { 866, 990 , 1169, 1346, -5.46, 0  };
        fit_info.guess = { 683, 804, 948 , 1164, 1363, 6.0,20 };
        // fit_info.guess = { 1, 1 , 1, 1  };
        // fit_info.guess = { 0.05 };
        mysprintf(namefit, NAMESIZE, "aMpi2_over_afpi2_A12_noC20_cor_cov_unitary");
        fit_result fit_aMpi2_over_afpi2 = fit_all_data(argv, jackextra, lhs_Mpi2_over_afpi2_max_twist, fit_info, namefit);
        fit_info.band_range = { 0.00001,0.0055 };

        print_fit_band(argv, jackextra, fit_info, fit_info, namefit, "amu", fit_aMpi2_over_afpi2, fit_aMpi2_over_afpi2, 0 /*amu */, 0, 0.0001, xcont);

        compute_amul_print_res(argv, namefit, fit_info, Njack, fit_aMpi2_over_afpi2, Mpi2_fpi2_phys, { "A","B","C","D","E" });

    }








}