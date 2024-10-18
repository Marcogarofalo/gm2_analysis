#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#include "global.hpp"
#include "resampling.hpp"
#include "resampling_new.hpp"
#include "read.hpp"
#include "m_eff.hpp"
#include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "linear_fit.hpp"
#include "various_fits.hpp"
#include "mutils.hpp"
#include "functions_gm2_analysis.hpp"
#include "correlators_analysis.hpp"
#include "non_linear_fit.hpp"
#include "tower.hpp"
#include "fit_all.hpp"
#include "functions_amu.hpp"
#include "gamma_analysis.hpp"
// #include "lhs_functions.hpp"

#include <string>
#include <cstring> 
#include <string>
#include <fstream>
#include <memory>

#include <gsl/gsl_integration.h>

// constexpr double muon_mass_MeV = 105.6583755;
// constexpr double alpha_em = 0.0072973525693;

struct params_integrand_amu {
    int nE;
    double* R;
    double* E; // GeV
    double* integrand_R;
    double t_fm;
};


// Energy in Gev
// t in fm
// integral will be given in [fm][GeV]^3, if you multiply the result by 1000/hbarc we get 
// [fm/hbarc] [GeV*1000] [GeV]^2 = [MeV]^-1 [MeV] [GeV]^2 =[GeV]^2
double integrand_amu(double x, void* params) {

    double r;
    double t = x;
    constexpr double d = 0.15;
    constexpr double t1_d = 0.4 / d;
    params_integrand_amu* p = (params_integrand_amu*)params;
    for (int nl = 0;nl < p->nE;nl++) {
        p->integrand_R[nl] = p->E[nl] * p->E[nl] * p->R[nl] * exp(-(t / hbarc) * p->E[nl] * 1000);
    }

    double Vt = integrate_simpson38(0, p->nE-1, p->integrand_R);
    Vt *= (p->E[1] - p->E[0]) / (12 * M_PI * M_PI);
    double z = muon_mass_MeV * (t / hbarc);
    double K;

    K = z * z * kernel_K(z);
    if ((p->E[1] - p->E[0]) < 0) printf("error dE  %g   %g\n", p->E[1], p->E[0]);
    if (K < 0) printf("error K\n");

    double th = gm2_step_function(t / 0.15, t1_d);
    if ((1 - th) < 0)printf("error th\n");

    r = K * Vt * (1 - th);

    return r;
}


int main(int argc, char** argv) {

    // FILE* frhad = open_file("r_had", "r+");
    std::ifstream infile;

    int num = 0; // num must start at 0
    infile.open("../r_had_formatted");// file containing numbers in 3 columns 
    std::string line;
    int nlines = 0;
    while (std::getline(infile, line)) {
        nlines++;
    }
    std::cout << "nlines: " << nlines << "\n";
    int nenergy = nlines - 1;
    std::vector<double> En(nenergy);
    std::vector<double> Ru(nenergy);
    std::vector<double> Rd(nenergy);
    std::vector<double> Rs(nenergy);
    std::vector<double> Rc(nenergy);
    std::vector<double> Rb(nenergy);
    std::vector<double> Rtot(nenergy);
    std::vector<double> integrand_Rs(nenergy);
    std::vector<double> integrand_Rc(nenergy);

    int nt = 200;

    infile.clear();
    infile.seekg(0); // rewind file to line 1
    std::getline(infile, line);
    int nl = 0;
    while (std::getline(infile, line)) {
        // nlines++;
        std::stringstream ss(line);
        ss >> En[nl] >> Ru[nl] >> Rd[nl] >> Rs[nl] >> Rc[nl] >> Rb[nl] >> Rtot[nl];
        nl++;
    }


    params_integrand_amu par;
    par.nE = nenergy;
    par.R = Rs.data();
    par.E = En.data();
    par.integrand_R = (double*)malloc(sizeof(double) * nenergy);

    int Maxiter = 1e+6;
    double epsrel = 1e-4;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(Maxiter);
    double result, error;
    gsl_function F;
    F.function = &integrand_amu;
    F.params = &par;

    static constexpr double muon_mass_GeV = muon_mass_MeV / 1000.0;
    double tmin_MeV = par.t_fm / hbarc;
    double tmin_GeV = tmin_MeV * 1000;
    double dt = 0.07951; // a(E)= 0.04891  fm


    std::string name_out_s("/home/garofalo/analysis/gm2_analysis/gm2_book/amu_SD_s_pert.txt");
    printf("writing in %s\n", name_out_s.c_str());
    FILE* out_s = open_file(name_out_s.c_str(), "w+");
    fprintf(out_s, "t_fm     amu_(s)_pert\n");

    par.t_fm = dt; //    , +0.07951, 0.07951*3/2.0,  0.07951*2
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 2.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 3.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 4.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    fclose(out_s);
    //////////////////////////////////////////////////////////////
    // charm
    //////////////////////////////////////////////////////////////

    par.R = Rc.data();

    name_out_s = "/home/garofalo/analysis/gm2_analysis/gm2_book/amu_SD_c_pert.txt";
    printf("writing in %s\n", name_out_s.c_str());
    out_s = open_file(name_out_s.c_str(), "w+");
    fprintf(out_s, "t_fm     amu_(s)_pert\n");

    par.t_fm = dt; //    , +0.07951, 0.07951*3/2.0,  0.07951*2
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 2.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 3.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 4.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    fclose(out_s);

    infile.close();
    //////////////////////////////////////////////////////////////
    // bar MS
    //////////////////////////////////////////////////////////////
    infile.open("../r_had_barMS");// file containing numbers in 3 columns 
    nlines = 0;
    while (std::getline(infile, line)) {
        nlines++;
    }
    std::cout << "nlines: " << nlines << "\n";
    nenergy = nlines - 1;
    En.resize(nenergy);
    Ru.resize(nenergy);
    Rd.resize(nenergy);
    Rs.resize(nenergy);
    Rc.resize(nenergy);
    Rb.resize(nenergy);
    Rtot.resize(nenergy);
    integrand_Rs.resize(nenergy);
    integrand_Rc.resize(nenergy);

    infile.clear();
    infile.seekg(0); // rewind file to line 1
    std::getline(infile, line);
    nl = 0;
    while (std::getline(infile, line)) {
        // nlines++;
        std::stringstream ss(line);
        ss >> En[nl] >> Ru[nl] >> Rd[nl] >> Rs[nl] >> Rc[nl] >> Rb[nl] >> Rtot[nl];
        nl++;
    }

    par.nE = nenergy;
    par.R = Rs.data();
    par.E = En.data();
    free(par.integrand_R);
    par.integrand_R = (double*)malloc(sizeof(double) * nenergy);


    name_out_s = "/home/garofalo/analysis/gm2_analysis/gm2_book/amu_SD_s_pert_barMS.txt";
    printf("writing in %s\n", name_out_s.c_str());
    out_s = open_file(name_out_s.c_str(), "w+");
    fprintf(out_s, "t_fm     amu_(s)_pert\n");

    par.t_fm = dt; //    , +0.07951, 0.07951*3/2.0,  0.07951*2
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 2.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 3.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 4.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    fclose(out_s);
    //////////////////////////////////////////////////////////////
    // charm
    //////////////////////////////////////////////////////////////

    par.R = Rc.data();

    name_out_s = "/home/garofalo/analysis/gm2_analysis/gm2_book/amu_SD_c_pert_barMS.txt";
    printf("writing in %s\n", name_out_s.c_str());
    out_s = open_file(name_out_s.c_str(), "w+");
    fprintf(out_s, "t_fm     amu_(s)_pert\n");

    par.t_fm = dt; //    , +0.07951, 0.07951*3/2.0,  0.07951*2
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 2.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 3.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    par.t_fm = dt + dt * 4.0 / 2.0;
    gsl_integration_qags(&F, 0, par.t_fm, 0, epsrel, Maxiter, w, &result, &error);
    result *= (1000 / hbarc) * 4 * alpha_em * alpha_em / (muon_mass_GeV * muon_mass_GeV);
    fprintf(out_s, "%-20.12g  %.12g\n", par.t_fm, result);

    fclose(out_s);

    infile.close();

    gsl_integration_workspace_free(w);



    return result;


}