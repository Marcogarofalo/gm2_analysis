#ifndef functions_gm2_analysis_H
#define functions_gm2_analysis_H
#include "non_linear_fit.hpp"

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

#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include "tower.hpp"
#include "fit_all.hpp"
extern "C" {
#include "dzeta_function.h"
#include "qzeta_function.h"
}
#include <algorithm>
#include "non_linear_fit.hpp"
static int once = 0;

using namespace std::complex_literals;
constexpr double hbarc = 197.326963;
constexpr double muon_mass_MeV = 105.6583755;
constexpr double alpha_em = 0.0072973525693;
constexpr double muon_mass_fm = muon_mass_MeV * 197.326963;
constexpr double Mpi_MeV = 135;
constexpr double Mpi_MeV_err = 0.2;

enum enum_ensembles {
    B72_64,
    B72_96,
    C06,
    D54,
    A53,
    A40,
    A30
};

double integrand_K(double x, void* params);
double kernel_K(double z, double epsrel = 1e-7);
double gm2_step_function(double t_d, double t1_d);
double integrand_K_W(double x, void* params);
double kernel_K_W(double z, double epsrel = 1e-4);
double integrate_simpson38(int lower, int upper, double* f);
double integrate_reinman(int lower, int upper, double* f);
double* compute_amu_full(double**** in, int id, int Njack, double* Z, double* a, double q2, double (*int_scheme)(int, int, double*), FILE* outfile, const char* description, const char* resampling, int isub = -1);
double* compute_amu_sd(double**** in, int id, int Njack, double* Z, double* a, double q2, double (*int_scheme)(int, int, double*), FILE* outfile, const char* description, const char* resampling, int isub = -1);
double* compute_amu_W(double**** in, int id, int Njack, double* Z, double* a, double q2, double (*int_scheme)(int, int, double*), FILE* outfile, const char* description, const char* resampling);
template<int idn, int idd>
double ZAl_lhs(int j, double**** in, int t, struct fit_type fit_info);
template<int idn, int idd>
double ZVl_lhs(int j, double**** in, int t, struct fit_type fit_info);
template<int id>
double GPS_OS_lhs(int j, double**** in, int t, struct fit_type fit_info);
template<int id>
double GPS_lhs(int j, double**** in, int t, struct fit_type fit_info);
template<int id>
double lhs_ct(int j, double**** in, int t, struct fit_type fit_info);
double* interpol_Z(int Nmus, int Njack, double** Meta, double** Z, double* aMetas_exp,
    FILE* outfile, const char* description, const char* resampling);
double compute_cotd(int Nvar, double* x);
double compute_delta(double omega, void* params);
double compute_deriv_delta(int Nvar, double* x);
double w_js(int j, int s, double q);
double Z00(double q);
double omega_QC2(int n, int Nvar, double* x, int Npar, double* P);
double compute_phi(double q, void* params);
double compute_deriv_phi(double q);
inline double compute_FGS2(int Nvar, double* x) {
    double omega = x[0];
    double mass = x[1];
    double L = x[2];
    double Mrho = x[3];
    double grhopipi = x[4];
    double a = x[5];
    double k = sqrt(omega * omega / 4. - mass * mass);
    // double q = k * (L * (a / 197.326963) / (2. * pi_greco));

    double ho = (grhopipi * grhopipi * k * k * k * 2) * log((omega + 2 * k) / (2 * mass)) / (6.0 * M_PI * M_PI * omega);
    double kM = sqrt(Mrho * Mrho / 4. - mass * mass);
    double hM = (grhopipi * grhopipi * kM * kM * kM * 2) * log((Mrho + 2 * kM) / (2 * mass)) / (6.0 * M_PI * M_PI * Mrho);
    double h1M = (grhopipi * grhopipi * kM * kM) * (1 + (1 + 2 * mass * mass / (Mrho * Mrho)) * (Mrho / kM) * log((Mrho + 2 * kM) / (2 * mass))) / (6.0 * M_PI * M_PI * Mrho);
    double gamma = (grhopipi * grhopipi * k * k * k) / (6.0 * M_PI * omega * omega);

    double Apipi0 = hM - Mrho * h1M / 2. + pow(grhopipi * mass, 2) / (6 * M_PI * M_PI);
    std::complex<double> Apipi = hM + (omega * omega - Mrho * Mrho) * h1M / (2 * Mrho) - ho + 1i * omega * gamma;
    std::complex<double> FGS = (Mrho * Mrho - Apipi0) / (Mrho * Mrho - omega * omega - Apipi);
    double FGS2 = real(FGS * conj(FGS));
    if (once == -1) {
        printf("omega=%.15g    Mpi=%.15g   Mrho=%.15g  g=%.15g\n", omega, mass, Mrho, grhopipi);
    }
    return FGS2;
}
double matrix_element_nuA2(int Nvar, double* x);
double integrand_V_infL(double omega, void* params);
double* compute_V_infL(int Nvar, double* x, int T);
double compute_V_infL_t(int Nvar, double* x, double t);
double* compute_VL(int n, double** omega, double** nuA2, double* a, int T, int j);
double compute_VL_t(int n, double* omega, double* nuA2, double t);
double** compute_DVt(int L, int Njack, double* Mpi, double* Mrho, double* a, double* grhopipi, FILE* outfile, const char* description, const char* resampling);
double integrand_DV(double t, void* params);
double* compute_DVt_and_integrate(int L, int Njack, double* Mpi, double* Mrho, double* a, double* grhopipi, FILE* outfile, const char* description, const char* resampling);
double rhs_amu_log_a4(int n, int Nvar, double* x, int Npar, double* P);
double linear_fit_mu_correction(int n, int Nvar, double* x, int Npar, double* P);
double exp_MpiL(int n, int Nvar, double* x, int Npar, double* P);
double const_A(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_RF(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_separate(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_a4(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_a4_charm(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_diff_ratio_charm(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_diff_ratio(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_cut(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_cut_charm(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_pade(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_FVE_RF(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_diff_RF(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_diff(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_ratio(int n, int Nvar, double* x, int Npar, double* P) ;
double rhs_amu_common(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a2_FVE(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a2_FVE_log_eq(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a2_FVE_log_op(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a2_FVE_log_eq_op(int n, int Nvar, double* x, int Npar, double* P);

double rhs_amu_common_a2_FVE_log_a4(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a2_FVE_log_eq_a4_eq(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a2_FVE_log_op_a4_eq(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a2_FVE_log_eq_op_a4_eq(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a2_FVE_a4_eq(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a2_FVE_a4_op(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a4(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a4_n0(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_a4_n1(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_log_a4_n0(int n, int Nvar, double* x, int Npar, double* P);
double rhs_amu_common_log_a4_n1(int n, int Nvar, double* x, int Npar, double* P);

template<int ieq, int iop>
double lhs_amu_common(int n, int e, int j, data_all gjack, struct fit_type fit_info);

double lhs_amu_common_GS(int n, int e, int j, data_all gjack, struct fit_type fit_info);



double lhs_amu_common_GS_diff(int n, int e, int j, data_all gjack, struct fit_type fit_info);

double lhs_amu_diff(int n, int e, int j, data_all gjack, struct fit_type fit_info);

double lhs_amu_ratio(int n, int e, int j, data_all gjack, struct fit_type fit_info);


double lhs_amu(int n, int e, int j, data_all gjack, struct fit_type fit_info);
double lhs_Acharm(int n, int e, int j, data_all gjack, struct fit_type fit_info);


double lhs_amu_separate(int n, int e, int j, data_all gjack, struct fit_type fit_info);


double lhs_amu_diff_ratio(int n, int e, int j, data_all gjack, struct fit_type fit_info);
double lhs_mu_corrections(int n, int e, int j, data_all gjack, struct fit_type fit_info);

double lhs_SD_mu_corrections(int n, int e, int j, data_all gjack, struct fit_type fit_info);



template<int ieq, int iop>
double lhs_amu_common_FVE(int n, int e, int j, data_all gjack, struct fit_type fit_info);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// print band
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void print_fit_band_amu_W_l(char** argv, data_all gjack, struct fit_type fit_info,
    struct fit_type fit_info_m0, const char* label, const char* dir_name,
    struct fit_result fit_out, struct fit_result fit_out_m0, int var, int en, double h, std::vector<double> xval,
    double lhs_fun(int, int, int, data_all, struct fit_type));

double rhs_2exp(int n, int Nvar, double* x, int Npar, double* P);

double rhs_1exp(int n, int Nvar, double* x, int Npar, double* P) ;
double rhs_poly(int n, int Nvar, double* x, int Npar, double* P);


double** corr_plus_dm(int j, double**** in, int t, struct fit_type fit_info);

double** corr_plus_dm_correlated(int j, double**** in, int t, struct fit_type fit_info);

double** mu_sea_correction(int j, double**** in, int t, struct fit_type fit_info);


void compute_syst_eq28(data_all in, const char* outpath, const char* filename) ;

#endif
