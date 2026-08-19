#define CONTROL

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "global.hpp"
#include "mutils.hpp"
#include "read.hpp"
#include "resampling.hpp"
// #include "m_eff.hpp"
// #include "gnuplot.hpp"
#include "eigensystem.hpp"
#include "functions_gm2_analysis.hpp"
#include "linear_fit.hpp"
#include "mutils.hpp"
#include "various_fits.hpp"
// #include "correlators_analysis.hpp"
// #include "eigensystem.hpp"
#include "fit_all.hpp"
#include "global.hpp"
#include "non_linear_fit.hpp"
#include "resampling_new.hpp"
#include "tower.hpp"

#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

double rhs_linear(int n, int Nvar, double *x, int Npar, double *P) {
  double r;
  double fpi = x[0];
  r = P[0] + fpi * P[1];
  return r;
}

double lhs_fun(int n, int e, int j, data_all gjack, struct fit_type fit_info) {
  //   printf("id = %d e=%d n=%d j=%d\n", fit_info.corr_id[0], e, n, j);
  //   printf("%g\n",gjack.en[e].jack[fit_info.corr_id[0]][j]);
  return gjack.en[e].jack[fit_info.corr_id[0]][j];
}

int main(int argc, char **argv) {
  error(argc != 1, 1, "main ", "usage:./program ");
  char namefile[NAMESIZE];

  char **options = (char **)malloc(sizeof(char *) * 4);
  options[0] = (char *)malloc(sizeof(char) * NAMESIZE);
  options[1] = (char *)malloc(sizeof(char) * NAMESIZE);
  options[2] = (char *)malloc(sizeof(char) * NAMESIZE);
  options[3] = (char *)malloc(sizeof(char) * NAMESIZE);

  int Njack = 51;
  myres = new resampling_jack(Njack - 1);
  std::vector<int> iWs = {0, 1, 2, 3, 5};
  //////////////////////////////////////////////////////////////
  //  jackall
  //////////////////////////////////////////////////////////////
  data_all jackall;
  jackall.resampling = "jack";
  mysprintf(options[1], NAMESIZE, "%s", jackall.resampling.c_str());
  mysprintf(options[3], NAMESIZE, "interpolation_fpi");

  // jackall->en = (data_single*)malloc(sizeof(data_single) * files.size());
  jackall.ens = 2;
  jackall.en = new data_single[jackall.ens];
  int count = 0;
  int Nobs = std::ranges::max(iWs)+1;
  std::vector<double> fpi = {130.5, 131.1};
  for (int e = 0; e < jackall.ens; e++) {

    data_single dj;
    // dj.header = read_header(f);
    dj.Nobs = Nobs;
    dj.Njack = Njack;
    dj.jack = malloc_2<double>(dj.Nobs, dj.Njack);

    //
    size_t i = 0;
    // for (int obs = 0; obs < dj.Nobs; obs++) {
    for (int obs : iWs) {
      mysprintf(namefile, NAMESIZE,
                "../../g-2_new_stat/%.1f/fit_all_charm//ave_BAIC_%d_%d.jack",
                fpi[e], obs, Njack);
      printf("reading %s\n", namefile);
      myres->read_jack_from_file(dj.jack[obs], namefile);
    }
    dj.resampling = jackall.resampling;

    jackall.en[e] = dj;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // fits
  /////////////////////////////////////////////////////////////////////////////////////////////////
  fit_type fit_info;
  double *fpi_wp25 = myres->create_fake_exact(fpi_wp25_val, fpi_wp25_err, 1);
  double *interp = myres->create_zero();

  // for (int obs = 0; obs < jackall.en[0].Nobs; obs++) {
  for (int obs : iWs) {
    printf("obs=%d\n", obs);
    fit_info.corr_id = std::vector<int>{obs};
    fit_info.Nxen = {{0, 1}};
    fit_info.init_N_etot_form_Nxen();
    fit_info.Nvar = 1;
    fit_info.Npar = 2;
    fit_info.function = rhs_linear;
    fit_info.Njack = Njack;

    fit_info.x = double_malloc_3(fit_info.Nvar, fit_info.entot, fit_info.Njack);
    count = 0;
    for (int n = 0; n < fit_info.N; n++) {
      for (int e : fit_info.Nxen[n]) {
        for (int j = 0; j < Njack; j++) {
          fit_info.x[0][count][j] = fpi[e];
        }
        count++;
      }
    }
    fit_info.linear_fit = true;
    fit_info.verbosity = 0;
    fit_info.covariancey = true;
    fit_info.compute_cov_fit(options, jackall, lhs_fun);
    fit_info.compute_cov1_fit();
    std::string namefit = std::string("fit_fpi_obs_") + std::to_string(obs);
    fit_result fit_out =
        fit_all_data(options, jackall, lhs_fun, fit_info, namefit.c_str());
    fit_info.band_range = {129.5, 132};
    print_fit_band(options, jackall, fit_info, fit_info, namefit.c_str(), "fpi",
                   fit_out, fit_out, 0, fit_info.Nxen[0][0], 0.0002);

    for (size_t j = 0; j < Njack; j++) {
      interp[j] = fit_out.P[0][j] + fpi_wp25[j] * fit_out.P[1][j];
      /* code */
    }
    printf("obs = %-4d  value = %-15.12g  %-15.12g\n", obs, myres->mean(interp),
           myres->comp_error(interp));
    fit_info.restore_default();
  }
}