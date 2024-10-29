library(Rose)
library(ggplot2)
library(plotly)
library(knitr)
source("/home/garofalo/programs/Rose/R/plot_routines.R")

dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all/"
namefit <- paste0(dir, "afpi_max_twist_A12_noC20_cor_cov")
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
cat("\n\n")
dt <- make_table_fit_result(fit)
# dt<- data.frame("par"=fit$P[,1] )
# dt$value<- mapply(mean_print, fit$P[,2], fit$P[,3])
# dt$no_max_twist<- mapply(mean_print, fit_no_max_twist$P[,2], fit_no_max_twist$P[,3])
#print(dt)
print(dt)
gg <- myggplot(repeat_color = 1)
gg <- plot_fit(
  basename = namefit,
  var = "xi",
  id_x = 5,
  data_type = c("A", "B", "C", "D", "E", "B", "C", "D"),
  width = 2e-4,
  gg = gg, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, size=1.5, alpha_ribbon = 0.2
)
gg <- gg + geom_vline(xintercept = 0.00678723)
fig <- myplotly(gg, "", "$\\xi$", "$af_\\pi$", to_print = FALSE, 
                save_pdf="fpi_vs_xi")


###############################################################################

dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all/"


namefit <- paste0(dir, "aMpi2_over_afpi2_A12_noC20_cor_cov")



gg <- plot_fit(
  basename = namefit,
  var = "amu",
  id_x = 1,
  data_type = c("A", "B", "C", "D", "E", "B", "C", "D"),
  width = 0.7e-4,
  gg = NULL, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, size=1.5, alpha_ribbon = 0.2
)

gg <- gg + geom_hline(yintercept = 1.07015457788)
fig <- myplotly(gg, "", "$a\\mu_\\ell$", "$(M_\\pi/f_\\pi)^2$", to_print = FALSE,
                save_pdf="Mpi_over_fpi" , yrange = c(0,6))
