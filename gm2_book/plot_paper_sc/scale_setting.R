library(Rose)
library(ggplot2)
library(plotly)
library(knitr)
source("/home/garofalo/programs/Rose/R/plot_routines.R")

dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all/"
namefit <- paste0(dir, "afpi_max_twist_A12_noC20_cor_cov_unitary")
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
  data_type = c("A,i=2", "B,i=2", "C,i=2", "D,i=2", "E,i=2", "B,i=2", "C,i=2", "D,i=2"),
  width = 2e-4,
  gg = gg, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, size=1.5, alpha_ribbon = 0.2
)
gg <- gg + geom_vline(xintercept = 6.77683303196e-3)
fig <- myplotly(gg, "", "$\\xi$", "$af_\\pi$", to_print = FALSE, 
                save_pdf="fpi_vs_xi",  xrange = c(0.007,0.034),
                legend_position = c(0.8,0.65))
###################################

gg <- myggplot(repeat_color = 1)
gg <- plot_fit(
  basename = namefit,
  var = "xi",
  id_x = 5,
  data_type = c("A,i=2", "B,i=2", "C,i=2", "D,i=2", "E,i=2"),
  width = 0.3e-4,
  gg = gg, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, size=1.5, alpha_ribbon = 0.2
)
namefit <- paste0(dir, "afpi_max_twist_A12_noC20_cov_unitary_iter0")
gg <- plot_fit(
  basename = namefit,
  var = "xi",
  id_x = 5,
  data_type = c("A,i=1", "B,i=1", "C,i=1", "D,i=1", "E,i=1"),
  width = 0.3e-4,
  gg = gg, 
  labelfit = "",
  nudge = 0.00004,
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, noribbon = TRUE, size=1.5, alpha_ribbon = 0.2
)
# gg <- plot_fit(
#   basename = namefit,
#   var = "xi",
#   id_x = 5,
#   data_type = c("A", "B", "C", "D", "E"),
#   width = 0.5e-4,
#   gg = gg, 
#   labelfit = "",
#   # , single_name_for_fit = "fit"
#   # , nolabel_for_fit = TRUE
#   noline = TRUE, size=1.5, alpha_ribbon = 0.2
# )
xiphys<-6.77683303196e-3
gg <- gg + geom_vline(xintercept = xiphys)

df<-read.table(paste0(dir, "afpi_max_twist_A12_noC20_cor_cov_fit_data.txt"))
df<-df[c(13:15),]
iy<-dim(df)[2]-2
lab<-c("B+direct correction","C+direct correction","D+direct correction")
# the real x should be df[,5]
gg<- gg + geom_point(aes(x=xiphys, y=df[,iy],
                               color=lab, shape=lab, fill=lab))
gg<- gg + geom_errorbar(aes(x=xiphys,
                              ymin=df[,iy]-df[,iy+1],
                              ymax=df[,iy]+df[,iy+1], color=lab, shape=lab, fill=lab),
                        width = 0.5e-4)

fig <- myplotly(gg, "", "$\\xi$", "$af_\\pi$", to_print = FALSE, 
                save_pdf="fpi_vs_xi_zoom", xrange = c(0.006,0.009),
                yrange = c(0.037,0.054),
                legend_position = c(0.8,0.65))

###############################################################################

dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all/"


namefit <- paste0(dir, "aMpi2_over_afpi2_A12_noC20_cor_cov_unitary")



gg <- plot_fit(
  basename = namefit,
  var = "amu",
  id_x = 1,
  data_type = c("A,i=2", "B,i=2", "C,i=2", "D,i=2", "E,i=2"),
  width = 0.7e-4,
  gg = NULL, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, size=1.5, alpha_ribbon = 0.2
)

gg <- gg + geom_hline(yintercept = 1.07015457788)
fig <- myplotly(gg, "", "$a\\mu_\\ell$", "$(M_\\pi/f_\\pi)^2$", to_print = FALSE,
                save_pdf="Mpi_over_fpi" , yrange = c(0,5.5),
                legend_position = c(-0.06,0.99))
#########


gg <- plot_fit(
  basename = namefit,
  var = "amu",
  id_x = 1,
  data_type = c("A,i=2", "B,i=2", "C,i=2", "D,i=2", "E,i=2", "B,i=2", "C,i=2", "D,i=2"),
  width = 0.1e-4,
  gg = NULL, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, size=1.5, alpha_ribbon = 0.2
)

namefit <- paste0(dir, "aMpi2_over_afpi2_A12_noC20_cov_unitary_iter0")
gg <- plot_fit(
  basename = namefit,
  var = "amu",
  id_x = 1,
  data_type = c("A,i=1", "B,i=1", "C,i=1", "D,i=1", "E,i=1"),
  width = 0.1e-4,
  gg = gg, 
  labelfit = "",
  nudge = 0.5e-5,
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, noribbon = TRUE, size=1.5, alpha_ribbon = 0.2
)
df<-read.table(paste0(dir, "aMpi2_over_afpi2_A12_noC20_cor_cov_fit_data.txt"))
df<-df[c(13:15),]
iy<-dim(df)[2]-2
lab<-c("B+direct correction","C+direct correction","D+direct correction")
Mpi_fpi_phys<-1.07015457788
# the real x should be df[,1]
amu_phys<-c(#0.000748072311508,
            0.000666857356526,
            0.000586436819304,
            0.000493352338666
            #0.000430647061354
            )
gg<- gg + geom_point(aes(x=amu_phys, y=Mpi_fpi_phys,
                         color=lab, shape=lab, fill=lab))
gg<- gg + geom_errorbar(aes(x=amu_phys,
                            ymin=Mpi_fpi_phys-df[,iy+1],
                            ymax=Mpi_fpi_phys+df[,iy+1], color=lab, shape=lab, fill=lab),
                        width = 0.1e-4)

gg <- gg + geom_hline(yintercept = 1.07015457788)
fig <- myplotly(gg, "", "$a\\mu_\\ell$", "$(M_\\pi/f_\\pi)^2$", to_print = FALSE,
                save_pdf="Mpi_over_fpi_zoom" , yrange = c(1,1.2),
                xrange = c(0.00035,0.00072),
                legend_position = c(0.06,0.99))
