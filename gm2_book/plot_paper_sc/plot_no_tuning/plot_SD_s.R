library(Rose)
library(ggplot2)
library(plotly)
library(knitr)



library("gridExtra")
library("ggpubr")
library(stringr)

source("../functions.R")
source("/home/garofalo/programs/Rose/R/plot_routines.R")
# df<-plot_all_fits("SDetasNoCor","SD","s",mybinwidth=5e-13,
#                   myyrange=c(8.75e-10,9.35e-10),
#                   quark_full_name="strange_no_fine_tuning_latt")
# df<-plot_all_fits("SD","SD","c",mybinwidth=5e-13,
#                   myyrange=c(11.2e-10,12.95e-10),
#                   quark_full_name="charm_no_fine_tuning_latt")
# #plot_all_fits("SD","SD","c",quark_full_name="charm_no_fine_tuning_latt")
# 
# df<-plot_all_fits("WetasNoCor","W","s",mybinwidth=5e-13,
#                   myyrange=c(26.78e-10,27.6e-10),
#                   quark_full_name="strange_no_fine_tuning_latt")
# df<-plot_all_fits("W","W","c",mybinwidth=5e-13,
#                   myyrange=c(2.3e-10,4.5e-10),
#                   quark_full_name="charm_no_fine_tuning_latt")
# 
# df<-plot_all_fits("LDetasNoCor","LD","s",mybinwidth=5e-13,
#                   myyrange=c(16.0e-10,17.8e-10),
#                   quark_full_name="strange_no_fine_tuning_latt")
# df<-plot_all_fits("LD","LD","c",  mybinwidth=1e-14,
#                   myyrange=c(1.25e-12,2.1e-12),
#                   quark_full_name="charm_no_fine_tuning_latt")
# 
# df<-plot_all_fits("SDpWpLDetasNoCor","","s",mybinwidth=5e-13,
#                   myyrange=c(51.5e-10,54.5e-10),
#                   leg_pos= c(0.1, 0.3),
#                   quark_full_name="strange_no_fine_tuning_latt")
# df<-plot_all_fits("SDpWpLD","","c",mybinwidth=5e-13,
#                   myyrange=c(14.1e-10,17.1e-10),
#                   quark_full_name="charm_no_fine_tuning_latt")
gg <- plot_fit(
  basename = "/home/garofalo/analysis/g-2_new_stat/fit_all_strange/amu_SDpWpLDcor_3b_BOS_BTM_a4OS_a4TM",
  var = "afm",
  id_x = 1,
  data_type = c("OS, i=2", "TM, i=2"),
  width = 0.5e-4,
  gg = NULL,
  single_name_for_fit = "",
  noline = TRUE,noribbon = TRUE, size=1.5
)
gg <- plot_fit(
  basename = "/home/garofalo/analysis/g-2_new_stat/fit_all_strange_no_fine_tuning_latt/amu_SDpWpLD_3b_BOS_BTM_a4OS_a4TM",
  var = "afm",
  id_x = 1,
  data_type = c("OS, i=1", "TM, i=1"),
  width = 0.5e-4,
  gg = gg,
  single_name_for_fit = "",
  noline = TRUE,noribbon = TRUE,
  nudge = 0.00005, size=1.5
)
gg <- gg + theme(text = element_text(size = 15))

myplotly(gg, "", "$a^2~\\mbox{[fm$^2$]}$", "$a_{\\mu}^{\\rm HVP}(s)$",
         save_pdf = "amu__s_tuning", legend_position = c(0.01,0.25))

######################################################
gg <- plot_fit(
  basename = "/home/garofalo/analysis/g-2_new_stat/fit_all_charm/amu_SDpWpLDcor_3b_BOS_BTM_a4OS_a4TM",
  var = "afm",
  id_x = 1,
  data_type = c("OS, i=2", "TM, i=2"),
  width = 0.5e-4,
  gg = NULL,
  single_name_for_fit = "",
  noline = TRUE,noribbon = TRUE, size=1.5
)
gg <- plot_fit(
  basename = "/home/garofalo/analysis/g-2_new_stat/fit_all_charm_no_fine_tuning_latt/amu_SDpWpLD_3b_BOS_BTM_a4OS_a4TM",
  var = "afm",
  id_x = 1,
  data_type = c("OS, i=1", "TM, i=1"),
  width = 0.5e-4,
  gg = gg,
  single_name_for_fit = "",
  noline = TRUE,noribbon = TRUE, size=1.5,
  nudge = 0.00005
)
gg <- gg + theme(text = element_text(size = 15))

myplotly(gg, "", "$a^2~\\mbox{[fm$^2$]}$", "$a_{\\mu}^{\\rm HVP}(c)$",
         save_pdf = "amu__c_tuning", legend_position = c(0.01,0.9))
