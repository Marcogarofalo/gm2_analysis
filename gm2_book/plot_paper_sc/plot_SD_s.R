library(Rose)
library(ggplot2)
library(plotly)
library(knitr)



library("gridExtra")
library("ggpubr")
library(stringr)

source("functions.R")
df<-plot_all_fits("SDcor","SD","s",mybinwidth=5e-13, myyrange=c(8.75e-10,9.35e-10))
df<-plot_all_fits("SDcor","SD","c",mybinwidth=5e-13, myyrange=c(11.2e-10,12.95e-10))
#plot_all_fits("SD","SD","c",quark_full_name="charm_no_fine_tuning_latt")

#df<-plot_all_fits("Wcor","W","s",mybinwidth=5e-13, myyrange=c(26.78e-10,27.6e-10))
df<-plot_all_fits("Wcor","W","s",mybinwidth=5e-13, myyrange=c(26e-10,28.25e-10))
df<-plot_all_fits("Wcor","W","c",mybinwidth=5e-13, myyrange=c(2.3e-10,4.5e-10))

df<-plot_all_fits("LDcor","LD","s",mybinwidth=5e-13, myyrange=c(16.0e-10,17.8e-10))
df<-plot_all_fits("LDcor","LD","c",  mybinwidth=1e-14, myyrange=c(1.25e-12,2.1e-12))

df<-plot_all_fits("SDpWpLDcor","","s",mybinwidth=5e-13, myyrange=c(51.5e-10,54.5e-10), leg_pos= c(0.1, 0.3))
df<-plot_all_fits("SDpWpLDcor","","c",mybinwidth=5e-13, myyrange=c(14.1e-10,17.1e-10))
