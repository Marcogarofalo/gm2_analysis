library(Rose)
library(ggplot2)
library(plotly)
library(knitr)



library("gridExtra")
library("ggpubr")
library(stringr)

source("functions.R")
plot_all_fits("SDetas","SD","s")
plot_all_fits("SDcor","SD","c")

plot_all_fits("Wetas","W","s")
plot_all_fits("Wcor","W","c")

plot_all_fits("LDetas","LD","s")
plot_all_fits("LDcor","LD","c")

plot_all_fits("SDpWpLDetas","","s")
plot_all_fits("SDpWpLDcor","","c")
