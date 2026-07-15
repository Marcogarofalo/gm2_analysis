library(Rose)
library(ggplot2)
library(plotly)
library(knitr)


library("gridExtra")
library("ggpubr")
library(stringr)

source("functions.R")
df <- plot_all_fits("SD", "SD", "s", mybinwidth = 5e-13, myyrange = c(8.75e-10, 9.35e-10), fpi = "130.5")

df <- plot_all_fits("W", "W", "s", mybinwidth = 5e-13, myyrange = c(26e-10, 28.25e-10), fpi = "130.5")


df <- plot_all_fits("LD", "LD", "s", mybinwidth = 5e-13, myyrange = c(16.0e-10, 17.8e-10), fpi = "130.5")
df <- plot_all_fits("SDpWpLD", "", "s", mybinwidth = 5e-13, myyrange = c(51.5e-10, 54.5e-10), fpi = "130.5", leg_pos = c(0.1, 0.3))


df <- plot_all_fits("SDpWpLD", "SDpWpLD", "c",
  mybinwidth = 5e-13,
  myyrange = c(14.1e-10, 17.1e-10),
  fpi = "130.5"
)

df <- plot_all_fits("SD", "SD", "c",
  mybinwidth = 5e-13,
  myyrange = c(11.2e-10, 12.95e-10),
  fpi = "130.5"
)

df <- plot_all_fits("LD", "LD", "c",
  mybinwidth = 5e-13,
  myyrange = c(1.25e-12, 2.1e-12),
  fpi = "130.5"
)

df <- plot_all_fits("W", "W", "c",
  mybinwidth = 5e-13,
  myyrange = c(2.3e-10, 4.5e-10),
  fpi = "130.5"
)
