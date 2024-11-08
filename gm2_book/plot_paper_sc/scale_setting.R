library(Rose)
library(ggplot2)
library(plotly)
library(knitr)
source("/home/garofalo/programs/Rose/R/plot_routines.R")

dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all/"
namefit <- paste0(dir, "afpi_max_twist_A12_noC20_cov_unitary_iter0")
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
cat("\n\n")
dt <- make_table_fit_result(fit)
# dt<- data.frame("par"=fit$P[,1] )
# dt$value<- mapply(mean_print, fit$P[,2], fit$P[,3])
# dt$no_max_twist<- mapply(mean_print, fit_no_max_twist$P[,2], fit_no_max_twist$P[,3])
# print(dt)
print(dt)
gg <- myggplot(repeat_color = 1)
gg <- plot_fit(
  basename = namefit,
  var = "xi",
  id_x = 5,
  data_type = c("A, i=1", "B, i=1", "C, i=1", "D, i=1", "E, i=1"),
  width = 2e-4,
  gg = gg,
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, size = 2, alpha_ribbon = 0.2
)
gg <- gg + theme(text = element_text(size = 15))

gg <- gg + geom_vline(xintercept = 6.77683303196e-3)
fig <- myplotly(gg, "", "$\\xi_\\pi$", "$aF_\\pi$",
  to_print = FALSE,
  save_pdf = "fpi_vs_xi", xrange = c(0.007, 0.034),
  legend_position = c(0.8, 0.65)
)
###################################

gg <- myggplot(repeat_color = 1)
# gg <- plot_fit(
#   basename = namefit,
#   var = "xi",
#   id_x = 5,
#   data_type = c("A,i=2", "B,i=2", "C,i=2", "D,i=2", "E,i=2"),
#   width = 0.3e-4,
#   gg = gg,
#   labelfit = "",
#   # , single_name_for_fit = "fit"
#   # , nolabel_for_fit = TRUE
#   noline = TRUE, size=1, alpha_ribbon = 0.2
# )

namefit <- paste0(dir, "afpi_max_twist_A12_noC20_cov_unitary_iter0")
gg <- plot_fit(
  basename = namefit,
  var = "xi",
  id_x = 5,
  data_type = c("A, i=1", "B, i=1", "C, i=1", "D, i=1", "E, i=1"),
  width = 0.3e-4,
  gg = gg,
  labelfit = "",
  nudge = 0,
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, noribbon = FALSE, size = 2, alpha_ribbon = 0.2
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
#   noline = TRUE, size=1, alpha_ribbon = 0.2
# )
xiphys <- 6.77683303196e-3
gg <- gg + geom_vline(xintercept = xiphys)
# df<-read.table(paste0(dir, "afpi_A12_noC20_cov_fit_data.txt"))
# df<-df[c(13:15),]
# iy<-dim(df)[2]-2
# lab<-c("B reweighting i=1","C reweighting i=1","D reweighting i=1")
# # the real x should be df[,5]
# gg<- gg + geom_point(aes(x=!!xiphys, y=!!df[,iy],
#                                color=!!lab, shape=!!lab, fill=!!lab))
# gg<- gg + geom_errorbar(aes(x=!!xiphys,
#                               ymin=!!(df[,iy]-df[,iy+1]),
#                               ymax=!!(df[,iy]+df[,iy+1]), color=!!lab, shape=!!lab, fill=!!lab),
#                         width = 0.5e-4)
#
#
#
# df<-read.table(paste0(dir, "afpi_max_twist_A12_noC20_cor_cov_fit_data.txt"))
# df<-df[c(13:15),]
# iy<-dim(df)[2]-2
# lab<-c("B reweighting i=2","C reweighting i=2","D reweighting i=2")
# # the real x should be df[,5]
# gg<- gg + geom_point(aes(x=xiphys+0.00002, y=df[,iy],
#                          color=lab, shape=lab, fill=lab))
# gg<- gg + geom_errorbar(aes(x=xiphys+0.00002,
#                             ymin=df[,iy]-df[,iy+1],
#                             ymax=df[,iy]+df[,iy+1], color=lab, shape=lab, fill=lab),
#                         width = 0.5e-4)
########################### à
namefit <- paste0(dir, "afpi_max_twist_A12_noC20_cor_cov_unitary")
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
afm <- fit$P[c(1, 2, 3, 4), c(2, 3)]
afpi <- afm * 130.5 / 197.326963
lab <- c("i=2", "i=2", "i=2", "i=2")
gg <- gg + geom_point(aes(
  x = !!xiphys, y = !!(afpi[, 1]),
  color = !!lab, shape = !!lab, fill = !!lab
), size=2)
gg <- gg + geom_errorbar(
  aes(
    x = !!xiphys,
    ymin = !!(afpi[, 1] - afpi[, 2]),
    ymax = !!(afpi[, 1] + afpi[, 2]), color = !!lab, shape = !!lab, fill = !!lab
  ),
  width = 0.3e-4, size=2
)

namefit <- paste0(dir, "afpi_max_twist_noA12_noC20_cor_cov_phys_only")
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
afm <- fit$P[c(1, 2, 3, 4), c(2, 3)]
afpi <- afm * 130.5 / 197.326963
lab <- c("direct approach", "direct approach", "direct approach", "direct approach")
xiphys <- 6.77683303196e-3 - 0.5e-4
gg <- gg + geom_point(aes(
  x = !!xiphys, y = !!(afpi[, 1]),
  color = !!lab, shape = !!lab, fill = !!lab
), size=2)
gg <- gg + geom_errorbar(
  aes(
    x = !!xiphys,
    ymin = !!(afpi[, 1] - afpi[, 2]),
    ymax = !!(afpi[, 1] + afpi[, 2]), color = !!lab, shape = !!lab, fill = !!lab
  ),
  width = 0.3e-4, size=2
)

####

gg <- gg + theme(text = element_text(size = 15))
fig <- myplotly(gg, "", "$\\xi_\\pi$", "$aF_\\pi$",
  to_print = FALSE,
  save_pdf = "fpi_vs_xi_zoom", xrange = c(0.0065, 0.008),
  yrange = c(0.037, 0.053),
  legend_position = c(0.8, 0.65)
)
################################################################################
#
###############################################################################

dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all/"


namefit <- paste0(dir, "aMpi2_over_afpi2_A12_noC20_cov_unitary_iter0")

amu_phys <- c( # 0.000748072311508,
  0.000666857356526,
  0.000586436819304,
  0.000493352338666
  # 0.000430647061354
)

gg <- plot_fit(
  basename = namefit,
  var = "amu",
  id_x = 1,
  data_type = c("A, i=1", "B, i=1", "C, i=1", "D, i=1", "E, i=1"), # i=1
  width = 0.7e-4,
  gg = NULL,
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, size = 2, alpha_ribbon = 0.2
)
gg <- gg + theme(text = element_text(size = 15))

gg <- gg + geom_hline(yintercept = 1.07015457788)
fig <- myplotly(gg, "", "$a\\mu_\\ell$", "$(M_\\pi/F_\\pi)^2$",
  to_print = FALSE,
  save_pdf = "Mpi_over_fpi", yrange = c(0, 5.5),
  legend_position = c(-0.0, 0.99)
)
#########

#
# gg <- plot_fit(
#   basename = namefit,
#   var = "amu",
#   id_x = 1,
#   data_type = c("A,i=2", "B,i=2", "C,i=2", "D,i=2", "E,i=2", "B,i=2", "C,i=2", "D,i=2"),
#   width = 0.1e-4,
#   gg = NULL,
#   labelfit = "",
#   # , single_name_for_fit = "fit"
#   # , nolabel_for_fit = TRUE
#   noline = TRUE, size=1, alpha_ribbon = 0.2
# )
gg <- NULL
# namefit <- paste0(dir, "aMpi2_over_afpi2_A12_noC20_cor_cov_unitary")
namefit <- paste0(dir, "aMpi2_over_afpi2_A12_noC20_cov_unitary_iter0")
gg <- plot_fit(
  basename = namefit,
  var = "amu",
  id_x = 1,
  data_type = c("A, i=1", "B, i=1", "C, i=1", "D, i=1", "E, i=1"), # i=1
  width = 0.08e-4,
  gg = gg,
  labelfit = "",
  nudge = 0,
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE, noribbon = FALSE, size = 2, alpha_ribbon = 0.2
)

##############################

# df <- read.table(paste0(dir, "aMpi2_over_afpi2_A12_noC20_cov_fit_data.txt"))
# lab <- c("B reweighting i=1", "C reweighting i=1", "D reweighting i=1")
# df1 <- df[c(8, 9, 11), ]
# df <- df[c(13:15), ]
# iy <- dim(df)[2] - 2
# Mpi_fpi_phys <- 1.07015457788
# # the real x should be df[,1]
# mu2 <- df[, 1]
# Mf2 <- df[, iy]
# 
# mu1 <- df1[, 1]
# Mf1 <- df1[, iy]
# deriv <- (Mf2 - Mf1) / (mu2 - mu1)
# Mpi_fpi_phys <- Mf2 - (mu2 - amu_phys) * deriv
# amu_phys<-mu2
# gg<- gg + geom_point(aes(x=!!amu_phys, y=!!Mpi_fpi_phys,
#                          color=!!lab, shape=!!lab, fill=!!lab),size=1)
# gg<- gg + geom_errorbar(aes(x=!!amu_phys,
#                             ymin=!!(Mpi_fpi_phys-df[,iy+1]),
#                             ymax=!!(Mpi_fpi_phys+df[,iy+1]),
#                             color=!!lab, shape=!!lab, fill=!!lab),
#                         width = 0.08e-4,size=1)

# df <- read.table(paste0(dir, "aMpi2_over_afpi2_A12_noC20_cor_cov_fit_data.txt"))
# df <- df[c(13:15), ]
# iy <- dim(df)[2] - 2
# lab <- c("B reweighting i=2", "C reweighting i=2", "D reweighting i=2")
# the real x should be df[,1]


# gg<- gg + geom_point(aes(x=!!df[,1], y=!!df[,iy],
#                          color=!!lab, shape=!!lab, fill=!!lab),size=1)
# gg<- gg + geom_errorbar(aes(x=!!df[,1],
#                             ymin=!!(df[,iy]-df[,iy+1]),
#                             ymax=!!(df[,iy]+df[,iy+1]), color=lab, shape=lab, fill=lab),
#                         width = 0.08e-4,size=1)
#

# gg<- gg + geom_point(aes(x=!!amu_phys+0.2e-5, y=!!Mpi_fpi_phys,
#                          color=lab, shape=lab, fill=lab),size=1)
# gg<- gg + geom_errorbar(aes(x=!!amu_phys+0.2e-5,
#                             ymin=!!(Mpi_fpi_phys-df[,iy+1]),
#                             ymax=!!(Mpi_fpi_phys+df[,iy+1]), color=lab, shape=lab, fill=lab),
#                         width = 0.08e-4,size=1)
Mpi_fpi_phys <- 1.07015457788
# 
name <- "aMpi2_over_afpi2_A12_noC20_cor_cov_unitary"
df <- read.table(paste0(dir, name, "_amul_res.txt"))
lab <- c("i=2","i=2", "i=2", "i=2", "i=2")
gg <- gg + geom_point(aes(
  x = !!df[,2] , y = !!Mpi_fpi_phys,
  color = !!lab, shape = !!lab, fill = !!lab
), size = 2)
gg <- gg + geom_errorbarh(
  aes(
    y = !!Mpi_fpi_phys ,
    xmin = !!(df[,2] - df[, 3]),
    xmax = !!(df[,2] + df[, 3]), color = !!lab, shape = !!lab, fill = !!lab
  ),
  height = 0.003, size = 2
)


name <- "aMpi2_over_afpi2_A12_noC20_cor_cov_phys_only"
df <- read.table(paste0(dir, name, "_amul_res.txt"))
lab <- c("direct approach", "direct approach", "direct approach", "direct approach")
gg <- gg + geom_point(aes(
  x = !!df[,2]-0.003 , y = !!Mpi_fpi_phys,
  color = lab, shape = lab, fill = lab
), size = 2)
gg <- gg + geom_errorbarh(
  aes(
    y = !!Mpi_fpi_phys-0.003 ,
    xmin = !!(df[,2] - df[, 3]),
    xmax = !!(df[,2] + df[, 3]), color = lab, shape = lab, fill = lab
  ),
  height = 0.003, size = 2
)




gg <- gg + theme(text = element_text(size = 15))
# colorlist <- c(
#   "#404040", "#4863A0", "#C04000",
#   "#228B22", "#8B008B", "#fc03e3", 
#   "#996600", "#999999", "#FFCC33",
#   "#FF6600", "#6633FF", "#9966FF",
#   "#006666", "#FFCCFF", "#fc0303",
#   "#03fc07", "#0335fc", "#fc03e3",
#   "#d7fc03", "#00CCFF"
# )
# colorlist <- rep(colorlist, each = 1)
# gg <- gg + scale_color_manual(values = colorlist)

gg <- gg + geom_hline(yintercept = 1.07015457788)
gg <- gg + ggplot2::theme(
  # panel.background     = ggplot2::element_rect(fill = alpha("white", 0), color = NA),
  legend.key           = element_rect(fill = NA, color = NA),
  legend.background    = element_rect(fill = "#ffffff80", color = NA, linewidth = 1),
  legend.justification = c(1, 1),
  legend.box.margin    = margin(1, 0, 0, 0),
)
legend_title <- NULL
gg <- gg + ggplot2::labs(
  color = legend_title,
  fill = legend_title,
  shape = legend_title,
  linewidth = legend_title
)
fig <- myplotly(gg, "", "$a\\mu_\\ell$", "$(M_\\pi/F_\\pi)^2$",
  to_print = FALSE,
  save_pdf = "Mpi_over_fpi_zoom", yrange = c(1.05, 1.159),
  xrange = c(0.00042, 0.00072),
  legend_position = c(0.12, 1.01), restyle = FALSE
)
