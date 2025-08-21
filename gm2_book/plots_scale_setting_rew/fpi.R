library(Rose)
library(ggplot2)
library(plotly)
library(knitr) 


#### iter 0

dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all_rew/"
namefit <- paste0(dir, "afpi_max_twist_A12_noC20_cov_unitary_iter0")
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
cat("\n\n")
dt <- make_table_fit_result(fit)
dt<- data.frame("par"=fit$P[,1], "value"= mapply(mean_print, fit$P[,2], fit$P[,3]),
                "percent"=fit$P[,3]/fit$P[,2])

kable(dt)
gg <- myggplot(repeat_color = 1)
gg <- plot_fit(
  basename = namefit,
  var = "xi",
  id_x = 5,
  data_type = c("A", "B", "C", "D", "E", "B", "C", "D"),
  width = 1e-4,
  gg = gg, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE
)

a_iter0 <- fit$P[c(1:5),2]
da_iter0 <- fit$P[c(1:5),3]
gg <- gg + geom_vline(xintercept = 0.00678723)
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$\\xi$", "$af_\\pi$", to_print = F, save_pdf = "fpi_iter0")


#################################################################################
#  s and c
#############################


dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all/"
namefit <- paste0(dir, "afpi_max_twist_A12_noC20_cor_cov_unitary")
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
cat("\n\n")
dt <- make_table_fit_result(fit)
dt<- data.frame("par"=fit$P[,1], "value"= mapply(mean_print, fit$P[,2], fit$P[,3]),
                "percent"=fit$P[,3]/fit$P[,2])

kable(dt)
gg <- myggplot(repeat_color = 1)
gg <- plot_fit(
  basename = namefit,
  var = "xi",
  id_x = 5,
  data_type = c("A", "B", "C", "D", "E", "B", "C", "D"),
  width = 1e-4,
  gg = gg, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE
)

a_gm2 <- fit$P[c(1:5),2]
da_gm2 <- fit$P[c(1:5),3]
gg <- gg + geom_vline(xintercept = 0.00678723)
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$\\xi$", "$af_\\pi$", to_print = F, save_pdf = "fpi_sc_paper")

####################################################################################
# rew analysis
######################################################################

dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all_rew/"
namefit <- paste0(dir, "afpi_max_twist_A12_noC20_rew")
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
cat("\n\n")
dt <- make_table_fit_result(fit)
# print(fit$P)
fit$P[1,1]<-"A"
fit$P[2,1]<-"B"
fit$P[3,1]<-"C"
fit$P[4,1]<-"D"
fit$P[5,1]<-"E"
dt<- data.frame("par"=fit$P[,1], "value"= fit$P[,2],
                "error"=fit$P[,3])

# dt<- data.frame("par"=fit$P[,1] )
# dt$value<- mapply(mean_print, fit$P[,2], fit$P[,3])
# dt$no_max_twist<- mapply(mean_print, fit_no_max_twist$P[,2], fit_no_max_twist$P[,3])
#print(dt)
kable(dt)
gg <- myggplot(repeat_color = 1)
gg <- plot_fit(
  basename = namefit,
  var = "xi",
  id_x = 5,
  data_type = c("A", "B", "C", "D", "E", "B", "C", "D"),
  width = 1e-4,
  gg = gg, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE
)
# filed <- paste0(namefit, "_fit_data.txt")
# df <- read.table(filed, header = FALSE, fill = TRUE)
# idy <- ncol(df) - 2
# kable(df[, c(1,2,3,4,idy,idy+1)],col.names =c("mu","aMpi","afpi","L","afpi-inf","err") )
gg <- gg + geom_vline(xintercept = 0.00678723)
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$\\xi$", "$af_\\pi$", to_print = F, save_pdf = "fpi_rew")

############## scaling with tau         ########################################
a_tau <- a_gm2
da_tau <- da_gm2#c(0.000535517, 0.00004, 0.00008, 0.00006, 0.00006)
gg <- myggplot()

diff <- (a_tau - fit$P[c(1:5), 2])
# gg <- gg + geom_pointrange(aes(
#   x = a_tau^2, y = diff,
#   ymin = diff - da_tau,
#   ymax = diff + da_tau,
#   color = "sc-paper-rew",
#   shape = "sc-paper-rew",
#   fill = "sc-paper-rew"
# ))
gg <- gg + geom_pointrange(aes(
  x = fit$P[c(1:5), 2]^2, y = a_iter0-a_iter0,
  ymin = a_iter0 - da_iter0-a_iter0,
  ymax = a_iter0 + da_iter0-a_iter0,
  color = "iter0",
  shape = "iter0",
  fill = "iter0",
  linewidth = "iter0",
  size = "iter0"
))

gg <- gg + geom_pointrange(aes(
  x = fit$P[c(1:5), 2]^2+2e-4, y = a_gm2-a_iter0,
  ymin = a_gm2 - da_gm2-a_iter0,
  ymax = a_gm2 + da_gm2-a_iter0,
  color = "sc-paper",
  shape = "sc-paper",
  fill = "sc-paper",
  linewidth = "sc-paper",
  size = "sc-paper"
))
gg <- gg + geom_pointrange(aes(
  x = fit$P[c(1:5), 2]^2+1e-4, y = fit$P[c(1:5), 2]-a_iter0,
  ymin = fit$P[c(1:5), 2] - fit$P[c(1:5), 3]-a_iter0,
  ymax = fit$P[c(1:5), 2] + fit$P[c(1:5), 3]-a_iter0,
  color = "rew",
  shape = "rew",
  fill = "rew",
  linewidth = "rew",
  size = "rew"
))
gg<- gg +  scale_size_manual(values = c(0.5, 0.5, 0.5))   
gg<- gg +guides(size = "none") 
gg<- gg +  scale_linewidth_manual(values = c(1.5, 1.5, 1.5))   

a_rew_all<-fit$P[c(1:5), 2]
da_rew_all<-fit$P[c(1:5), 3]

gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$a^2(rew)$", "$a(X)-a(iter0)$", to_print = F, 
                save_pdf = "scaling_a", xrange = c(0, 0.0085))