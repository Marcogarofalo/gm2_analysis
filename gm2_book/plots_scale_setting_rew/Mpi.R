library(Rose)
library(ggplot2)
library(plotly)
library(knitr) 


dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all_rew/"
name<-"aMpi2_over_afpi2_A12_noC20_cor_cov_unitary"
namefit <- paste0(dir, name)
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
dt <- make_table_fit_result(fit)
print(dt)


gg <- plot_fit(
  basename = namefit,
  var = "amu",
  id_x = 1,
  data_type = c("A", "B", "C", "D", "E", "B", "C", "D"),
  width = 1e-4,
  gg = NULL, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE
)

# filed <- paste0(namefit, "_fit_data.txt")
# df <- read.table(filed, header = FALSE, fill = TRUE)
# idy <- ncol(df) - 2
# kable(df[, c(1,2,3,4,idy,idy+1)],col.names =c("mu","aMpi","afpi","L","afpi-inf","err") )
gg <- gg + geom_hline(yintercept = 1.07015457788)
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$a\\mu_\\ell$", "$(M_\\pi/f_\\pi)^2$", to_print = F, save_pdf = "Mpi_sc_paper")
df<-read.table(paste0(dir,name,"_amul_res.txt"))

# df<-read.table(paste0(dir,"aMpi2_over_afpi2_a2_A_cov_amul_res.txt"))
# df1<-read.table(paste0(dir,"aMpi2_over_afpi2_a2_A_no_max_twist_cov_amul_res.txt"))
dft<-df[,c(1,2)]
dft[,2] <- mapply(mean_print, df[,2],df[,3])
loop<- df

# dft[,3] <- mapply(mean_print, df1[,2],df1[,3])
kable(dft, col.names = c("Ens.","$a\\mu_\\ell$"))

#############################################
# rew
##############################


dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all_rew/"
name<-"aMpi2_over_afpi2_A12_noC20_rew"
namefit <- paste0(dir, name)
file <- paste0(namefit, "_fit_P.dat")
fit <- read_fit_P_file(file)
dt <- make_table_fit_result(fit)
print(dt)


gg <- plot_fit(
  basename = namefit,
  var = "amu",
  id_x = 1,
  data_type = c("A", "B", "C", "D", "E", "B", "C", "D"),
  width = 1e-4,
  gg = NULL, 
  labelfit = "",
  # , single_name_for_fit = "fit"
  # , nolabel_for_fit = TRUE
  noline = TRUE
)

# filed <- paste0(namefit, "_fit_data.txt")
# df <- read.table(filed, header = FALSE, fill = TRUE)
# idy <- ncol(df) - 2
# kable(df[, c(1,2,3,4,idy,idy+1)],col.names =c("mu","aMpi","afpi","L","afpi-inf","err") )
gg <- gg + geom_hline(yintercept = 1.07015457788)
gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$a\\mu_\\ell$", "$(M_\\pi/f_\\pi)^2$", to_print = F, save_pdf = "Mpi_rew")
df<-read.table(paste0(dir,name,"_amul_res.txt"))

# df<-read.table(paste0(dir,"aMpi2_over_afpi2_a2_A_cov_amul_res.txt"))
# df1<-read.table(paste0(dir,"aMpi2_over_afpi2_a2_A_no_max_twist_cov_amul_res.txt"))
dft<-df[,c(1,2)]
dft[,2] <- mapply(mean_print, df[,2],df[,3])
rew<-df
# dft[,3] <- mapply(mean_print, df1[,2],df1[,3])
kable(dft, col.names = c("Ens.","$a\\mu_\\ell$"))


########################################################################


gg <- myggplot()
dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all_rew/"
namefit <- paste0(dir, "afpi_max_twist_A12_noC20_rew")
file <- paste0(namefit, "_fit_P.dat")
fit_a <- read_fit_P_file(file)

gg <- gg + geom_pointrange(aes(
  x = fit_a$P[c(1:5), 2]^2+1e-4, y = loop[,2]-loop[,2],
  ymin = loop[,2] - loop[,3]-loop[,2],
  ymax = loop[,2] + loop[,3]-loop[,2],
  color = "sc-paper",
  shape = "sc-paper",
  fill = "sc-paper",
  linewidth = "sc-paper",
  size = "sc-paper"
))

gg <- gg + geom_pointrange(aes(
  x = fit_a$P[c(1:5), 2]^2, y = rew[,2]-loop[,2],
  ymin = rew[,2] - rew[,3]-loop[,2],
  ymax = rew[,2] + rew[,3]-loop[,2],
  color = "rew",
  shape = "rew",
  fill = "rew",
  linewidth = "rew",
  size = "rew"
))

colorlist <- c(
  # "#404040",
  "#4863A0", "#C04000",
  "#228B22", "#8B008B", "#00CCFF",
  "#996600", "#999999", "#FFCC33",
  "#FF6600", "#6633FF", "#9966FF",
  "#006666", "#FFCCFF", "#fc0303",
  "#03fc07", "#0335fc", "#fc03e3",
  "#d7fc03"
)

gg <- gg + scale_color_manual(values = colorlist)
gg<- gg +  scale_size_manual(values = c(0.5, 0.5, 0.5))   
gg<- gg +guides(size = "none") 
gg<- gg +  scale_linewidth_manual(values = c(1.5, 1.5, 1.5))   

gg<-gg+ theme(text = element_text(size = 20))
fig <- myplotly(gg, "", "$a^2(rew)$", "$m_\\ell(X)-m_\\ell(\\mbox{sc-paper})$", to_print = F,
                save_pdf = "ml_scaling", xrange = c(0, 0.0085))