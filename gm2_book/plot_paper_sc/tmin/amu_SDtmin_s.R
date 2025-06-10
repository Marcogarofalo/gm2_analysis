library(Rose)
library(ggplot2)
library(plotly)
library(knitr)
library(stringr)

eq_op <- function(x) {
  if (x == 0) {
    return("OS")
  } # equal
  else {
    return("TM")
  } # opposite
}
eq_op1 <- function(x) {
  if (x == 0) {
    return("OS-a(tau)")
  } # equal
  else {
    return("TM-a(tau)")
  } # opposite
}
AIC <- function(v, err, chi2dof, dof, npar, multiplicity = 1) {
  l <- list()
  Nmeas <- dof + npar
  AIC <- exp(-0.5 * (chi2dof * dof + 2 * npar - Nmeas)) / multiplicity
  N <- sum(AIC)
  AIC <- AIC / N
  l$AIC <- AIC
  l$m <- sum(v * AIC)
  stat <- sum(AIC * err^2)
  syst <- sum(AIC * (v - l$m)^2)
  l$dm <- sqrt(stat + syst)
  l$stat <- sqrt(stat)
  l$syst <- sqrt(syst)
  return(l)
}
AIC2 <- function(v, err, chi2dof, dof, npar, multiplicity = 1) {
  l <- list()
  Nmeas <- dof + npar
  AIC <- exp(-0.5 * (chi2dof * dof + 2 * npar - 2 * Nmeas)) / multiplicity
  N <- sum(AIC)
  AIC <- AIC / N
  l$AIC <- AIC
  l$m <- sum(v * AIC)
  stat <- sum(AIC * err^2)
  syst <- sum(AIC * (v - l$m)^2)
  l$dm <- sqrt(stat + syst)
  l$stat <- sqrt(stat)
  l$syst <- sqrt(syst)
  return(l)
}
plot_the_Ambaradam <- function(name) {
  list_lat <- c(
    "_3b", "_3b_BOS", "_3b_BTM", "_3b_BOS_BTM",
    "_3b_noC", "_3b_noC_BOS", "_3b_noC_BTM",
    "_3b_onlyOS", "_3b_onlyTM", "_4b_onlyOS", "_4b_onlyTM",
    "_3b_noC_onlyOS", "_3b_noC_onlyTM"
  )
  list_a <- c(
    "", "_a4OS", "_a4TM", "_a4OS_a4TM",
    "_alogOS", "_alogTM", "_alogOS_alogTM",
    "_alog2OS", "_alog2TM", "_alog2OS_alog2TM",
    "_alog3OS", "_alog3TM", "_alog3OS_alog3TM",
    "_a4",
    "_alog", "_alog2", "_alog3"
  )
  source("/home/garofalo/programs/Rose/R/read_block.R")
  dir <- "/home/garofalo/analysis/g-2_new_stat/fit_all_strange/"
  count <- 0
  for (lat in list_lat) {
    for (a2 in list_a) {
      namefit <- paste0("amu_", name, lat, a2)
      namefile <- paste0(dir, namefit)
      
      file <- paste0(namefile, "_fit_P.dat")
      if (file.exists(file)) {
        count <- count + 1
      }
    }
  }
  df <- data.frame(
    "fit" = rep("", count),
    "res" = rep(0, count),
    "err" = rep(0, count),
    "chi2dof" = rep(0, count),
    "dof" = rep(0, count),
    "Npar" = rep(0, count),
    "Ndat" = rep(0, count),
    "mult" = rep(0, count)
  )
  count <- 1
  for (lat in list_lat) {
    # for (lat in c("_3b_BOS_BTM" )){
    for (a2 in list_a) {
      # for (a2 in c( "_a4OS_a4TM")) {
      namefit <- paste0("amu_", name, lat, a2)
      namefile <- paste0(dir, namefit)
      
      file <- paste0(namefile, "_fit_P.dat")
      if (file.exists(file)) {
        fit <- read_fit_P_file(file)
        mydata <- c("OS", "TM")
        if (str_detect(namefit, "onlyOS")) mydata <- c("OS")
        if (str_detect(namefit, "onlyTM")) mydata <- c("TM")
        # if (fit$dof == 1) next
        df[count, 1] <- namefit
        df[count, c(2, 3)] <- fit$P[1, c(2, 3)]
        df[count, 4] <- fit$chi2dof
        df[count, 5] <- fit$dof
        df[count, 6] <- fit$npar
        df[count, 7] <- fit$ndata
        if (str_detect(namefit, "log")) {
          df[count, 8] <- 3
        } else {
          df[count, 8] <- 1
        }
        
        count <- count + 1
      }
    }
  }
  a <- which(df$fit == "")
  if (length(a) > 0) df <- df[-a, ]
  ave_AIC <- AIC2(v = df$res, err = df$err, chi2dof = df$chi2dof, dof = df$dof, npar = df$Npar, multiplicity = df$mult)
  
  df$AIC <- ave_AIC$AIC
  return(ave_AIC)
}

AIC_res <- plot_the_Ambaradam("SD")
df_ref <- data.frame("t_ref" = rep(0, 5), "amu_ref" = rep(0, 5), "damu_ref" = rep(0, 5))

df_sub <- data.frame("t_ref" = 0, "amu_ref" = AIC_res$m, "damu_ref" =AIC_res$dm)
t<-1
df_sub$t_ref[t] <- (t-1) * 0.07951
df_sub$amu_ref[t] <- AIC_res$m
df_sub$damu_ref[t] <- AIC_res$dm

for (t in c(1:5)) {
  AIC_res <- plot_the_Ambaradam(paste0("SDtmin", t - 1))
  df_ref$t_ref[t] <- (t-1) * 0.07951/2 +0.07951
  df_ref$amu_ref[t] <- AIC_res$m
  df_ref$damu_ref[t] <- AIC_res$dm
}
df_tilde <- df_ref
df_tilde$amu_ref<- df_ref$amu_ref

gg <- myggplot(fill = FALSE) 
gg <- gg + geom_point(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref,
  color = "LQCD", shape="LQCD"
), stroke=2)
gg<- gg + geom_errorbar(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref, ymin = !!(df_tilde$amu_ref - df_tilde$damu_ref),
  ymax = !!(df_tilde$amu_ref + df_tilde$damu_ref), color = "LQCD", shape="LQCD"
), width=0.005, linewidth=1.5)


df_pert<-read.table("/home/garofalo/analysis/gm2_analysis/gm2_book/amu_SD_s_pert_0.65.txt",header = TRUE)
df_tilde$amu_ref<- df_pert$amu_.s._pert
gg <- gg + geom_point(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref,
  color = "pQCD 0.65", shape="pQCD 0.65"
), stroke=2)
gg<- gg + geom_errorbar(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref, ymin = !!(df_tilde$amu_ref - df_tilde$damu_ref),
  ymax = !!(df_tilde$amu_ref + df_tilde$damu_ref), color = "pQCD 0.65", shape="pQCD 0.65"
), width=0.005, linewidth=1.5)

df_tilde$amu_ref<- df_ref$amu_ref+df_pert$amu_.s._pert
gg <- gg + geom_point(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref,
  color = "LQCD + pQCD 0.65", shape="LQCD + pQCD 0.65"
), stroke=2)
gg<- gg + geom_errorbar(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref, ymin = !!(df_tilde$amu_ref - df_tilde$damu_ref),
  ymax = !!(df_tilde$amu_ref + df_tilde$damu_ref), color = "LQCD + pQCD 0.65", shape="LQCD + pQCD 0.65"
), width=0.005, linewidth=1.5)

df_pert<-read.table("/home/garofalo/analysis/gm2_analysis/gm2_book/amu_SD_s_pert.txt",header = TRUE)
df_tilde$amu_ref<- df_pert$amu_.s._pert
gg <- gg + geom_point(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref,
  color = "pQCD", shape="pQCD"
), stroke=2)
gg<- gg + geom_errorbar(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref, ymin = !!(df_tilde$amu_ref - df_tilde$damu_ref),
  ymax = !!(df_tilde$amu_ref + df_tilde$damu_ref), color = "pQCD", shape="pQCD"
), width=0.005, linewidth=1.5)


df_tilde$amu_ref<- df_ref$amu_ref+df_pert$amu_.s._pert
gg <- gg + geom_point(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref,
  color = "LQCD + pQCD", shape="LQCD + pQCD"
), stroke=2)
gg<- gg + geom_errorbar(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref, ymin = !!(df_tilde$amu_ref - df_tilde$damu_ref),
  ymax = !!(df_tilde$amu_ref + df_tilde$damu_ref), color = "LQCD + pQCD", shape="LQCD + pQCD"
), width=0.005, linewidth=1.5)

df_pert<-read.table("/home/garofalo/analysis/gm2_analysis/gm2_book/amu_SD_s_pert_barMS.txt",header = TRUE)
df_tilde <- df_ref

df_tilde$amu_ref<- df_pert$amu_.s._pert
gg <- gg + geom_point(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref,
  color = "pQCD $\\overline{\\mbox{MS}}$", shape="pQCD $\\overline{\\mbox{MS}}$"
), stroke=2)
gg<- gg + geom_errorbar(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref, ymin = !!(df_tilde$amu_ref - df_tilde$damu_ref),
  ymax = !!(df_tilde$amu_ref + df_tilde$damu_ref), color = "pQCD $\\overline{\\mbox{MS}}$", shape="pQCD $\\overline{\\mbox{MS}}$"
), width=0.005, linewidth=1.5)


df_tilde$amu_ref<- df_ref$amu_ref+df_pert$amu_.s._pert
gg <- gg + geom_point(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref,
  color = "LQCD + pQCD $\\overline{\\mbox{MS}}$", shape="LQCD + pQCD $\\overline{\\mbox{MS}}$"
), stroke=2)
gg<- gg + geom_errorbar(aes(
  x = !!df_tilde$t_ref, y = !!df_tilde$amu_ref, ymin = !!(df_tilde$amu_ref - df_tilde$damu_ref),
  ymax = !!(df_tilde$amu_ref + df_tilde$damu_ref), color = "LQCD + pQCD $\\overline{\\mbox{MS}}$", shape="LQCD + pQCD $\\overline{\\mbox{MS}}$"
), width=0.005, linewidth=1.5)


gg<- gg +  geom_pointrange(aes(
  x = df_sub$t_ref, y = df_sub$amu_ref, ymin = df_sub$amu_ref - df_sub$damu_ref,
  ymax = df_sub$amu_ref + df_sub$damu_ref, color = "SD ",shape ="SD "
), stroke=2)
gg<- gg +  geom_errorbar(aes(
  x = df_sub$t_ref, y = df_sub$amu_ref, ymin = df_sub$amu_ref - df_sub$damu_ref,
  ymax = df_sub$amu_ref + df_sub$damu_ref, color = "SD ",shape ="SD "
), width=0.005, linewidth=1.5)





gg<- gg +  geom_pointrange(aes(
  x = df_sub$t_ref, y = df_sub$amu_ref, ymin = df_sub$amu_ref - df_sub$damu_ref,
  ymax = df_sub$amu_ref + df_sub$damu_ref, color = "SD ",shape ="SD "
), stroke=2)
gg<- gg +  geom_errorbar(aes(
  x = df_sub$t_ref, y = df_sub$amu_ref, ymin = df_sub$amu_ref - df_sub$damu_ref,
  ymax = df_sub$amu_ref + df_sub$damu_ref, color = "SD ",shape ="SD "
), width=0.005, linewidth=1.5)
gg<- gg + geom_hline(yintercept = df_sub$amu_ref - df_sub$damu_ref,linetype="dashed",color="blue")
gg<- gg + geom_hline(yintercept = df_sub$amu_ref + df_sub$damu_ref,linetype="dashed",color="blue")
gg<- gg + theme(legend.position = "none")
gg <- gg + theme(text = element_text(size = 15))

fig <- myplotly(gg, "", "$t_\\mathrm{min}$~[fm]", "$a_\\mu^{\\rm HVP, SD}(s,t_\\mathrm{min})$", to_print = FALSE,
                save_pdf = "amu_SDtmin_s", legend_position = c(0.2,0.5),
                yrange = c(0,1.6)*10e-10 )
