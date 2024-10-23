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
plot_all_fits<- function(name,  window,quark){
  
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
  source("/home/garofalo/programs/Rose/R/plot_routines.R")
  quark_full_name<-""
  if (quark=="s")quark_full_name<-"strange"
  else if (quark=="c")quark_full_name<-"charm"
  else stop("no quark selected")
  dir <- paste0("/home/garofalo/analysis/g-2_new_stat/fit_all_",quark_full_name,"/")
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
  gg <- myggplot(repeat_color = 1)
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
        
        namelegend <- "fits"
        # if (str_detect(namefit, "log")) {
        #   namelegend <- paste0(namelegend, "-log")
        # }
        # if (str_detect(namefit, "_3b_BOS_BTM")) {
        #   namelegend <- paste0(namelegend, "-4beta")
        # }
        
        # if (
        #   # (lat == "_3b_BOS_BTM" & a2 == "_a4OS_a4TM") |
        #   #   (lat == "_3b" & a2 == "") |
        #   #   (lat == "_3b" & a2 == "_a4OS_a4TM") |
        #     (lat == "_3b_BOS_BTM" & a2 == "_alogOS_alogTM") |
        #     (lat == "_3b_onlyTM" & a2 == "")
        # ) {
        gg <- plot_fit(
          basename = namefile,
          var = "afm",
          id_x = 1,
          data_type = mydata,
          width = 1e-4,
          gg = gg,
          single_name_for_fit = namelegend,
          noribbon = TRUE
        )
        # }
        count <- count + 1
      }
    }
  }
  a <- which(df$fit == "")
  if (length(a) > 0) df <- df[-a, ]
  ave_AIC <- AIC(v = df$res, err = df$err, chi2dof = df$chi2dof, dof = df$dof, npar = df$Npar, multiplicity = df$mult)
  cat("AIC (chi2 + 2Npar - Ndat)/2 = ", mean_print(ave_AIC$m, ave_AIC$dm), "\n\n")
  
  gg <- gg + geom_pointrange(aes(
    x = 0, y = ave_AIC$m, ymin = ave_AIC$m - ave_AIC$dm,
    ymax = ave_AIC$m + ave_AIC$dm, color = "AIC", shape = "AIC", fill = "AIC"
  ),stroke=1.5,linewidth=1.5)
  gg <- gg + geom_errorbar(aes(
    x = 0, y = ave_AIC$m, ymin = ave_AIC$m - ave_AIC$dm,
    ymax = ave_AIC$m + ave_AIC$dm, color = "AIC", shape = "AIC", fill = "AIC"
  ),width=0.0001)
  
  
  scientific_10 <- function(x) {
    paste0(as.character(x * 10^10), "$ \\times 10^{-10}$")
  }
  gg <- gg + scale_y_continuous(label = scientific_10)
  gg<- gg + theme(text = element_text(size = 15))
  ylab<-paste0("$a_{\\mu}^{\\rm HVP, ",window,"}(",quark,")$")
  if (window=="")
    ylab<-paste0("$a_{\\mu}^{\\rm HVP}(",quark,")$")
  nameout<-paste0("amu_",window,"_",quark)
  fig <- myplotly(gg, "", "$a^2$", ylab, to_print = FALSE,
                  save_pdf = nameout,
                  output = "PDF", legend_position = c(0, 1))
  # fig <- fig + ylim(5.1e-9, 5.43e-9) +xlim(0,0.008)+ theme(axis.title.y=element_blank(),
  #         axis.text.y=element_blank(),
  #         axis.ticks.y=element_blank())+theme(plot.margin = margin(0,0,0,0, "cm"))
  gg <- ggplot() +
    geom_histogram(aes(y = df$res, weight = ave_AIC$AIC)) +
    theme_bw()
  
  gg <- gg +
    # ylim(5.1e-9, 5.43e-9) +
    ylab("$a_{\\mu}(s)$") + theme(plot.margin = margin(0, 0, 0, 0, "cm")) + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
    )
  print(gg)
  # ggarrange(gg, fig,
  #   align = "h",
  #   # heights = c(2, 0.7),
  #   widths = c(1, 3.5),
  #   ncol = 2, nrow = 1
  # )
  df$AIC <- ave_AIC$AIC
} 
