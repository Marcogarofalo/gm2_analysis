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
plot_all_fits <- function(name, Wind, quark, myyrange = NULL, quark_full_name = NULL, to_plot = TRUE,
                          mybinwidth = 1e-12, leg_pos = c(0.1, 1)) {
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
  source("/home/garofalo/programs/Rose/R/plot_routines.R")
  if (is.null(quark_full_name)) {
    if (quark == "s") quark_full_name <- "strange"
    if (quark == "c") quark_full_name <- "charm"
  }
  dir <- paste0("/home/garofalo/analysis/g-2_new_stat/fit_all_", quark_full_name, "/")
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
  gg <- myggplot()
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

        namelegend <- "fit"
        if (str_detect(namefit, "log")) {
          namelegend <- paste0(namelegend, "-log")
        }
        if (str_detect(namefit, "_3b_BOS_BTM")) {
          namelegend <- paste0(namelegend, "-4beta")
        }

        count <- count + 1
      }
    }
  }
  a <- which(df$fit == "")
  if (length(a) > 0) {
    print("eliminating:")
    print(df$fit[a])
    df <- df[-a, ]
  }
  ave_AIC <- AIC2(v = df$res, err = df$err, chi2dof = df$chi2dof, dof = df$dof, npar = df$Npar, multiplicity = df$mult)
  if (to_plot) {
    df$AIC <- ave_AIC$AIC
    # for (i in c(1, 2, 50)) {
    for (i in seq_along(df$fit)) {
      namefile <- paste0(dir, df$fit[i])
      mydata <- c("OS", "TM")
      if (str_detect(df$fit[i], "onlyOS")) mydata <- c("OS")
      if (str_detect(df$fit[i], "onlyTM")) mydata <- c("TM")
      namelegend <- "fit"
      gg <- plot_fit(
        basename = namefile,
        var = "afm",
        id_x = 1,
        data_type = mydata,
        width = 1e-4,
        gg = gg,
        single_name_for_fit = namelegend,
        noribbon = TRUE, alpha_line = 0.05
      )
    }

    ###########################

    gg <- gg + geom_pointrange(aes(
      x = 0, y = ave_AIC$m, ymin = ave_AIC$m - ave_AIC$dm,
      ymax = ave_AIC$m + ave_AIC$dm, color = "AIC", shape = "AIC", fill = "AIC"
    ), linewidth = 1.5)

    # fig <- myplotly(gg, "", "$a^2$", "$a_{\\mu}(s)$",
    #   to_print = FALSE, output = "PDF",
    #   save_pdf = "tmp",
    #   legend_position = c(0.1, 1)
    # )
    # fig <- fig + ylim(5.1e-9, 5.43e-9) +xlim(0,0.008)+ theme(axis.title.y=element_blank(),
    #         axis.text.y=element_blank(),
    #         axis.ticks.y=element_blank())+theme(plot.margin = margin(0,0,0,0, "cm"))
    # lim<-print(layer_scales(gg)$y$range$range)
    dfres<-df$res
    #### add point to send the density to zero
    dfres<- c(dfres,min(dfres)/2,max(dfres)*2 )
    myaic<-c(ave_AIC$AIC,0,0)
    gg1 <- myggplot() +ylim(0,1e-9)+
      geom_density(adjust = 0.3,aes(
        y = dfres, weight = myaic,
        # fill="fits", shape="fits", colour="fits"
      ))
    # +
    # geom_density(aes(y = df$res, weight = ave_AIC$AIC))

    gg1 <- gg1 +
      ylab("$a_{\\mu}(s)$") + theme(plot.margin = margin(0.2, 0, 0, 0.5, "cm")) + theme(
        # axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1.5, hjust = 1.5)
      )
    # gg <- gg +
    #    theme(plot.margin = margin(-20, 0, 0, 0, "cm")) + theme(
    #      axis.title.y = element_blank()
    #     # ,axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
    #   )
    # gg2<-ggarrange(gg1, gg,
    #   align = "h",
    #   # heights = c(2, 0.7),
    #   widths = c(1, 3.5),
    #   ncol = 2, nrow = 1  )
    scientific_10 <- function(x) {
      # gsub("$e", " x 10^", scientific_format(scale = 10^{-10})(x),"$")
      paste0(as.character(x * 10^10), "$ \\times 10^{-10}$")
    }
    get_limit <- function(gg) {
      lim <- (layer_scales(gg)$y$range$range)
      return(lim)
    }
    lim1 <- c(get_limit(gg)[1], get_limit(gg)[2])
    if (is.null(myyrange)) {
      # lim <- (layer_scales(gg)$y$range$range)
      gg <- gg + scale_y_continuous(label = scientific_10)
      gg1 <- gg1 + scale_y_continuous(label = scientific_10)
      gg1 <- gg1 + ylim(lim1)
      # gg <- gg + ylim(lim1)
      # gg1 <- gg1 + scale_y_continuous(label = scientific_10)
      # gg1 <- gg1 + ylim(lim)
      # gg1 <- gg1 + scale_y_continuous(label = scientific_10)
      # gg1 <- gg1 + ylim(lim)
    } else {
      #gg1 <- gg1 + ylim(myyrange)
      gg <- gg + scale_y_continuous(label = scientific_10)
      gg1 <- gg1 + scale_y_continuous(label = scientific_10)
      #### xscale
      p<-ggplot_build(gg1)
      xmin<-min(p$data[[1]]$x)
      xmax<-max(p$data[[1]]$x)
      gg1 <- gg1 + coord_cartesian(ylim = myyrange, xlim = c(xmin,xmax))
      ####
      gg <- gg + coord_cartesian(ylim = myyrange)
    }
    gg <- gg + theme(text = element_text(size = 15))
    gg1 <- gg1 + theme(text = element_text(size = 15))
    

    name_out <- paste0("amu_", Wind, "_", quark)
    myylab <- paste0("$a_{\\mu}^{HPV,", Wind, "}(", quark, ")$")
    if (Wind == "") myylab <- paste0("$a_{\\mu}^{HPV}(", quark, ")$")
    fig1 <- myplotly(gg1, "", "$\\,^{\\,}$", myylab,
      to_print = FALSE, output = "PDF",
      legend_position = leg_pos,
      save_pdf = FALSE,#paste0(name_out, "_hist"),
      # width = 300
    )
    fig <- myplotly(gg, "", "$a^2$", myylab,
      to_print = FALSE, output = "PDF",
      save_pdf = FALSE, #name_out,
      legend_position = leg_pos,
    )
    fig <- fig + theme(
      plot.margin = margin(
        t = 0, # Top margin
        r = 0, # Right margin
        b = 0, # Bottom margin
        l = 0
      ), # Left margin
      # axis.title.y = element_text(
      #   margin = margin(t = 0, r = -90, b = 0, l = 0),
      # ),
      axis.title.y = element_blank(),
      axis.text.y = element_blank() 
    )
    fig1 <- fig1 + theme(
      plot.margin = margin(
        t = 0, # Top margin
        r = -1, # Right margin
        b = 0, # Bottom margin
        l = 10
      ), # Left margin
      axis.title.y = element_text(margin = margin(t = 0, r = -60, b = 0, l = 0))
    )
    legend_title <- NULL
    gg <- gg + ggplot2::labs(
      color = legend_title,
      fill = legend_title,
      shape = legend_title,
      linewidth = legend_title
    )

    # gg <- gg + ggplot2::theme(
    #   panel.background     = ggplot2::element_rect(fill = "white", color = NA),
    #   legend.key           = element_rect(fill = NA, color = NA),
    #   legend.background    = element_rect(fill = "#ffffff80", color = NA, linewidth = 1),
    #   legend.justification = c(1, 1),
    #   legend.box.margin    = margin(1, 0, 0, 0),
    # )
    # gg<-gg+theme(
    #   plot.margin = margin(t = 0,  # Top margin
    #                        r = 0,  # Right margin
    #                        b = 0,  # Bottom margin
    #                        l = 100)) # Left margin
    gg2 <- ggarrange(fig1, fig,
      align = "h",
      # heights = c(2, 0.7),
      widths = c(120, 300),
      ncol = 2, nrow = 1
    )
    # gg2 <- gg2 + theme(text = element_text(size = 15))
    
    # fig <- myplotly(gg2, "", "$a^2$", myylab,
    #                 to_print = FALSE, output = "PDF",
    #                 save_pdf = paste0(name_out, "_combined")
    #                 # legend_position = leg_pos,
    # )
    texfile <- paste0(name_out, "_combined", ".tex")
    width <- 680
    height <- 480
    tikzDevice::tikz(texfile,
      standAlone = TRUE,
      width = width / 100,
      height = height / 100
    )
    plot(gg2)
    dev.off()
    tools::texi2dvi(texfile, texi2dvi = "pdflatex", pdf = TRUE)
  }

  return(ave_AIC)
}
