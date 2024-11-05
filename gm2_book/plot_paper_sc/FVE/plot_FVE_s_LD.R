library(Rose)
library(ggplot2)
library(plotly)
library(knitr)



library("gridExtra")
library("ggpubr")
library(stringr)

#source("functions.R")
mydata <- c("OS", "TM")
# gg <- plot_fit(
basename <- "/home/garofalo/analysis/g-2_new_stat/fit_all_strange/amu_LDetas_3b_BOS_BTM_FVE_a4OS_a4TM"
var <- "afm"
id_x <- 1
data_type <- mydata
width <- 1e-4
gg <- NULL
labelfit <- c("", "")
# single_name_for_fit = "",
noribbon <- TRUE
noline <- TRUE
# )
width = 0.00002
size = 1.5
filed <- paste0(basename, "_fit_data.txt")
df1 <- read.table(filed, header = FALSE, fill = TRUE)
df<-NULL
df<- df1[c(1:4,7:10),]
mydata<-rep(c("$L\\sim 5.1$", "$L\\sim 7.6$", "$L\\sim 5.4$", "$L\\sim 7.6$"),2)
idy <- ncol(df) - 2

############ relative precision
smallV<- df[c(1,3,5,7),]
bigV<- df[c(2,4,6,8),]
myen<- c("B","C","B ","C ")
reg<- c("OS","OS","TM","TM")
nudge<- rep(c(0.00002,0,0.00002,0),2)
ratio<-smallV[, idy]/bigV[,idy]
pool<- (smallV[,idy]-bigV[,idy])/sqrt(smallV[,idy+1]^2+bigV[,idy+1]^2)
pool <- abs(smallV[,idy]-bigV[,idy])*erf(abs(pool)/sqrt(2))
error<- (smallV[, idy+1]/bigV[,idy]) +(smallV[, idy]/bigV[,idy]^2)*bigV[,idy+1]
gg <- myggplot()
gg<- gg + geom_hline(yintercept = 1, color="black")
gg <- gg + geom_point(
  mapping = aes(
    x = myen, y = ratio,
    color = reg,
    shape = reg,
    fill = reg
  ),
  size = size*1.5,stroke=2
)
gg <- gg + geom_errorbar(
  mapping = aes(
    x = myen , y = ratio,
    ymin = ratio - error,
    ymax = ratio + error,
    color = reg,
    shape = reg,
    fill = reg
  ),
  width = 0.1, size = size
)
ylab <- paste0("$a_{\\mu}^{\\rm HVP, LD}(s,L)/a_{\\mu}^{\\rm HVP, LD}(s,L')$")
nameout <- paste0("amu_FVE_s_LD_ratio")
# scientific_10 <- function(x) {
#   paste0(as.character(x * 10^10), "$ \\times 10^{-10}$")
# }
# gg <- gg + scale_y_continuous(label = scientific_10)
gg<- gg + theme(text = element_text(size = 15),
                # legend.position = "none",
                axis.title.x = element_blank())
gg<- gg+ scale_color_manual(values=c("#4863A0","#C04000","#4863A0","#C04000"))
fig <- myplotly(gg, "", "", ylab,
                to_print = FALSE,
                save_pdf = nameout,
                output = "PDF", 
                legend_position = c(0.1,0.2),
)


#########################################

nudge<- rep(c(0.00002,0,0.00002,0),2)
gg <- myggplot()
gg <- gg + geom_point(
  data = df,
  mapping = aes(
    x = df[, id_x] + nudge, y = df[, idy],
    color = mydata,
    shape = mydata,
    fill = mydata
  ),
  size = size,stroke=2
)

gg <- gg + geom_errorbar(
  data = df,
  mapping = aes(
    x = df[, id_x] + nudge, y = df[, idy],
    ymin = df[, idy] - df[, idy + 1],
    ymax = df[, idy] + df[, idy + 1],
    color = mydata,
    shape = mydata,
    fill = mydata
  ),
  width = width, size = size
)
ylab <- paste0("$a_{\\mu}^{\\rm HVP,LD}(s)$")
nameout <- paste0("amu_FVE_s_LD")
scientific_10 <- function(x) {
  paste0(as.character(x * 10^10), "$ \\times 10^{-10}$")
}
gg <- gg + scale_y_continuous(label = scientific_10)
gg<- gg + theme(text = element_text(size = 15))
fig <- myplotly(gg, "", "$a^2$", ylab,
  to_print = FALSE,
  save_pdf = nameout,
  output = "PDF", legend_position = c(0, 0.2)
)
