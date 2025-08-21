library(Rose)
library(ggplot2)
library(plotly)
library(knitr) 

df<-read.table("corr_difference_a.txt", header = T)

gg<- myggplot()

gg <- gg + geom_pointrange(aes(
  x = df[, 2]^2, y = df[,6],
  ymin = df[,6] - df[,7],
  ymax = df[,6] + df[,7],
  color = "rew-(sc-paper)",
  shape = "rew-(sc-paper)",
  fill = "rew-(sc-paper)",
  linewidth = "rew-(sc-paper)",
  size = "rew-(sc-paper)"
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
fig <- myplotly(gg, "", "$a^2(rew)$", "$a(rew)-a(\\mbox{sc-paper})[\\mbox{fm}]$", to_print = F,
                save_pdf = "a_corr_diff", xrange = c(0, 0.0085), legend_position = c(0.2,0.9))