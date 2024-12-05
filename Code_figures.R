## This script process and generate figures
## Date: 241130

## Inialize ####
rm(list = ls())
setwd(r"(D:\OneDrive - Imperial College London\PhD_220506\LAI_framework\final_revision)")
library(matrixStats)
library(raster)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(ncdf4)
library(maps)
library(reshape2)
library(scales)
library(RColorBrewer)
library(qpcR)
library(abind)
library(patchwork)
library(stringr)
library(rnaturalearth)
library(bentcableAR)
library(ggsci)

# base map
coast <- ne_coastline(scale = "medium", returnclass = "sf")
ggplot(data = coast) + geom_sf() + theme_classic()
world <- ne_countries(scale = "medium", returnclass = "sf")
mp = ggplot(data = world) +
  geom_sf(color = "white", fill = "white") +
  geom_sf(data = coast) + theme_classic()

## Figure 2 ####
load("./Data/obs_multi_avg_fapar.rda")
load("./Data/mod_multi_avg_fapar.rda")
load("./Data/limiting_region.rda")

# plot limiting region
slp_t = t(limiting_reg)
longData<-melt(slp_t)
row = 360
col=720
res=0.5
longData$Var1 <- rep(seq(-180 + res/2, 180, res),row)
longData$Var2 <- rep(seq(-90 + res/2, 90, res),each = col)
longData <- longData[is.na(longData$value)==F,]

mp2 <- mp + geom_tile(data = longData, aes(x = Var1, y = -Var2,fill=factor(value)),inherit.aes = T, show.legend = T) + 
  scale_fill_manual(name = "",
                    values =c("-1"="#87ceeb","1"="#556b2f"),
                    breaks = c("-1","1"),
                    labels=c("Water","Energy")) +
  labs(x="Longitude", y="Latitude") +
  # labs(x=NULL, y=NULL)+
  labs(title = "Energy/water-limited region") +
  scale_y_continuous(breaks=seq(90,-90,-30))+
  scale_x_continuous(breaks=seq(-180,180,60)) +
  theme_minimal() + 
  theme(legend.position = c(0.2,0.4),legend.text = element_text(size=12)) 
mp2
# plot fAPAR comparison
slp_t = t(obs_multi_avg_fapar)
longData<-melt(slp_t)
res = 0.5
row = nrow(obs_multi_avg_fapar)
col = ncol(obs_multi_avg_fapar)
longData$Var1 <- rep(seq(-180 + res/2, 180, res),row)
longData$Var2 <- rep(seq(-90 + res/2, 90, res),each = col)
longData <- longData[longData$value!=0,]

mp3 <- mp + geom_tile(data = longData, aes(x = Var1, y = -Var2,fill=value),inherit.aes = T, show.legend = T) + 
  scale_fill_stepsn(limits=c(0, 1), breaks = seq(0, 1,0.1), colors = brewer.pal(8,"GnBu")) +  # scale for fAPAR
  labs(x="Longitude", y="Latitude", fill="fAPAR") +
  labs(title = "MODIS fAPAR") +
  scale_y_continuous(breaks=seq(90,-90,-30))+
  scale_x_continuous(breaks=seq(-180,180,60)) +
  theme_minimal() + 
  theme(legend.position = "bottom", legend.key.width= unit(2, 'cm')) 

slp_t = t(mod_multi_avg_fapar)
longData<-melt(slp_t)
longData$Var1 <- rep(seq(-180 + res/2, 180, res),row)
longData$Var2 <- rep(seq(-90 + res/2, 90, res),each = col)
longData<-longData[longData$value!=0,]

mp4 <- mp + geom_tile(data = longData, aes(x = Var1, y = -Var2,fill=value),inherit.aes = T, show.legend = T) + 
  scale_fill_stepsn(limits=c(0, 1), breaks = seq(0, 1,0.1), colors = brewer.pal(8,"GnBu")) +
  labs(x="Longitude", y="Latitude", fill="fAPAR") +
  labs(title = "Modelled fAPAR") +
  scale_y_continuous(breaks=seq(90,-90,-30))+
  scale_x_continuous(breaks=seq(-180,180,60)) +
  theme_minimal() + 
  theme(legend.position = "bottom")+
  theme(legend.position = "bottom", legend.key.width= unit(2, 'cm')) 

mp5 <- mp3+mp4 +
  plot_layout(guides = 'collect', ncol=1) & theme(legend.position = 'bottom',legend.key.width= unit(1.5, 'cm'))


lt <- "
AA##
AABB
AABB
AA##
"
wrap_plots(mp5, mp2, design = lt)
# ggsave(filename = "./final_revision/figure/figure2.png",width = 9.18, height = 8.22, dpi = 300)

## Figure 3 panel 1: correlation between MODIS and modelled fAPARmax ####
fapar_comp <- data.frame(fapar_obs = as.vector(obs_multi_avg_fapar),
                         fapar_sim = as.vector(mod_multi_avg_fapar))
fapar_comp <- na.omit(fapar_comp)

mod1 <-lm(fapar_obs~fapar_sim-1, data = fapar_comp)

mp1 = ggplot(fapar_comp ,aes(x=fapar_sim, y=fapar_obs))+
  geom_pointdensity() +
  geom_smooth(method = "lm",formula=y~x-1, na.rm=T,col = "black",se = F) +
  stat_poly_eq(formula = y~x-1,method = "lm",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE,
               label.x = "right", label.y = "bottom") +
  scale_color_viridis(option = "inferno",direction = -1) +
  labs(x="Modelled fAPAR",y = "MODIS fAPAR") +
  geom_abline(slope=1, intercept = 0,col='darkgrey', linetype = "dashed")+
  theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 10))+
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_light()+
  annotate('text', x = 0.15, y = 0.9, label = paste("RMSE = ",round(RMSE(mod1),2),sep = ""), size = 4.5)
mp1


## Figure 4 ####
load("./Data/mod_fapar_trend_masked.rda")
load("./Data/obs_fapar_trend_masked.rda")

slp_t = t(obs_fapar_trend_masked)
longData<-melt(slp_t)
res = 1.5
longData$Var1 <- rep(seq(-180 + res/2, 180, res),120)
longData$Var2 <- rep(seq(-90 + res/2, 90, res),each = 240)
longData <- longData[longData$value!=0,]

td1 <- mp + geom_tile(data = longData, aes(x = Var1, y = -Var2,fill=value),inherit.aes = T, show.legend = T) + 
  scale_fill_stepsn(limits=c(-1, 1), breaks = seq(-1, 1,0.2), colors = brewer.pal(8,"BrBG")) +
  labs(x="Longitude", y="Latitude", fill='trend') +
  labs(title = "MODIS fAPAR") +
  scale_y_continuous(breaks=seq(90,-90,-30))+
  scale_x_continuous(breaks=seq(-180,180,60)) +
  theme_minimal() + 
  theme(legend.position = "bottom")+
  theme(legend.position = "bottom", legend.key.width= unit(2, 'cm')) 

slp_t = t(mod_fapar_trend_masked)
longData<-melt(slp_t)
longData$Var1 <- rep(seq(-180 + res/2, 180, res),120)
longData$Var2 <- rep(seq(-90 + res/2, 90, res),each = 240)
longData<-longData[longData$value!=0,]

td2 <- mp + geom_tile(data = longData, aes(x = Var1, y = -Var2,fill=value),inherit.aes = T, show.legend = T) + 
  scale_fill_stepsn(limits=c(-1, 1), breaks = seq(-1, 1,0.2), colors = brewer.pal(8,"BrBG")) +
  labs(x="Longitude", y="Latitude", fill='trend') +
  labs(title = "Modelled fAPAR") +
  scale_y_continuous(breaks=seq(90,-90,-30))+
  scale_x_continuous(breaks=seq(-180,180,60)) +
  theme_minimal() + 
  theme(legend.position = "bottom")+
  theme(legend.position = "bottom", legend.key.width= unit(2, 'cm')) 

td1+td2 +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

## Figure 5 first panel: correlation between P model predicted fAPARmax trend and MODIS ones ####
slp_comp <- data.frame(modis_pc = as.vector(obs_fapar_trend_masked),
                       model_pc = as.vector(mod_fapar_trend_masked))
slp_comp <- na.omit(slp_comp)

mod1 <-lm(modis_pc~model_pc-1, data = slp_comp)

ggplot(slp_comp ,aes(x=model_pc, y=modis_pc))+
  geom_pointdensity() +
  geom_smooth(method = "lm",formula=y~x-1, na.rm=T,col = "red",se = F, linetype = "dashed") +
  stat_poly_eq(formula = y~x-1,method = "lm",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE,
               label.x = "right", label.y = "bottom") +
  scale_color_viridis(option = "inferno",direction = -1) +
  labs(xlab="P model",ylab = "MODIS") +
  geom_abline(slope=1, intercept = 0,col='black')+
  theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 10))+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"))+
  scale_x_continuous(limits = c(-1,1), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1,1), expand = c(0,0)) +
  theme_light()+
  annotate('text', x = -0.55, y = 0.85, label = paste("RMSE = ",round(RMSE(mod1),2),sep = ""), size = 3.5)

## Figure 6 ####
load("./Data/fac_imp_masked.rda")
slp_t = t(imp_masked)
longData<-melt(slp_t)
res = 1.5
longData$Var1 <- rep(seq(-180 + res/2, 180, res),120)
longData$Var2 <- rep(seq(-90 + res/2, 90, res),each = 240)
longData<-longData[longData$value!=0,]

mypal = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

mp + geom_tile(data = longData, aes(x = Var1, y = -Var2,fill=factor(value)),inherit.aes = T, show.legend = T) + 
  scale_fill_manual(name = "",
                    values =mypal,
                    breaks = seq(0.9,6.3,0.9),
                    labels=c("Precipitation (+)",expression(paste(CO[2], " (+)")), 
                             "Radiation (+)", "Temperature (+)",
                             "Temperature (-)", "Radiation (-)", "Precipitation (-)")) +
  labs(x="Longitude", y="Latitude") +
  scale_y_continuous(breaks=seq(90,-90,-30))+
  scale_x_continuous(breaks=seq(-180,180,60)) +
  theme_minimal() + 
  theme(legend.position = "right", legend.key.width= unit(0.3, 'cm'),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),
        legend.text = element_text(size=10)) 

## supp figure 1 ####
bi <- mod_multi_avg_fapar - obs_multi_avg_fapar
# hist(bi)
bi_na <- na.omit(as.vector(bi))

# plot bias (0.15 as a bin)
slp_t = t(bi)
longData<-melt(slp_t)
res=0.5
row = nrow(obs_multi_avg_fapar)
col = ncol(obs_multi_avg_fapar)
longData$Var1 <- rep(seq(-180 + res/2, 180, res),row)
longData$Var2 <- rep(seq(-90 + res/2, 90, res),each = col)
longData<-longData[longData$value!=0,]

bi1 <- mp + geom_tile(data = longData, aes(x = Var1, y = -Var2,fill=value),inherit.aes = T, show.legend = T) + 
   scale_fill_stepsn(colors = brewer.pal(8,"RdBu"), oob = scales::oob_squish,
                    limits=c(-0.45, 0.5), breaks = seq(-0.45, 0.5,0.15),labels = scales::comma) +
  labs(x="Longitude", y="Latitude", fill="Differences") +
  labs(title = "fAPAR differences") +
  scale_y_continuous(breaks=seq(90,-90,-30))+
  scale_x_continuous(breaks=seq(-180,180,60)) +
  theme_minimal() + 
  theme(legend.position = "bottom")+
  theme(legend.position = "bottom", legend.key.width= unit(2, 'cm')) 
bi1

bi2 <- ggplot(data.frame(bi = bi_na), aes(x=bi)) + 
  geom_histogram(aes(y=..density..),breaks = seq(-0.45, 0.5, 0.05),binwidth = 0.05 , colour="black", fill="white")+
  geom_density(fill="darkgrey", alpha=.3) +
  geom_vline(xintercept=c(-0.15,0.15),
             color="red", linetype="dashed", size=.8)+
  geom_vline(xintercept=c(-0.1,0.1),
             color="blue", linetype="dashed", size=0.8)+
  scale_x_continuous(limits = c(-0.45, 0.5),breaks = seq(-0.45, 0.5, 0.3),
                     oob = scales::oob_squish,labels = scales::comma)+
  xlab("")+
  ylab("Density")+
  theme_minimal()+
  theme(axis.text = element_text(size=10), axis.title = element_text(size=10))
bi2

bi1 + inset_element(bi2, 0, 0.08, 0.32, 0.55)

## supp figure 2
load("./Data/mod_multi_avg_fapar_variedf0nz.rda")
load("./Data/obs_multi_avg_fapar.rda")

fapar_comp <- data.frame(fapar_obs = as.vector(obs_multi_avg_fapar),
                         fapar_sim = as.vector(mod_multi_avg_fapar_variedf0nz))
fapar_comp <- na.omit(fapar_comp)

mod1 <-lm(fapar_obs~fapar_sim-1, data = fapar_comp)

ggplot(fapar_comp ,aes(x=fapar_sim, y=fapar_obs))+
  geom_pointdensity() +
  geom_smooth(method = "lm",formula=y~x-1, na.rm=T,col = "black",se = F) +
  stat_poly_eq(formula = y~x-1,method = "lm",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE,
               label.x = "right", label.y = "bottom") +
  scale_color_viridis(option = "inferno",direction = -1) +
  labs(x="Modelled fAPAR",y = "MODIS fAPAR") +
  geom_abline(slope=1, intercept = 0,col='darkgrey', linetype = "dashed")+
  theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 10))+
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_light()+
  annotate('text', x = 0.15, y = 0.9, label = paste("RMSE = ",round(RMSE(mod1),2),sep = ""), size = 4.5)

## supp figure 3 ####
load("./Data/MODIS_derived_k.rda")
plot(k, zlim=c(0,1),xlab="Longitude", ylab="Latitude",xaxt="n",yaxt="n")

## supp figure 4 ####
load("./Data/obs_multi_avg_fapar.rda")
load("./Data/mod_multi_avg_fapar.rda")
load("./Data/mod_multi_avg_fapar_modis_pftk.rda")
load("./Data/mod_multi_avg_fapar_insitu_pftk.rda")
load("./Data/mod_multi_avg_fapar_modisk.rda")

mod <- list()
mod[[1]] <- mod_multi_avg_fapar
mod[[2]] <- mod_multi_avg_fapar_modisk
mod[[3]] <- mod_multi_avg_fapar_modis_pftk
mod[[4]] <- mod_multi_avg_fapar_insitu_pftk

p <- list()
for (i in 1:4) {
  slp_comp <- data.frame(obs_fapar = as.vector(obs_multi_avg_fapar),
                         mod_fapar = as.vector(mod[[i]]))
  slp_comp <- na.omit(slp_comp)
  
  mod1 <-lm(obs_fapar~mod_fapar-1, data = slp_comp)
  
  p[[i]] = ggplot(slp_comp ,aes(x=mod_fapar, y=obs_fapar))+
    geom_pointdensity() +
    geom_smooth(method = "lm",formula=y~x-1, na.rm=T,col = "black",se = F) +
    stat_poly_eq(formula = y~x-1,method = "lm",
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                 parse = TRUE,
                 label.x = "right", label.y = "bottom") +
    scale_color_viridis(option = "inferno",direction = -1) +
    labs(x="Modelled fAPAR",y = "MODIS fAPAR") +
    geom_abline(slope=1, intercept = 0,col='darkgrey', linetype = "dashed")+
    theme(axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 10))+
    scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
    theme_light()+
    annotate('text', x = 0.15, y = 0.9, label = paste("RMSE = ",round(RMSE(mod1),2),sep = ""), size = 4.5)
}

wrap_plots(p, ncol = 2)+plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a') & 
  scale_color_viridis(option = "inferno",direction = -1, limits = c(0,32000)) 

## supp figure 5 & 6 ####
# region where poly fit were better (from ANOVA test) and significance of poly fit in these region
load("./Data/polyfit_better_region.rda")

slp_q = t(qp_sig)
longDataQ<-melt(slp_q)
res = 1.5
longDataQ$Var1 <- rep(seq(-180 + res/2, 180, res),120)
longDataQ$Var2 <- rep(seq(-90 + res/2, 90, res),each = 240)
longDataQ <- longDataQ[longDataQ$value!=0,]

mp + geom_tile(data = longDataQ, aes(x = Var1, y = -Var2,fill=value),
                      inherit.aes = T, show.legend = T) + 
  scale_fill_stepsn(limits=c(0,1), 
                    labels=NULL,
                    breaks = c(0,0.05,1),
                    colors = brewer.pal(8,"BrBG"),
                    show.limits=T) +
  labs(x="Longitude", y="Latitude", fill='trend') +
  labs(title = "Significance of quadratic polynomial regression") +
  scale_y_continuous(breaks=seq(90,-90,-30))+
  scale_x_continuous(breaks=seq(-180,180,60)) +
  theme_minimal() + 
  theme(legend.position = "none")

# bent cable analysis
load("./Data/bent_cable_analysis.rda")

# turning point
tp <- bc_trend$tau
# only select middle 70% (2003-2014)
tp[tp<2003] <- NA
tp[tp>2014] <- NA
# generate tp mask for selected period
tp_mask <- tp
tp_mask[is.na(tp_mask)==F] <- 1
# break years into specific time slots
tp_period <- tp
tp_period[tp_period>=2003 & tp_period<2006] <- 1
tp_period[tp_period>=2006 & tp_period<2009] <- 2
tp_period[tp_period>=2009 & tp_period<2012] <- 3
tp_period[tp_period>=2012 & tp_period<=2014] <- 4

# trend shape in selected period
status_lim <- bc_trend$status * tp_mask

# show in bar chart 
df <- data.frame(reponse = as.vector(status_lim),
                 tau = as.vector(tp_period))
table(df$reponse,df$tau)
df_count <- data.frame(response = rep(c("Continued increase","Hat shaped", "Cup shaped", "Continued decrease"), each=4),
                       tau = rep(c("2003-2005","2006-2008","2009-2011","2012-2014"),4),
                       count = c(t(table(df$reponse,df$tau))))
tp1 = ggplot(data=df_count, aes(x=response, y=count, fill=tau)) +
  geom_bar(stat="identity", position=position_dodge())+
  xlab("Type of trend shift")+
  ylab("Number of gridcells")+
  labs(fill="Turning point")+
  scale_fill_jco()+
  theme_minimal()
tp1

# only separate into before 2008 and after 2009 (early/late period)
tp[tp<2009] <- -1
tp[tp>=2009] <- 1
status_lim_tau <- status_lim * tp

slp_t = t(status_lim_tau)
longData<-melt(slp_t)
res = 1.5
longData$Var1 <- rep(seq(-180 + res/2, 180, res),120)
longData$Var2 <- rep(seq(-90 + res/2, 90, res),each = 240)
longData <- longData[longData$value!=0,]
mypal = brewer.pal(n = 8, name = "RdBu")

tp2 <- mp + geom_tile(data = longData, aes(x = Var1, y = -Var2,fill=factor(value)),
                      inherit.aes = T, show.legend = T) + 
  scale_fill_manual(name = "Type of shift and turning point",
                    values =c("-4"=mypal[1],"-3"=mypal[2], "-2"=mypal[3], "-1"=mypal[4],
                              "1"=mypal[5], "2"=mypal[6], "3"=mypal[7],"4"=mypal[8]),
                    breaks = c("-4","-3","-2","-1","1","2","3","4"),
                    labels=c("Continued decrease (E)","Cup shaped (E)", "Hat shaped (E)", "Continued increase (E)",
                             "Continued increase (L)","Hat shaped (L)", "Cup shaped (L)", "Continued decrease (L)")) +
  labs(x="Longitude", y="Latitude", fill='trend') +
  scale_y_continuous(breaks=seq(90,-90,-30))+
  scale_x_continuous(breaks=seq(-180,180,60)) +
  theme_minimal() + 
  theme(legend.position = "right", legend.key.width= unit(0.4, 'cm')) 
tp2

tp1/tp2 + plot_layout(heights = c(1,1.5)) +plot_annotation(tag_levels = 'a')

## supp figure 8 ####
val <- read.csv("./Data/mean_gstemp_f0nz.csv")

mypal = pal_uchicago()(5)
show_col(mypal)
cols <- c('z'= mypal[4], 'f0' = mypal[5])

ggplot(val, aes(x=mean_temp)) +
  geom_line(aes(y=z, color='z'), size = 1.5) + 
  geom_line(aes(y=f0*25, color = 'f0'), size = 1.5) + 
  scale_color_manual(values = cols,
                     labels=c(expression(f[0]),"z")) +
  labs(color = NULL) +
  xlab("Mean growing season temperature [°C]")+
  scale_x_continuous(limits=c(2.5,28.5), breaks=seq(2.5,28.5,2)) +
  scale_y_continuous(
    limits = c(0,30),breaks = c(seq(0,30,5)),
    # Features of the first axis
    name = "z",
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans = ~./25,name=expression(f[0]))
  ) +
  theme_bw()+
  theme(legend.position = c(0.1,0.9),legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 15),axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10))



