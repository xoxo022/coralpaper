##---------------------------------------------
## Explore Density of recruits over time  
##---------------------------------------------

rm(list=ls())


library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(patchwork)## Arrange multiple plots for compound figure
library(ggpubr)
library(here)
library(corrplot)
library(data.table)
library(wesanderson)
library(ggpmisc)
library(mgcv)

myplots=paste(here(), "Plots", sep="/")

recs<-readRDS("juveniles.Rdata")


juv=5 ## number of replicated quadrats
se <- function(x) sd(x) / sqrt(length(x))

seasoncols=rev(wes_palette(17,2, type=c("discrete")))

p3<-recs%>%filter(!is.na(juvperQ))%>%group_by(Season,Site,Quadrat)%>%
  summarise(meanperm2=mean(juvperm2)) %>%
  ggplot(aes(x = Season, y=meanperm2, fill = factor(Season))) +
  geom_boxplot()+
  facet_grid(~factor(Site))+
  ylab(bquote(Juveniles~per~m^2~"(<"~.(juv)~ "cm)"))+
  xlab("")+
  scale_x_discrete() +
  scale_fill_manual(values=seasoncols)+
  theme_bw() +  
  theme(legend.position = "none",
        legend.title =element_blank(), text = element_text(size=7), 
        legend.key.size = unit(0.5,"cm"),
        legend.text = element_text(size=6))

##--------------------------------
## Explore cover data
##--------------------------------

benthic.Quadrats<-readRDS("benthic_data.Rdata")

cols<-wes_palette(9,8, type=c("continuous"))


## Total cover per season
cover.seasons<-benthic.Quadrats%>%
  dplyr::select(Site:Season)%>%
  group_by(Site,Season,group.code) %>%
  summarise(COUNT = mean(COUNT),
            TOTAL = mean(TOTAL)
  ) %>% 
  ungroup() %>% 
  droplevels()%>%arrange(Site,Season)%>%
  filter(!group.code %in% c("MI"))%>%##remove inverts 
  ggplot(aes(x=Season, y=COUNT/TOTAL,fill=group.code))+
  scale_fill_manual(values=cols)+
  geom_bar(position = "fill",stat="identity")+
  ylab("Proportion of benthic substrate")+
  xlab("")+
  facet_wrap(~Site, scales = "free")+
  theme_bw()+theme(legend.position = "top",
                   legend.title =element_blank(), text = element_text(size=7), 
                   legend.key.size = unit(0.3,"cm"),
                   legend.text = element_text(size=4))+
  guides(fill=guide_legend(nrow=1))



## Explore SFD corals at initial survey
source("SizeFreqDist.R")


sitecolors <- viridis::viridis(2, alpha = 0.6, end = 0.7, option="D")

summer.SFD<-mean.sfd.area.total%>%
  ggplot(aes(x = class, y= TotalFreq, fill = Site)) +
  # Layer for Site 1 with narrower bars
  geom_bar(data = subset(mean.sfd.area.total, Site == "Bundegi"), position = "dodge", width = 0.7, alpha = 0.3, stat="identity") +
  # Layer for Site 2 with wider bars
  geom_bar(data = subset(mean.sfd.area.total, Site == "Tantabiddi"), position = "dodge", width = 0.4, alpha = 0.7, stat="identity") +
  ylab(expression("Average colonies (50 m"^2~")"))+
  xlab("Size-class (cm)")+
  scale_fill_viridis_d(alpha=0.6,end=0.7)+ 
  theme_bw()+
  theme(legend.position="inside",legend.position.inside = c(0.2,0.8),
        legend.title =element_blank(), text = element_text(size=7), 
        legend.key.size = unit(0.3,"cm"),
        legend.text = element_text(size=7),
        legend.background = element_rect(fill ="transparent"))


##Plot relative contribution per taxa
taxa.cols<-wes_palette(6,4, type=c("continuous"))
taxa.cols[2]<-"#A94F6F"

taxa.sfd<-mean.sfd.area.taxa%>%
  ggplot(aes(x = class, y= meanFreq, fill = `Coral group`)) +
  geom_bar(position = "stack", width = 0.5, alpha = 0.5, stat="identity") +
  geom_bar(position = "stack", width = 0.5, alpha = 0.5, stat="identity") +
  scale_fill_manual(values=taxa.cols)+
  ylab(expression("Average colonies (50 m"^2~")"))+
  xlab("Size-class (cm)")+
  facet_grid(~Site)+ 
  theme_bw()+
  theme(legend.position="inside",legend.position.inside = c(0.7,0.8),
        legend.title =element_blank(), text = element_text(size=7), 
        legend.key.size = unit(0.3,"cm"),
        legend.text = element_text(size=7),
        legend.background = element_rect(fill ="transparent"))


##---------------------------------------------
## Explore borad-benthic as predictors of Juvenile density
##---------------------------------------------

recs$survey <- format(recs$Date, "%m-%Y")
recs.benthic=recs%>%left_join(benthic.Quadrats)

##Correlate variables per site
library(energy)## calculate distance correlations

cors<-recs.benthic%>%
  filter(!is.na(juvperQ))%>%
  filter(group.code %in% c("HC", "MA","SS","TA"))%>%
  mutate(cover=(COUNT/TOTAL)*100)%>%
  dplyr::select(Site,Quadrat,juvperm2,group.code,cover)## density standardised to m2



## Plot Associations per cover group- based on gam models with random effects for Quadrat/site

cors$Site=as.factor(cors$Site)
cors$Quadrat=as.factor(cors$Quadrat)

hc.mod=cors%>%filter(group.code=="HC")
hcmod<-gam(juvperm2 ~s(sqrt(cover),k=3)+
             s(Site,bs='re')+## Random intercept 
             s(0+Quadrat,Site,bs='re'), ## Random slope of plot within site
           family=tw(link="log"),data=hc.mod,method = "REML")
summary(hcmod)

##predict on new data
hc_data <- expand.grid(cover = seq(min(hc.mod$cover), max(hc.mod$cover), 
                                   length.out = 100),
                       Site= hc.mod$Site[1],
                       Quadrat= hc.mod$Quadrat[1])

hcpreds<-predict.gam(hcmod,newdata=hc_data,type="response",se.fit=TRUE,exclude= c('s(Site)','s(Quadrat,Site)'))

hc = hc_data %>%data.frame(hcpreds)%>%
  group_by(cover)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  mutate(benthic="HC")%>%
  ungroup()

##Macro
ma.mod=cors%>%filter(group.code=="MA")

mamod<-gam(juvperm2 ~s(sqrt(cover),k=3)+
             s(Site,bs='re')+## Random intercept 
             s(0+Quadrat,Site,bs='re'), ## Random slope of plot within site
           family=tw(link="log"),data=ma.mod,method = "REML")
summary(mamod)

ma_data <- expand.grid(cover = seq(min(ma.mod$cover), max(ma.mod$cover), 
                                   length.out = 100),
                       Site= ma.mod$Site[1],
                       Quadrat= ma.mod$Quadrat[1])
mapreds<-predict(mamod,newdata=ma_data,type="response",se.fit=TRUE,exclude= c('s(Site)','s(Quadrat,Site)'))

ma = ma_data %>%data.frame(mapreds)%>%
  group_by(cover)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  mutate(benthic="MA")%>%
  ungroup()


ss.mod=cors%>%filter(group.code=="SS")
ssmod<-gam(juvperm2 ~s(sqrt(cover),k=3)+
             s(Site,bs='re')+## Random intercept 
             s(0+Quadrat,Site,bs='re'), ## Random slope of plot within site
           family=tw(link="log"),data=ss.mod,method = "REML")
summary(ssmod)

ss_data <- expand.grid(cover = seq(min(ss.mod$cover), max(ss.mod$cover), 
                                   length.out = 100),
                       Site= ss.mod$Site[1],
                       Quadrat= ss.mod$Quadrat[1])
sspreds<-predict(ssmod,newdata=ss_data,type="response",se.fit=TRUE,exclude= c('s(Site)','s(Quadrat,Site)'))

ss = ss_data %>%data.frame(sspreds)%>%
  group_by(cover)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  mutate(benthic="SS")%>%
  ungroup()


ta.mod=cors%>%filter(group.code=="TA")
tamod<-gam(juvperm2 ~s(sqrt(cover),k=3)+
             s(Site,bs='re')+## Random intercept 
             s(0+Quadrat,Site,bs='re'), ## Random slope of plot within site
           family=tw(link="log"),data=ta.mod,method = "REML")
summary(tamod)

ta_data <- expand.grid(cover = seq(min(ta.mod$cover), max(ta.mod$cover), 
                                   length.out = 100),
                       Site= ta.mod$Site[1],
                       Quadrat= ta.mod$Quadrat[1])
tapreds<-predict(tamod,newdata=ta_data,type="response",se.fit=TRUE,exclude= c('s(Site)','s(Quadrat,Site)'))

ta = ta_data %>%data.frame(tapreds)%>%
  group_by(cover)%>% 
  summarise(response=mean(fit), lower=mean(fit-(1.96*se.fit)),
            upper=mean(fit+(1.96*se.fit)))%>%
  mutate(benthic="TA")%>%
  ungroup()


##Plot each
hr2<-paste0("R^2==", round(1-(hcmod$deviance/hcmod$null.deviance),2))

HC.cor<-
  ggplot()+
  geom_point(data=subset(cors,group.code=="HC"), aes(x=cover, y=juvperm2, color=Site))+
  geom_line(data = hc, aes(x= cover, y=response), color = "grey", size = 1)  +
  geom_ribbon(data = hc, aes(x=cover, y=response,ymin= lower, ymax = upper), alpha = 0.3) +
  scale_color_viridis_d(alpha=0.6,end=0.7)+
  facet_wrap(~benthic, scales="free_x",
             labeller = labeller(benthic = c("HC" = "Hard coral (HC)"))) +
  ylab(bquote(Juveniles~per~m^2~"(<"~.(juv)~ "cm)"))+
  xlab("% cover")+
  theme_bw() +  
  theme(legend.position="none",
        text = element_text(size=7))+
  geom_text(data = data.frame(
    x = c(50,50),  # x positions for annotations
    y = c(1.5,0.5),   # y positions for annotations
    label = c("Bundegi", "Tantabiddi"),  # Site labels
    Site = c("Bundegi", "Tantabiddi"),  # Site names
    group.code = "HC"  # Only one facet
  ), aes(x = x, y = y, label = label, color = Site), 
  size = 3)+
  annotate(geom="text", y=15, x=30, label= hr2, parse =T)


mr2<-paste0("R^2==", round(1-(mamod$deviance/mamod$null.deviance),2))
MA.cor<- ggplot()+
  geom_point(data=subset(cors,group.code=="MA"), aes(x=cover, y=juvperm2, color=Site))+
  geom_line(data = ma, aes(x= cover, y=response), color = "grey", size = 1)  +
  geom_ribbon(data = ma, aes(x=cover, y=response,ymin= lower, ymax = upper), alpha = 0.3) +
  scale_color_viridis_d(alpha=0.6,end=0.7)+
  facet_wrap(~group.code,
             labeller = labeller(group.code = c("MA" = "Macroalgae (MA)"))) +
  ylab(bquote(Juveniles~per~m^2~"(<"~.(juv)~ "cm)"))+
  xlab("% cover")+ 
  theme_bw() +  
  theme(legend.position="none",
        text = element_text(size=7))+
  annotate(geom="text", y=15, x=10, label= mr2, parse =T)

tr2<-paste0("R^2==", round(1-(tamod$deviance/tamod$null.deviance),2))
TA.cor<-ggplot()+
  geom_point(data=subset(cors,group.code=="TA"), aes(x=cover, y=juvperm2, color=Site))+
  geom_line(data = ta, aes(x= cover, y=response), color = "grey", size = 1)  +
  geom_ribbon(data = ta, aes(x=cover, y=response,ymin= lower, ymax = upper), alpha = 0.3) +  facet_wrap(~group.code, scales="free_x",
                                                                                                        labeller = labeller(group.code = c("TA" = "Turf (TA)"))) +
  ylab(bquote(Juveniles~per~m^2~"(<"~.(juv)~ "cm)"))+
  xlab("% cover")+
  scale_color_viridis_d(alpha=0.6,end=0.7)+ 
  theme_bw() +  
  theme(legend.position="none",
        text = element_text(size=7))+
  annotate(geom="text", y=15, x=35, label= tr2, parse =T)


sr2<-paste0("R^2==", round(1-(ssmod$deviance/ssmod$null.deviance),2))
SS.cor<-ggplot()+
  geom_point(data=subset(cors,group.code=="SS"), aes(x=cover, y=juvperm2, color=Site))+
  geom_line(data = ss, aes(x= cover, y=response), color = "grey", size = 1)  +
  geom_ribbon(data = ss, aes(x=cover, y=response,ymin= lower, ymax = upper), alpha = 0.3) +  facet_wrap(~group.code, scales="free_x",
                                                                                                        labeller = labeller(group.code = c("SS" = "Sediments (SS)"))) +
  ylab(bquote(Juveniles~per~m^2~"(<"~.(juv)~ "cm)"))+
  xlab("% cover")+
  scale_color_viridis_d(alpha=0.6,end=0.7)+ 
  theme_bw() +  
  theme(legend.position="none",
        text = element_text(size=7))+
  annotate(geom="text", y=15, x=10, label= sr2, parse =T)


cors.plot=ggarrange(HC.cor,
                    MA.cor+ylab("")+theme(axis.text.y = element_blank(), axis.ticks.y =element_blank()),
                    TA.cor+ylab(""),
                    SS.cor+ylab("")+theme(axis.text.y = element_blank(), axis.ticks.y =element_blank()), 
                    align="hv",
                    ncol=2,
                    nrow=2)


##Load Map
install.packages("imager")
library(imager)
ningaloo <- load.image("Sites.png")
plot(ningaloo,axes = FALSE)
map_df <- as.data.frame(ningaloo, wide = "c")

map <- ggplot(map_df, aes(x, y)) +
  geom_raster(aes(fill = rgb(c.1, c.2, c.3))) +
  scale_fill_identity() +
  scale_y_reverse() +  # Reverse y-axis to match image orientation
  theme_void() 


Fig1<-map/((summer.SFD+p3)/
             (taxa.sfd+cover.seasons))+ plot_annotation(tag_levels = list(c('a','b', 'd', 'c','e')))


ggsave(Fig1, filename =paste(myplots, "/Cover-Recruits.png",sep = "/"), dpi=300, width = 6, height = 9)
ggsave(cors.plot, filename =paste(myplots, "/SuppFig_Cover-Recruits.png",sep = "/"), dpi=300, width = 6, height = 6)


