   library(survival)
   library(survminer)
   library(magrittr)
   library(reshape2)
   library(dplyr)
   library(tidyr)
   library(ggplot2)
   library(broom)
   library(GGally)
   library(compositions) #clr
   library(patchwork)


   #folder_path = ""

   exposure_all <- read.csv(file=paste0(folder_path,"ESigs_integrated_measure_4146_20210501.csv")) 
   types=unique(exposure_all$cancertype)
   
   
   ### PFI
   for (i in types) {
     print(i)
     exposure <- exposure_all %>%
       filter(PFI.time<1200,PFI.time >0, !is.na(PFI.x), cancertype==i) %>%
       mutate(sig1_cat=ifelse(sig_1<= quantile(sig_1,prob=c(0.33,0.66))[1],"low",ifelse(sig_1>=quantile(sig_1,prob=c(0.33,0.66))[2],"high","medium")),
              age.cat = car::recode(as.numeric(age_at_initial_pathologic_diagnosis), "lo:40=1; 40:50=2; 50:60=3; 60:70=4; 70:80=5; 80:hi=6")) %>%
       mutate(sig1_cat=factor(sig1_cat,levels=c("low","medium","high")),
              age.cat= factor(age.cat),
              subtype=factor(subtype,levels=c("ESig1","EHybrid","ESig3","ESig4"))) 
      
              

       exposure[,1:4] = apply(exposure[,1:4],2,clr)
     
       mod1 =coxph(Surv(PFI.time,PFI.x==1) ~ sig_1+sig_2+sig_3+sig_4,data=exposure,robust=TRUE)
       
       mod2 =coxph(Surv(PFI.time,PFI.x==1) ~ sig1_cat,data=exposure,robust=TRUE)
       
       fit = survfit(Surv(PFI.time,PFI.x==1) ~ sig1_cat,data=exposure,robust=TRUE)
       km = ggsurvplot(fit,
                  ggtheme = theme_pubr()+ theme(plot.title = element_text(hjust = 0.5, face = "bold"))+ # theme_minimal will give a white background with grid lines on the plot
                    theme(plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic")),
                  data=exposure,
                  xlab = "Days", ylab = "Survival probability",
                  title = paste0(i," - PFI"),
                  font.x=c(18,"bold"), # changes x axis labels
                  font.y=c(18,"bold"), # changes y axis labels
                  font.xtickslab=c(14,"plain"), # changes the tick label on x axis
                  font.ytickslab=c(14,"plain"),
                  surv.median.line = "hv",
                  ######## Format Legend #######
                  legend = c(0.15,0.3),
                  legend.title = "Sig1(clr)",
                  legend.labs = c("Sig1_low","Sig1_medium","Sig1_high"),
                  pval = TRUE,
                  pval.size = 4,
                  pval.coord = c(1,0),
                  censor.shape="|",
                  censor.size = 4)
       
       
       p1=ggforest(mod1,data=exposure)
       p2=ggforest(mod2,data=exposure)
       
      
       pdf(paste0(i,"_Survival_PFI_sig1.pdf"),height=6,width=15)
       print(km$plot|(p1/p2))
       dev.off()
   }
   
   ## os
   for (i in types) {
     #i="COADREAD"
     print(i)
     exposure <- exposure_all %>%
       filter(OS.time<1200,OS.time >0, !is.na(OS.x), cancertype==i) %>%
       mutate(sig1_cat=ifelse(sig_1<= quantile(sig_1,prob=c(0.33,0.66))[1],"low",ifelse(sig_1>=quantile(sig_1,prob=c(0.33,0.66))[2],"high","medium")),
              age.cat = car::recode(as.numeric(age_at_initial_pathologic_diagnosis), "lo:40=1; 40:50=2; 50:60=3; 60:70=4; 70:80=5; 80:hi=6")) %>%
       mutate(sig1_cat=factor(sig1_cat,levels=c("low","medium","high")),
              age.cat= factor(age.cat),
              subtype=factor(subtype,levels=c("ESig1","EHybrid","ESig3","ESig4"))) 
     
     
     
     exposure[,1:4] = apply(exposure[,1:4],2,clr)
     
     mod1 =coxph(Surv(OS.time,OS.x==1) ~ sig_1+sig_2+sig_3+sig_4,data=exposure,robust=TRUE)
     
     mod2 =coxph(Surv(OS.time,OS.x==1) ~ sig1_cat,data=exposure,robust=TRUE)
     
     fit = survfit(Surv(OS.time,OS.x==1) ~ sig1_cat,data=exposure,robust=TRUE)
     km = ggsurvplot(fit,
                     ggtheme = theme_pubr()+ theme(plot.title = element_text(hjust = 0.5, face = "bold"))+ # theme_minimal will give a white background with grid lines on the plot
                       theme(plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic")),
                     data=exposure,
                     xlab = "Days", ylab = "Survival probability",
                     title = paste0(i," - OS"),
                     font.x=c(18,"bold"), # changes x axis labels
                     font.y=c(18,"bold"), # changes y axis labels
                     font.xtickslab=c(14,"plain"), # changes the tick label on x axis
                     font.ytickslab=c(14,"plain"),
                     surv.median.line = "hv",
                     ######## Format Legend #######
                     legend = c(0.15,0.3),
                     legend.title = "Sig1(clr)",
                     legend.labs = c("Sig1_low","Sig1_medium","Sig1_high"),
                     pval = TRUE,
                     pval.size = 4,
                     pval.coord = c(1,0),
                     censor.shape="|",
                     censor.size = 4)
     
     
     p1=ggforest(mod1,data=exposure)
     p2=ggforest(mod2,data=exposure)
     
     
     pdf(paste0(i,"_Survival_OS_sig1.pdf"),height=6,width=15)
     print(km$plot|(p1/p2))
     dev.off()
   }
   
   