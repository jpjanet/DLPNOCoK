require(plyr)
require(dplyr)
require(ggplot2)
require(reshape2)
source("plotPrefs.R")
dat<-read.csv(file="dlpnores.csv",sep = ',')
dat$T1 <- as.numeric(as.character(dat$T1))
dat<-dat[!is.na(dat$T1),]
dat$time <- as.numeric(as.character(dat$time))
dat$ccsdEnergy <- as.numeric(as.character(dat$ccsdEnergy))

## unit converstion 
HF_to_kcalmol= 627.5095 

## run times
graphics.off()
x11()
g <- ggplot(data=dat, aes(x=time)) +
  stat_ecdf(aes(color=basis,linetype=basis),size=2.5)+
  theme_normal() + xlab('time, [h]') + ylab('fraction complete')+
  scale_x_continuous(breaks = seq(from=0,by=5,to=25),
                     limits = c(0,25),
                     labels=seq(from=0,by=5,to=25))+
  scale_color_manual(values=c("cc-pVDZ"='gray46',
                       "cc-pVTZ"='black'))+
  theme(legend.position = c(0.75,0.25))
print(g)

cairo_pdf(file="timeVsBasis.pdf",width=3.33,height=2.95)
print(g)
dev.off()
break
## assemble rows by spliting energy
findpartners <- function(row,df){
   this_metal <- row[["metal"]]
   this_ox <- row[["ox"]]
   this_lig <- row[["lig"]]
   this_basis <- row[["basis"]]
   this_spin <- row[["spin"]]
   this_energy <- row[["ccsdEnergy"]]
   print(c(this_metal,this_ox,this_lig,this_basis))
   partner = df[df$metal == this_metal &
                 df$ox == this_ox &
                 df$lig == this_lig &
                 df$spin != this_spin &
                 as.character(df$basis) == as.character(this_basis), ]
   print(partner)
   if (partner$spin >  this_spin)
     {
      this_split <- HF_to_kcalmol*(partner$ccsdEnergy - this_energy)
      this_spin_cat <- "LS"
     }
   else{
     this_split <- HF_to_kcalmol*(-partner$ccsdEnergy + this_energy)
     this_spin_cat <- "HS"
   }
   rl  = list()
   rl["split"] = this_split
   rl["sc"] = this_spin_cat
   return(rl)
}
dat$spin_cat <- "undef"
dat$split <- "undef"
for (i in seq(1,nrow(dat))){
  rl = findpartners(row=dat[i,],df=dat)
  dat[i,]$spin_cat <- rl$sc
  dat[i,]$split <- rl$split
}

dat=dat[order(dat$split),]
