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
dat$scfEnergy <-as.numeric(as.character(dat$scfEnergy))
dat$spinContam <-as.numeric(as.character(dat$spinContam))
dat$corEnergy <-as.numeric(as.character(dat$corEnergy))
dat$Ecorr <- dat$ccsdEnergy-  dat$scfEnergy 
dat$SCFact <- dat$ccsdEnergy-  dat$corEnergy 
dat$metal <- factor(dat$metal)
dat$metal <- revalue(dat$metal,c("2"="Fe","3"="Co"))

## unit converstion 
HF_to_kcalmol= 627.5095 


corr_CBS_extrap <- function(Ecorr2,Ecorr3){
  ## calc CBS based on DZ/TZ extrap
  ## based on
  ## https://pubs.acs.org/doi/abs/10.1021/ct100396y
  beta <- 2.46
  corr_inf <- ((2^beta)*Ecorr2 + (3^beta)*Ecorr3)/((2^beta + 3^beta))
  return(corr_inf)
}
SCF_CBS_extrap <- function(Escf2,Escf3){
  ## calc CBS based on DZ/TZ extrap
  ## based on
  ## https://aip.scitation.org/doi/abs/10.1063/1.3009651
  ## and
  ## https://pubs.acs.org/doi/abs/10.1021/ct100396y
  alpha <-  4.42
  deltaSCF <- Escf2 - Escf3
  A <- deltaSCF/(exp(-alpha*sqrt(2))-exp(-alpha*sqrt(3)))
  scf_inf <- Escf3 - A*exp(-alpha*sqrt(3))
  return(scf_inf)
}





## run times
{
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
}

# check T1 amplitudes
{
  graphics.off()
  x11()
  dat$lig<-factor(dat$lig)
  g <- ggplot(data=dat, aes(x=lig,y=T1,fill=basis,group=basis)) +
    geom_col(position=position_dodge())+
    theme_normal() + xlab('ligand fields') + ylab('T1 diagnostic')+
    scale_fill_manual(values=c("cc-pVDZ"='gray46',
                               "cc-pVTZ"='black'))+
    geom_hline(yintercept = 0.02,color='red',size=1.5)+
    theme(legend.position = c(0.85,0.05),legend.text = element_text(size=9),axis.text.x = element_blank()) +
    facet_grid(spin ~ .) 
  print(g)
  ## save this extra tall 
  cairo_pdf(file="T1VsBasis.pdf",width=1.25*3.33,height=2.5*2.95)
  print(g)
  dev.off() 
}

# check T1 vs spin error 
{
  graphics.off()
  x11()
  dat$lig<-factor(dat$lig)
  g <- ggplot(data=dat, aes(y=T1,x=spinContam,fill=basis,group=basis)) +
    geom_point(shape=21,size=4)+
    theme_normal() + xlab('spin contamination') + ylab('T1 diagnostic')+
    scale_fill_manual(values=c("cc-pVDZ"='gray46',
                               "cc-pVTZ"='black'))+
    geom_hline(yintercept = 0.02,color='red',size=1.5)+
    theme(legend.position = c(0.85,0.15),legend.text = element_text(size=9))# +
    #facet_grid(spin ~ .) 
  print(g)
  ## save this extra tall 
  cairo_pdf(file="T1VsSC.pdf",width=1.25*3.33,height=2.5*2.95)
  print(g)
  dev.off() 
}

## CBS extrapolation
{
# let's do a cast to put the different basis set energies in the same rows (this will delete times, T1 etc)
require(data.table)
dat.c <- dcast.data.table(as.data.table(dat),  slot + metal + ox + lig + spin  ~ basis,value.var = c('scfEnergy','SCFact','corEnergy','ccsdEnergy','T1'))
dat.c[,"corEnergy-CBS"] <- apply(dat.c[,c("corEnergy_cc-pVDZ","corEnergy_cc-pVTZ")],1,function(x) corr_CBS_extrap(x[["corEnergy_cc-pVDZ"]],x[["corEnergy_cc-pVTZ"]]))
dat.c[,"SCFact-CBS"] <- apply(dat.c[,c("SCFact_cc-pVDZ","SCFact_cc-pVTZ")],1,function(x) SCF_CBS_extrap(x[["SCFact_cc-pVDZ"]],x[["SCFact_cc-pVTZ"]]))
dat.c[,'cbsEnergy']<-dat.c[,"SCFact-CBS"]+dat.c[,"corEnergy-CBS"]
}


## correaltion energies 
dat.m <- melt(dat.c,measure.vars = c('corEnergy-CBS','corEnergy_cc-pVTZ','corEnergy_cc-pVDZ'))
dat.m$basis <- revalue(dat.m$variable,c('corEnergy-CBS'='CBS','corEnergy_cc-pVTZ'='cc-pVTZ','corEnergy_cc-pVDZ'='cc-pVDZ'))
dat.m$corrE <- dat.m$value
dat.m$value<-NULL
dat.m$variable<-NULL
dat.m$lig<-factor(dat.m$lig)
dat.m=dat.m[order(dat.m$corrE),]

{
  graphics.off()
  g <- ggplot(data=dat.m[dat.m$spin<=2], aes(x=lig,y=corrE*HF_to_kcalmol,fill=basis,group=basis)) +
    geom_col(position=position_dodge())+
    theme_normal() + xlab('ligand fields') + ylab('correlation energy, [kcal/mol]')+
    scale_fill_manual(values=c("cc-pVDZ"='gray46',
                               "cc-pVTZ"='black',
                               "CBS" = 'firebrick'))+
    theme(legend.position = c(0.205,0.90),
          legend.text = element_text(size=8),axis.text.x = element_blank(),
          legend.key.size = unit(x = 0.25,units = "cm"),
          legend.background = element_rect(fill='white'),
          legend.margin = margin(0.05,0.05,0.05,0.05,unit='cm')
          ) +
    facet_grid(spin ~ .) 
  print(g)
  cairo_pdf(file="corrVsBasis-CBS-low-spin.pdf",width=3.33,height=2.95)
  print(g)
  dev.off() 
  ## ## ## 
  graphics.off()
  g <- ggplot(data=dat.m[dat.m$spin>2], aes(x=lig,y=corrE*HF_to_kcalmol,fill=basis,group=basis)) +
    geom_col(position=position_dodge())+
    theme_normal() + xlab('ligand fields') + ylab('correlation energy, [kcal/mol]')+
    scale_fill_manual(values=c("cc-pVDZ"='gray46',
                               "cc-pVTZ"='black',
                               "CBS" = 'firebrick'))+
    theme(legend.position = c(0.205,0.90),
          legend.text = element_text(size=8),axis.text.x = element_blank(),
          legend.key.size = unit(x = 0.25,units = "cm"),
          legend.background = element_rect(fill='white'),
          legend.margin = margin(0.05,0.05,0.05,0.05,unit='cm')
    ) +
    facet_grid(spin ~ .) 
  print(g)
  cairo_pdf(file="corrVsBasis-CBS-high-spin.pdf",width=3.33,height=2.95)
  print(g)
  dev.off() 
}


## total e correction
dat.m <- melt(dat.c,measure.vars = c('cbsEnergy','ccsdEnergy_cc-pVTZ','ccsdEnergy_cc-pVDZ'))
dat.m$basis <- revalue(dat.m$variable,c('cbsEnergy'='CBS','ccsdEnergy_cc-pVTZ'='cc-pVTZ','ccsdEnergy_cc-pVDZ'='cc-pVDZ'))
dat.m$totalE <- dat.m$value
dat.m$value<-NULL
dat.m$variable<-NULL
dat.m$lig<-factor(dat.m$lig)

## assemble rows by spliting energy
findpartners <- function(row,df){
   this_metal <- row[["metal"]]
   this_ox <- row[["ox"]]
   this_lig <- row[["lig"]]
   this_basis <- row[["basis"]]
   this_spin <- row[["spin"]]
   this_energy <- row[["totalE"]]
   print(c(this_metal,this_ox,this_lig,this_basis))
   partner = df[df$metal == this_metal &
                 df$ox == this_ox &
                 df$lig == this_lig &
                 df$spin != this_spin &
                 as.character(df$basis) == as.character(this_basis), ]
   print(partner)
   if (partner$spin >  this_spin)
     {
      this_split <- HF_to_kcalmol*(partner$totalE - this_energy)
      this_spin_cat <- "LS"
     }
   else{
     this_split <- HF_to_kcalmol*(-partner$totalE + this_energy)
     this_spin_cat <- "HS"
   }
   rl  = list()
   rl["split"] = this_split
   rl["sc"] = this_spin_cat
   return(rl)
}
dat.m$spin_cat <- "undef"
dat.m$split <- "undef"
for (i in seq(1,nrow(dat.m))){
  rl = findpartners(row=dat.m[i,],df=dat.m)
  dat.m[i,]$spin_cat <- rl$sc
  dat.m[i,]$split <- rl$split
}

## plot basis set dependence of spin splitting
require(data.table)
dat.c <- dcast.data.table(as.data.table(dat.m),
                          slot + metal + ox + lig + spin + spin_cat  ~ basis,
                          value.var = c('split'))

dat.c <- dat.c[!is.na(dat.c$`cc-pVDZ`),]
dat.c<-as.data.frame(dat.c)
dat.c$`cc-pVDZ` <-as.numeric(as.character(dat.c$`cc-pVDZ`))
dat.c$`cc-pVTZ` <-as.numeric(as.character(dat.c$`cc-pVTZ`))
dat.c$CBS <-as.numeric(as.character(dat.c$CBS))
dat.c$ox<-factor(dat.c$ox)

dat.cm <- melt(dat.c,measure.vars=c('cc-pVTZ','cc-pVDZ','CBS'))
dat.cm$value <- as.numeric(as.character(dat.cm$value ))
dat.cm$lig <- factor(dat.cm$lig)
dat.cm$slot <- factor(dat.cm$slot)
dat.cm<-dat.cm[order(dat.cm$ox,dat.cm$value),]
{
  graphics.off()
  x11()

  scay = seq(from=-15,to=55,by=10)
  g <- ggplot(data=dat.c[dat.c$spin_cat == "LS", ], aes(y=`cc-pVTZ`,x=`cc-pVDZ`,
                                                        color=metal,shape=ox)) +
    geom_abline(size=1.5,color='gray',lty='dashed')
  g <- g+ geom_point(size=4)+
    theme_normal() + xlab(expression(paste(Delta*E[H-L],', cc-pVDZ [kcal/mol]' ))) + coord_fixed()+
    ylab(expression(paste(Delta*E[H-L],', cc-pVTZ [kcal/mol]' ))) +
    scale_color_manual(values=c("Fe"="orange",
                               "Co"='#0047AB'))+
    theme(legend.position = c(0.85,0.40),legend.text = element_text(size=9))+
    scale_y_continuous(sec.axis = dup_axis(),breaks = scay,labels = scay)+
    scale_x_continuous(sec.axis = dup_axis(),breaks = scay,labels = scay)
    
  print(g)

  cairo_pdf(file="splitTVsDZ.pdf",width=3.33,height=2.95)
  print(g)
  dev.off() 
  
  g <- ggplot(data=dat.c[dat.c$spin_cat == "LS", ], aes(y=`CBS`,x=`cc-pVDZ`,
                                                        color=metal,shape=ox)) +
    geom_abline(size=1.5,color='gray',lty='dashed')
  g <- g+ geom_point(size=4)+
    theme_normal() + xlab(expression(paste(Delta*E[H-L],', cc-pVDZ [kcal/mol]' ))) + coord_fixed()+
    ylab(expression(paste(Delta*E[H-L],', CBS [kcal/mol]' ))) +
    scale_color_manual(values=c("Fe"="orange",
                                "Co"='#0047AB'))+
    theme(legend.position = c(0.85,0.40),legend.text = element_text(size=9))+
    scale_y_continuous(sec.axis = dup_axis(),breaks = scay,labels = scay)+
    scale_x_continuous(sec.axis = dup_axis(),breaks = scay,labels = scay)
  
  print(g)
  
  cairo_pdf(file="splitCBSsDZ.pdf",width=3.33,height=2.95)
  print(g)
  dev.off() 
  
  
  
  g <- ggplot(data=dat.cm[dat.cm$spin_cat == "LS", ],aes(x=slot,y=value,fill=variable)) +
       geom_col(position = position_dodge())+ theme_normal()+xlab('complex ID') + 
       ylab(expression(paste(Delta*E[H-L],', [kcal/mol]' ))) +
       scale_fill_manual(values=c("cc-pVDZ"='gray46',
                               "cc-pVTZ"='black',
                               "CBS" = 'firebrick')) +
       theme(legend.position = c(0.225,0.80))

  print(g)
  cairo_pdf(file="splitAllcomp.pdf",width=3.33,height=2.95)
  print(g)
  dev.off() 
}

## compare DFT
dat.c$DFTsplit = NA

load(file="srm.Rdata")

dim(srm)

sr.m$ox <- sr.m$variable
sr.m$ox<- revalue(sr.m$ox,c("ox_2_split"="2","ox_3_split"="3"))

for (j in seq(1,nrow(sr.m))){


  ll<-as.character(strsplit(as.character(sr.m$gene),split = "_")[[j]][2])
  this_lig<-gsub("smi","",  ll)
  this_ox = sr.m[j,"ox"]
  this_split = sr.m[j,"value"]*HF_to_kcalmol
  
  dat.c[dat.c$ox == this_ox & dat.c$lig == this_lig &
        dat.c$spin_cat == "LS",]$DFTsplit <- this_split
}


{
  
  graphics.off()
  x11()
  
  scay = seq(from=-15,to=55,by=10)
  g <- ggplot(data=dat.c[dat.c$spin_cat == "LS", ], aes(y=CBS,x=DFTsplit,
                                                        color=metal,shape=ox)) +
    geom_abline(size=1.5,color='gray',lty='dashed')
  g <- g+ geom_point(size=4)+
    theme_normal() + xlab(expression(paste(Delta*E[H-L],', cc-pVDZ [kcal/mol]' ))) + #coord_fixed()+
    ylab(expression(paste(Delta*E[H-L],', cc-pVTZ [kcal/mol]' ))) +
    scale_color_manual(values=c("Fe"="orange",
                                "Co"='#0047AB'))+
    theme(legend.position = c(0.85,0.40),legend.text = element_text(size=9))+
    scale_y_continuous(sec.axis = dup_axis(),breaks = scay,labels = scay)+
    scale_x_continuous(sec.axis = dup_axis(),breaks = scay,labels = scay)
  
  print(g)
  
  cairo_pdf(file="splitCBSsDFT.pdf",width=3.33,height=2.95)
  print(g)
  dev.off() 
}

