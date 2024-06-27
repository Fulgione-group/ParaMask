#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


library(ggplot2)
library(VGAM)
library(data.table)
library(patchwork)

#default parameters
hetpath <- NULL
outpath <- getwd()
ID<-""
missingness <- 0.1
startline <- 2
endline <- 0
verbose <- FALSE
tolerance<- 0.001
num_iterations<-1000
dist_em_rep<-1000
chr<-0
cdist<- FALSE
min_dist <- 50
max_dist <- 2000
boundfit <- FALSE
useRRD <- TRUE

for(i in 1:length(args)){
	print(args[i])
}

for(i in 1:length(args)){
  if (args[i]=="--het" | args[i]=="-h" ){
     hetpath <- args[(i+1)]
     i <- i + 1
  } else if (args[i]=="--chrom" | args[i]=="-c" ){
    chr <- args[(i+1)]
     i <- i + 1
  } else if (args[i]=="--startline" | args[i]=="-s" ){
    startline <- as.numeric(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--endline" | args[i]=="-e" ){
    endline <- as.numeric(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--missingness" | args[i]=="-m" ){
    missingness <- as.numeric(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--outdir" | args[i]=="-o" ){
    outpath <- as.character(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--ID" | args[i]=="-id" ){
    ID <- as.character(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--verbose" | args[i]=="-v" ){
    verbose <- TRUE
  } else if (args[i]=="--tolerance" | args[i]=="-t" ){
    tolerance <- as.numeric(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--num_it" | args[i]=="-ni" ){
    num_iterations <- as.numeric(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--dist_em_rep" | args[i]=="-dr" ){
    dist_em_rep <- as.numeric(args[(i+1)])
    i <- i + 1
  } else if (args[i]=="--boundary" | args[i]=="-b" ){
    boundary <- as.character(args[(i+1)])
    lboundary <- as.numeric(strsplit(boundary,split = ",")[[1]][1])
    uboundary <- as.numeric(strsplit(boundary,split = ",")[[1]][2])
    boundfit <- TRUE
    i <- i + 1
  } else if (args[i]=="--cdist" | args[i]=="-cd" ){
    cdist <- TRUE
  } else if (args[i]=="--noRRD"){
    useRRD <- FALSE
  }
}

if(!endsWith(x = outpath, suffix = "/")){
  outpath<- paste(outpath, "/", sep="")
}
  

## for local testing
# hetpath<- "/media/btjeng/Elements/ParaMask/ParaMask_bastiaan/SeDuS5_0.2_reps/Simulations_0.2_SV_SC_rep3_mat_cpos.vcf.het.stat.txt"
# outpath<- "/home/btjeng/Data/Paramask/Cluster_test/"
# startline <- 724
# endline <- 100000
# missingness <- 0.1
# iteration<-0
# ID<-"testrun"
################functions
# define function for results

plot_EM_it <- function(het,fit, iteration, outpath,ID){

  predMeans1<-data.frame(MAF=seq((0.5/10000),0.5, by=(0.5/10000)), Z="K")
  predMeans1$MAF<- (1-predMeans1$MAF)
  predMeans2<-data.frame(MAF=rep(0,nrow(predMeans1)), Z=rep("D",nrow(predMeans1)))
  
  
  predMeans1$Mean <- predictvglm(fit3,newdata = predMeans1, type = "response")
  predMeans2$Mean <- predictvglm(fit3,newdata = predMeans2, type = "response")
  
  predMeans1$MAFT <- (1-predMeans1$MAF)
  predMeans2$MAFT <- predMeans1$MAFT
  
  predMeans1$MeanUT<- predMeans1$Mean
  predMeans2$MeanUT<- predMeans2$Mean 
  
  predMeans1$MeanT<- predMeans1$Mean *predMeans1$MAFT*2
  predMeans2$MeanT<- predMeans2$Mean *predMeans2$MAFT*2 
  
  
  predMeans1<-as.data.frame(predMeans1)
  predMeans2<-as.data.frame(predMeans2)
  
  #plotting
  rplot1 <- ggplot(data = het[rs,])+
    geom_point(aes(x = Minor.allele.freq, y = Heterozygous.geno.freq/(2*Minor.allele.freq), color = weight1[rs]))+
    geom_path(data = predMeans1, aes(x = MAFT, y = MeanUT, linetype="Paralogs"), color="black", size=1.2)+
    geom_path(data = predMeans2, aes(x = MAFT, y = MeanUT, linetype="SNPs"), color="black", size=1.2)+
    scale_x_continuous(lim = c(-0.025,0.525), breaks = c(0,0.5), expand=c(0,0))+
    scale_y_continuous(lim = c(-0.05,1.05), breaks = c(0,0.5,1), expand=c(0,0))+
    scale_linetype_manual(name="SNP", labels=c("single-copy SNP","multicopy SNP"),values=c("solid","twodash"))+
    scale_color_gradient(name="SNP",limits=c(0,1),breaks=c(1,0), labels=c("single-copy SNP","multicopy SNP"), "",high="blue3", low="brown2")+
    labs(x= "minor allele frequency (maf)", y = "heterozygote frequency / (2*maf)")+
    coord_fixed(ratio = 0.5)+
    geom_segment(x=-0.025, y=-0.003, xend=-0.025, yend=1.003, size=1.2)+
    geom_segment(x=-0.0013, y=-0.05, xend=0.5013, yend=-0.05, size=1.2)+
    theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=14) , axis.ticks=element_line(size=1),axis.ticks.length = unit(0.25,"cm"), axis.title = element_text(size=16),legend.text = element_text(size=12), legend.position = "right",legend.title=element_blank(), legend.key = element_rect(color="white", fill=NA), legend.key.size = unit(1.5,"cm"), panel.background = element_rect(fill=NA, color="white"))+
    guides(color = guide_legend(override.aes = list(size = 4)))
  
    
  rplot2<-ggplot(data = het[rs,])+
    geom_point(aes(x = Minor.allele.freq, y = Heterozygous.geno.freq, color = weight1[rs]))+
    geom_path(data = predMeans1, aes(x = MAFT, y = MeanT, linetype="Paralogs"), color="black", size=1.2)+
    geom_path(data = predMeans2, aes(x = MAFT, y = MeanT, linetype="SNPs"), color="black", size=1.2)+
    scale_x_continuous(lim = c(-0.025,0.525), breaks = c(0,0.5), expand=c(0,0))+
    scale_y_continuous(lim = c(-0.05,1.05), breaks = c(0,0.5,1), expand=c(0,0))+
    scale_linetype_manual(name="SNP", labels=c("single-copy SNP","multicopy SNP"),values=c("solid","twodash"))+
    scale_color_gradient(name="SNP",limits=c(0,1),breaks=c(1,0), labels=c("single-copy SNP","multicopy SNP"), "",high="blue3", low="brown2")+
    labs(x= "minor allele frequency (maf)", y = "heterozygote frequency")+
    coord_fixed(ratio = 0.5)+
    geom_segment(x=-0.025, y=-0.003, xend=-0.025, yend=1.003, size=1.2)+
    geom_segment(x=-0.0013, y=-0.05, xend=0.5013, yend=-0.05, size=1.2)+
    theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=14) , axis.ticks=element_line(size=1),axis.ticks.length = unit(0.25,"cm"), axis.title = element_text(size=16),legend.text = element_text(size=12), legend.position = c(0.2, 0.6),legend.title=element_blank(), legend.key = element_rect(color="white", fill=NA), legend.key.size = unit(1.5,"cm") ,panel.background = element_rect(fill=NA, color="white"))+
    guides(color = guide_legend(override.aes = list(size = 4)))
  pdf(width = 7, height = 7, file = paste(outpath, ID,"_EM_iteration", as.character(iteration),".pdf", sep = ""))
  print(rplot1)
  print(rplot2)
  dev.off()  
}

##em algorithm for distances 
em_algorithm_geom <- function(data, max_iter = 1000, tol = 1e-6) {
  n <- length(data)
  k <- 2  # Number of components
  # Initialize parameters
  min<- 10^-5
  prob1 <- runif(min =min , max =1-min, n = 1 )
  prob2 <- runif(min =min , max =1-min, n = 1 )
  weights <- runif(min =0.1 , max =0.9, n = n)
  # weights <- rep(0.5, n)
  
  for (iter in 1:max_iter) {
    # E-step
    likelihood1 <- dgeom(data, prob = prob1, log = F)
    likelihood2 <- dgeom(data, prob = prob2, log = F)
    likelihood1[likelihood1==0] <- 10^-300
    likelihood2[likelihood2==0] <- 10^-300
    # Calculate responsibilities
    resp1 <- weights * likelihood1
    resp2 <- (1 - weights) * likelihood2
    total_resp <- resp1 + resp2
    
    weights <- resp1 / total_resp
    
    # M-step
    
    # prob1_new <- sum(weights>=0.5) / sum(data[weights>=0.5])
    # prob2_new <- sum(weights<0.5) / sum(data[weights<0.5]) 
    
    prob1_new <- min(0.9, sum(weights) / sum(data * weights))
    prob2_new <- min(0.9, sum(1-weights) / sum(data * (1-weights)))
    # 
    
    # Check convergence
    if (abs(prob1_new - prob1) < tol && abs(prob2_new - prob2) < tol) {
      break
    }
    
    # Update parameters
    prob1 <- prob1_new
    prob2 <- prob2_new
  }
  
  # Return the estimated parameters
  return(list(prob1 = prob1, prob2 = prob2, weights = weights))
}
###############################

#read hetfile 
if(verbose){
  print("reading hetfile.....")
}
tryCatch({
  if(endline == 0){
    het <- fread(file = hetpath,header = T)
  } else {
    header <- fread(file = hetpath, nrows =1, header = F)
    het <- fread(file = hetpath, skip =(startline-1), nrows = (endline-startline+1),header = F )
    colnames(het) <- as.character(header)
  }
  if(chr!=0){
    het <- het[het$Chromosome==chr,]
  }
}, warning = function(war){
  print(paste("Execution halted, WARNING while loading hetfile:  ",war))
  q(save = "no", status = 1, runLast = F)
  
}, error = function(err){
  print(paste("Execution halted, ERROR while loading hetfile:  ",err))
})


#convert to dataframe
het <- data.frame(het, stringsAsFactors = F)

het <- het[het$Non.missing >= (1-missingness)*max(het$Non.missing) & het$Minor.allele.freq> 0,]

#subsample for plotting
subsamplesize <- min(10000, nrow(het))
rs<-sample(1:nrow(het), size=subsamplesize, replace = F)

#assemble regression data
regData <- cbind((1-het$Minor.allele.freq),  round(het$Heterozygous.geno.freq* het$Non.missing), round(het$Minor.allele.freq*2*het$Non.missing), sample(c(0,1), size = nrow(het),replace = T))
regData <- as.data.frame(regData)  

regData2<- rbind(regData,regData)
regData2$Z<- c(rep("K",nrow(regData)),rep("D",nrow(regData)))
regData2 <- as.data.frame(regData2)  

colnames(regData2) <- c("MAF", "het", "N", "cluster", "Z")
colnames(regData) <- c("MAF", "het", "N", "cluster")
regData2$MAF[regData2$Z=="D"] <- 0 

#initialize SNPs
regData$InCutoff<-apply(X = regData, MARGIN = 1, FUN = function(x){qbinom(p =0.999, size = x[3],prob = x[1], lower.tail = T)})
regData$InCutoff2<-apply(X = regData, MARGIN = 1, FUN = function(x){qbinom(p =0.999, size = x[3],prob = x[1], lower.tail = F)})

predDF<-data.frame(MAF=regData2$MAF, Z=regData2$Z)
weight<-0
weight<- runif(n = nrow(regData), min = 0.01, max = 0.99)
# weight[(regData$het> regData$InCutoff) & (regData$MAF <0.7)] <- 0.01
weight[regData$het> regData$InCutoff] <- 0.01
weight[regData$het==regData$N] <- 0.01
weight[regData$het< regData$InCutoff2] <-0.99

weight1<- weight
weight2<- 1 - weight

if(verbose){
  fit3<-vglm(cbind(regData2$het, regData2$N - regData2$het) ~ Z +  MAF, data=regData2, family = betabinomial(lmu="clogloglink", lrho = "logitlink"),trace=TRUE, subset=(regData2$N > 1), weights = c(weight1,weight2))
}else{
  suppressWarnings(  fit3<-vglm(cbind(regData2$het, regData2$N - regData2$het) ~ Z + MAF, data=regData2, family = betabinomial(lmu="clogloglink", lrho = "logitlink"),trace=FALSE, subset=(regData2$N > 1), weights = c(weight1,weight2))
)
}
plot_EM_it(het = het, fit = fit3, iteration = 0, outpath = outpath, ID = ID)
coef_fit3 <- coefficients(fit3)

###############################################################run EM

offsetfit <- FALSE

if(boundfit){
  if(coef_fit3[4]> uboundary | coef_fit3[4] < lboundary){
    if(coef_fit3[4]> uboundary){
      cboundary <- uboundary
    }else{
      cboundary <- lboundary
    }
    if(verbose){
      print(paste("taking a modified step. Boundary exeeded, slope: ",as.character(coef_fit3[4]), ", with boundary: ", as.character(cboundary),". Using boundary as offset on minor allele frequency to refit...", sep = ""))
      tryCatch({
        fit3<-vglm(cbind(regData2$het, regData2$N - regData2$het) ~ Z + offset(MAF*cboundary), data=regData2, family = betabinomial(lmu="clogloglink", lrho = "logitlink"),start=c(2,1,0,2),trace=TRUE, subset=(regData2$N > 1), weights = c(weight1,weight2))
      }, error = function(err){
        print(paste("ERROR in Mstep during VGAM fitting:  ",err))
        q(save = "no", status = 1, runLast = F)
      })
      coef_fit3<- coefficients(fit3)[1:3]
      coef_fit3[4]<-cboundary
      offsetfit<-TRUE
    }else{
      tryCatch({
        suppressWarnings(fit3<-vglm(cbind(regData2$het, regData2$N - regData2$het) ~ Z + offset(MAF*cboundary), data=regData2, family = betabinomial(lmu="clogloglink", lrho = "logitlink"),start=c(2,1,0,2),trace=TRUE, subset=(regData2$N > 1), weights = c(weight1,weight2)))
    }, error = function(err){
      print(paste("ERROR in Mstep during VGAM fitting:  ",err))
      q(save = "no", status = 1, runLast = F)
    })
      coef_fit3<- coefficients(fit3)[1:3]
      coef(fit3, matrix=T)
      coef_fit3[4]<-cboundary
      offsetfit<-TRUE
    }
  }  
}

pars3<-data.frame(mu=rep(0,nrow(predDF)), rho=rep(0,nrow(predDF)))
log_likelihood <- c()
iteration <-1
for(iteration in 1:num_iterations){
  ##E step to estimate Zj
  if(verbose){
    print(iteration)
    print("Estep: calulating weights.....")
  }
  if(boundfit & offsetfit){
      pars3$mu<- NA
      pars3$rho<- NA
      pars3$mu<-apply(X = predDF, MARGIN = 1, FUN = function(x){z<-ifelse(x[2]=="K", 1,0);lp<-  coef_fit3[1] + z*coef_fit3[3]  + as.numeric(x[1])*z*coef_fit3[4]; return(cloglog(lp, inverse = T))})
      pars3$rho<- logitlink(coef_fit3[2], inverse = T)
      offsetfit <- FALSE
  }else{
    pars3<- predictvglm(fit3,newdata = predDF, type = "link", untransform=TRUE)
    pars3<- as.data.frame(pars3)
    colnames(pars3)<- c("mu", "rho")
  }
  regData2$pred_fit3_mu<- pars3$mu
  regData2$pred_fit3_rho<- pars3$rho
  probs<-NA
  probs<-apply(X = regData2, MARGIN = 1, FUN = function(x){dbetabinom(as.numeric(x[2]), size = as.numeric(x[3]), prob = as.numeric(x[6]), rho  = as.numeric(x[7]))})
  probs1<-probs[regData2$Z=="K"]
  probs2<-probs[regData2$Z=="D"]
  ##store log_likelihood
  # responsibilities<-ifelse(weight1>=0.5, TRUE,FALSE)
  # probs1_ll<-apply(X = regData2[regData2$Z=="K",][responsibilities,], MARGIN = 1, FUN = function(x){dbetabinom(as.numeric(x[2]), size = as.numeric(x[3]), prob = as.numeric(x[6]), rho  = as.numeric(x[7]))})
  # probs2_ll<-apply(X = regData2[regData2$Z=="D",][!responsibilities,], MARGIN = 1, FUN = function(x){dbetabinom(as.numeric(x[2]), size = as.numeric(x[3]), prob = as.numeric(x[6]), rho  = as.numeric(x[7]))})
  # log_likelihood <- rbind(log_likelihood,c(iteration,(sum(log(probs1_ll))+ sum(log(probs2_ll)))))
  log_likelihood <- rbind(log_likelihood,c(iteration,logLik(fit3)))
  ## probs = posterior, weight prior, (weight*probs1 + (1-weight)*probs2) = sum over probabilties
  weight1 <- (weight1*probs1)/(weight1*probs1 + weight2*probs2) 
  weight2 <- (weight2*probs2)/(weight1*probs1 + weight2*probs2) 
  ## make sure weights do not drop to zero or 1
  weight1[weight1<1e-100]<-1e-100
  weight2[weight2<1e-100]<-1e-100
  ##M step
  if(verbose){
    print("Mstep: Model refitting.....")
    tryCatch({
      fit3<-vglm(cbind(regData2$het, regData2$N - regData2$het) ~ Z +  MAF, data=regData2, family = betabinomial(lmu=clogloglink(bvalue = 1e-10),lrho = "logitlink"),trace=TRUE, subset=(regData2$N > 1), weights = c(weight1,weight2))
    }, error = function(err){
      print(paste("ERROR in Mstep during VGAM fitting:  ",err))
      q(save = "no", status = 1, runLast = F)
    })
  }else{
    tryCatch({
      suppressWarnings(fit3<-vglm(cbind(regData2$het, regData2$N - regData2$het) ~ Z +  MAF, data=regData2, family = betabinomial(lmu=clogloglink(bvalue = 1e-10),lrho = "logitlink"),trace=FALSE, subset=(regData2$N > 1), weights = c(weight1,weight2)))
    }, error = function(err){
      print(paste("ERROR in Mstep during VGAM fitting:  ",err))
      q(save = "no", status = 1, runLast = F)
    })
  }
  coef_fit3_new <- coefficients(fit3)
  if(boundfit){
    if(coef_fit3_new[4]> uboundary | coef_fit3_new[4] < lboundary){
      if(coef_fit3_new[4]> uboundary){
        cboundary <- uboundary
      }else{
        cboundary <- lboundary
      }      
      if(verbose){
        print(paste("taking a modified step. Boundary exeeded, slope: ",as.character(coef_fit3_new[4]), ", with boundary: ", as.character(cboundary),". Using boundary as offset on minor allele frequency to refit...", sep = ""))
        tryCatch({
          fit3<-vglm(cbind(regData2$het, regData2$N - regData2$het) ~ Z + offset(MAF*cboundary), data=regData2, family = betabinomial(lmu="clogloglink", lrho = "logitlink"),start=c(2,1,0,2),trace=TRUE, subset=(regData2$N > 1), weights = c(weight1,weight2))
        }, error = function(err){
          print(paste("ERROR in Mstep during VGAM fitting:  ",err))
          q(save = "no", status = 1, runLast = F)
        })
        coef_fit3_new<- coefficients(fit3)[1:3]
        coef_fit3_new[4]<-cboundary
        offsetfit<-TRUE
      }else{
        tryCatch({
          suppressWarnings(fit3<-vglm(cbind(regData2$het, regData2$N - regData2$het) ~ Z + offset(MAF*cboundary), data=regData2, family = betabinomial(lmu="clogloglink", lrho = "logitlink"),start=c(2,1,0,2),trace=TRUE, subset=(regData2$N > 1), weights = c(weight1,weight2)))
        }, error = function(err){
          print(paste("ERROR in Mstep during VGAM fitting:  ",err))
          q(save = "no", status = 1, runLast = F)
        })
        coef_fit3_new<- coefficients(fit3)[1:3]
        #coef(fit3_new, matrix=T)
        coef_fit3_new[4]<-cboundary
        offsetfit<-TRUE
      }
    }  
  }
  if(verbose){
    print("old Parameters")
    print(coef_fit3)
    print("new Parameters")
    print(coef_fit3_new)
  }
  if(all(abs(coef_fit3_new-coef_fit3)<tolerance)){
    break
  }
  #update model params
  coef_fit3 <- coef_fit3_new 
  plot_EM_it(het = het, fit = fit3, iteration = iteration, outpath = outpath, ID = ID)
}

##stop script if EM converges at boundary
if(offsetfit){
  print("aborted excecution, convergence at boundary")
  stop()
}
#coef(fit3, mat=T)


###plot likelihood
colnames(log_likelihood) <- c("iteration", "logLikelihood")
log_likelihood <- as.data.frame(log_likelihood)


minLL<- min(log_likelihood$logLikelihood)
maxLL<-max(log_likelihood$logLikelihood)
minLL<- minLL - abs(minLL)/10
minLL<-signif(minLL, digits = 2)
maxLL <- maxLL +abs(minLL)/10
maxLL<-signif(maxLL, digits = 2)


ll_plot<-ggplot(data = log_likelihood)+
  geom_point(aes(x = iteration, y = logLikelihood), size=2)+
  scale_x_continuous(limits=c(0,(max(log_likelihood$iteration)+2)), breaks = seq(1,(max(log_likelihood$iteration)+1),by=2), expand = c(0,0))+
  scale_y_continuous(limits=c(minLL,maxLL+((maxLL-minLL)/10)),breaks = seq(minLL,maxLL,by=((maxLL-minLL)/5)), expand=c(0,0))+
  geom_segment(x=1, y=minLL, xend=(max(log_likelihood$iteration)+1), yend=minLL, size=1.2)+
  geom_segment(x=0, y=minLL, xend=0, yend=maxLL, size=1.2)+
  coord_fixed(ratio = 17/(maxLL-minLL))+
  theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(size=14) , axis.ticks=element_line(size=1),axis.ticks.length = unit(0.25,"cm"), axis.title = element_text(size=16),legend.text = element_text(size=12), legend.position = c(0.2, 0.6),legend.title=element_blank(), legend.key = element_rect(color="white", fill=NA), legend.key.size = unit(1.5,"cm") ,panel.background = element_rect(fill=NA, color="white"))+
  guides(color = guide_legend(override.aes = list(size = 4)))

pdf(width = 16, height =8 , file = paste(outpath, ID,"_LL.pdf", sep = ""))
print(ll_plot)
dev.off()  

##Loglikelihood ratio

het$L1<-probs1
het$L2<-probs2
het$LLR<-log(het$L1/het$L2)


LLR_CO1<- log(0.99/0.01)
LLR_CO2<- log(0.01/0.99)

het$EM_class<- 0
het$EM_class[het$LLR>LLR_CO1]<-0
het$EM_class[het$LLR<LLR_CO1 &het$LLR>LLR_CO2]<-1
het$EM_class[het$LLR<LLR_CO2]<-2

het$EM_class<- factor(het$EM_class)


minLLR<-min(het$LLR)- abs((max(het$LLR)-min(het$LLR))/20)
minLLR<- signif(minLLR,2)
maxLLR<-max(het$LLR) + abs((max(het$LLR)-min(het$LLR))/20)
maxLLR<- signif(maxLLR,2)

breakstep<-floor((maxLLR-minLLR)/10)
minLLR<- minLLR-minLLR%%breakstep
maxLLR<- maxLLR + breakstep-maxLLR%%breakstep

xlabels<-c("",paste(seq(0, 0.45, by =0.05), seq(0.05, 0.5, by =0.05),sep = "-"),"")


LLR_plot <- ggplot(data = het[rs,], aes(x=Minor.allele.freq, y=LLR, color=EM_class))+
  geom_point(size=4)+
  scale_x_continuous(expand=c(0,0), limits=c(-0.015,0.515), breaks = seq(0, 0.5, by=0.25), labels=c(0,0.25,0.5))+
  scale_y_continuous(limits = c((minLLR-(maxLLR-minLLR)/20), (maxLLR+(maxLLR-minLLR)/20)), expand=c(0,0), breaks =seq(minLLR,maxLLR, by=breakstep))+
  theme(panel.background = element_rect(fill=NA, color="white"), axis.title = element_text(size=22), axis.text.y = element_text(size=18), axis.text.x=element_text(size=18,vjust=0.5), legend.position = c(0.15,0.8), legend.key = element_rect(color="white", fill=NA), legend.title = element_blank(), legend.key.size = unit(1,"cm"), legend.text = element_text(size=16))+
  geom_segment(x=0, y=(minLLR-(maxLLR-minLLR)/20), xend=0.5, yend=(minLLR-(maxLLR-minLLR)/20), size=2.4, color="black")+
  geom_segment(x=-0.015, y=minLLR, xend=-0.015, yend=maxLLR, size=2.4, color="black")+
  scale_color_manual(labels =c("high confidence single copy SNP","Uncertain", "high confidence paralogous SNP"), values=c("blue3", "darkgreen"," brown2"))+
  labs(y="logLikelihood-ratio", x="minor allele frequency (maf)")+
  geom_hline(yintercept = LLR_CO1, size=1.4)+
  annotate(geom = "text", x=0.45, y = LLR_CO2-4, label="LLR lower cutoff", size=6)+
  geom_hline(yintercept = LLR_CO2, size=1.4)+
  annotate(geom = "text", x=0.45, y = LLR_CO1+4, label="LLR upper cutoff",  size=6)

pdf(width = 16, height =8 , file = paste(outpath, ID,"_LLR.pdf", sep = ""))
print(LLR_plot)
dev.off()  

###########allelic ratios
het$allele.deviation.seed<-0
if(useRRD){
  if(verbose){
    print("testing allele ratios.....")
  }
  
  ad_mu <- mean(het$Het.allele.deviation[het$EM_class==0], na.rm=T)
  ad_sd <- max(sd(het$Het.allele.deviation[het$EM_class==0], na.rm=T),1)
  q1<-qnorm(lower.tail = F,0.025, mean = ad_mu, sd = ad_sd)
  q2<-qnorm(lower.tail = T,0.025, mean = ad_mu, sd = ad_sd)
  
  head(het)
  dplot<-ggplot(het[rs,])+
    geom_point(data = het[rs,][het$EM_class[rs] %in% c(0,1),],aes(x=Heterozygous.geno.freq, y=Het.allele.ratio, color=EM_class, shape=EM_class, fill=EM_class),stroke=1, alpha=0.2,size=0.7)+
    geom_point(data = het[rs,][het$EM_class[rs]==2, ],aes(x=Heterozygous.geno.freq, y=Het.allele.ratio, color=EM_class, shape=EM_class, fill=EM_class),stroke=1, alpha=0.2,size=0.7)+
    scale_x_continuous(limits=c(-0.05,1.05), breaks = seq(0,1,by =0.25), expand=c(0,0))+
    labs(x="Heterozygote frequency", y="Read ratio deviation (D)")+
    # geom_vline(xintercept = 0.1)+
    scale_y_continuous(limits = c(-0.05,1.05),expand=c(0,0), breaks =seq(0,1, by=0.25))+
    theme(panel.background = element_rect(fill=NA, color="white"), axis.title = element_text(size=16), axis.text = element_text(size=14), legend.title = element_blank(), legend.position = c(0.5,0.9), legend.direction = "horizontal", legend.background = element_rect(fill=NA, color="white"), legend.key = element_rect(fill=NA, color="white"))+
    scale_color_manual( values=c("blue3", "darkgreen"," brown2"), labels=c("single copy","unclassified", "seed"))+
    scale_fill_manual(values=c("blue3", "darkgreen"," brown2"), labels=c("single copy","unclassified", "seed"))+
    scale_shape_manual(values = c(0,1,2), labels=c("single copy","unclassified", "seed"))+
    geom_segment(x=-0.05, xend=-0.05, y=0, yend=1, color="black")+
    geom_segment(x=0, xend=1, y=-0.05, yend=-0.05, color="black")+
    guides(color = guide_legend(override.aes = list(size=5, alpha=1)))
  
  dens1<-ggplot(het[rs,], aes(x =Het.allele.ratio , color = EM_class, fill=EM_class)) + 
    geom_density(alpha = 0.4) + 
    scale_x_continuous(limits = c(-0.05,1.05),expand=c(0,0), breaks =seq(0,1, by=0.25))+
    scale_color_manual( values=c("blue3", "darkgreen"," brown2"), labels=c("single copy","unclassified", "seed"))+
    scale_fill_manual(values=c("blue3", "darkgreen"," brown2"), labels=c("single copy","unclassified", "seed"))+
    theme_void() + 
    theme(legend.position = "none")+
    coord_flip()
  
  
  
  ARplot<- dplot + dens1 + 
    plot_layout(ncol = 2, nrow = 1, widths = c(4, 1), heights = c(1, 4))
  
  pdf(width = 16, height =8 , file = paste(outpath, ID,"_AR.pdf", sep = ""))
  print(ARplot)
  dev.off()  
  
  # head(het)
  dplot<-ggplot(het[rs,])+
    geom_point(data = het[rs,][het$EM_class[rs] %in% c(0,1),],aes(x=Heterozygous.geno.freq, y=Het.allele.deviation, color=EM_class, shape=EM_class, fill=EM_class),stroke=1, alpha=0.2,size=0.7)+
    geom_point(data = het[rs,][het$EM_class[rs]==2, ],aes(x=Heterozygous.geno.freq, y=Het.allele.deviation, color=EM_class, shape=EM_class, fill=EM_class),stroke=1, alpha=0.2,size=0.7)+
    scale_x_continuous(limits=c(-0.05,1.05), breaks = seq(0,1,by =0.25), expand=c(0,0))+
    labs(x="Heterozygote frequency", y="Read ratio deviation (D)")+
    # geom_vline(xintercept = 0.1)+
    geom_hline(yintercept = q1)+
    geom_hline(yintercept = q2)+
    scale_y_continuous(limits = c(-85,85),expand=c(0,0), breaks =seq(-80,80, by=20))+
    theme(panel.background = element_rect(fill=NA, color="white"), axis.title = element_text(size=16), axis.text = element_text(size=14), legend.title = element_blank(), legend.position = c(0.5,0.9), legend.direction = "horizontal", legend.background = element_rect(fill=NA, color="white"), legend.key = element_rect(fill=NA, color="white"))+
    scale_color_manual( values=c("blue3", "darkgreen"," brown2"), labels=c("single copy","unclassified", "seed"))+
    scale_fill_manual(values=c("blue3", "darkgreen"," brown2"), labels=c("single copy","unclassified", "seed"))+
    scale_shape_manual(values = c(0,1,2), labels=c("single copy","unclassified", "seed"))+
    geom_segment(x=-0.05, xend=-0.05, y=-80, yend=80, color="black")+
    geom_segment(x=0, xend=1, y=-85, yend=-85, color="black")+
    guides(color = guide_legend(override.aes = list(size=5, alpha=1)))
  
  dens1<-ggplot(het[rs,], aes(x =Het.allele.deviation , color = EM_class, fill=EM_class)) + 
    geom_density(alpha = 0.4) + 
    scale_x_continuous(limits = c(-85,85),expand=c(0,0), breaks =seq(-80,80, by=20))+
    scale_color_manual( values=c("blue3", "darkgreen"," brown2"), labels=c("single copy","unclassified", "seed"))+
    scale_fill_manual(values=c("blue3", "darkgreen"," brown2"), labels=c("single copy","unclassified", "seed"))+
    theme_void() + 
    theme(legend.position = "none")+
    coord_flip()
  
  
  RRDplot<- dplot + dens1 + 
    plot_layout(ncol = 2, nrow = 1, widths = c(4, 1), heights = c(1, 4))
  
  
  pdf(width = 16, height =8 , file = paste(outpath, ID,"_RRD.pdf", sep = ""))
  print(RRDplot)
  dev.off()  
  
  if(verbose){
    print("creating normal CI with mu;")
    print(ad_mu)
    print("and sd;")
    print(ad_sd)
    print("upper CI limit;")
    print(q1)
    print("lower CI limit;")
    print(q2)
  }
  het$allele.deviation.seed[het$EM_class==1 & (het$Het.allele.deviation > q1 | het$Het.allele.deviation < q2) ]<-1
  het$EM_class[het$EM_class==1 & (het$Het.allele.deviation > q1 | het$Het.allele.deviation < q2) ] <- 2
}



write.table(x = het,file =paste(outpath, ID,"_EMresults.het", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
############### distance part
if(verbose){
  print("calculate distance cutoff.....")
}

chrlist<-unique(het$Chromosome)

distances<-c()
for(j in chrlist){
  chr<- j
  het_par<- het[het$EM_class==2 &het$Chromosome==chr ,]
  for(i in 2:nrow(het_par)){
    distances<- c(distances,het_par$Position[i]-het_par$Position[(i-1)])
  }
}
ldist<- length(distances)
distances<-distances-1

##sample median cutoff from
dist_cutoff_samples <- c()

for(i in 1:dist_em_rep){
  cutoff<-NA
  result_dist <- em_algorithm_geom(distances)
  p1 <- result_dist[[1]]
  p2 <- result_dist[[2]]
  w <- result_dist[[3]]
  # sum(w>0.5) / length(w)
  # sum(w<=0.5) / length(w)
  i<-0
  if(p1<p2){
    tmp1<-p1
    p1<-p2
    p2<-tmp1
  }
  
  md <- max(distances)
  x_dgeom <- seq(0, md, by = 1)
  y1_dgeom <- dgeom(x_dgeom, prob = p1)    
  y2_dgeom <- dgeom(x_dgeom, prob = p2)    
  
  CF<-TRUE
  for(i in 0:1000){
    d1 <- dgeom(i, prob=p1)
    d2 <- dgeom(i, prob=p2)
    if((d2/d1)> 0.05){
      if(CF){
        conservative_cutoff<- i
        CF<-FALSE
      }
      if((d2/d1)> 1){
        cutoff<- i
        break
      }
    }
  }

  if(cdist){
    dist_cutoff <- conservative_cutoff
  }else{
    dist_cutoff <- cutoff
  }
  dist_cutoff_samples <- rbind(dist_cutoff_samples, c(dist_cutoff,p1,p2 ))
}
dist_cutoff_samples <- as.data.frame(dist_cutoff_samples)
dist_cutoff <- round(median(dist_cutoff_samples$V1, na.rm = T))

closestDist<- abs(dist_cutoff_samples$V1-dist_cutoff)
closestDist[is.na(closestDist)]<- Inf
closestDistIndex <- which(min(closestDist)==closestDist)[1]
p1<- dist_cutoff_samples$V2[closestDistIndex]
p2<- dist_cutoff_samples$V3[closestDistIndex]

dist_cutoff <- max(dist_cutoff, min_dist)
dist_cutoff <- min(dist_cutoff, max_dist)

y1_dgeom <- dgeom(x_dgeom, prob = p1)    
y2_dgeom <- dgeom(x_dgeom, prob = p2)    




distances_factor<- factor(distances, levels=c(0:md))
distfreq<-tabulate(distances_factor)/length(distances_factor)

disttab<-as.data.frame(rbind(cbind(c(0:md),distfreq, rep("empirical", (md+1))),cbind(c(0:md),y1_dgeom, rep("within haplotype", (md+1))), cbind(c(0:md),y2_dgeom, rep("between haplotype", (md+1)))))


colnames(disttab)<- c("bin", "freq", "dist")
disttab$bin<- as.numeric(disttab$bin)
disttab$freq<- as.numeric(disttab$freq)


#prepare scale for plotting
mdf<- max(disttab$freq)
for(i in 1:100000){
  if(10^(i)*mdf >=1){
    nks<-i
    break
  }
}
mdf_y <- mdf*1.1
mdf_y <- round(mdf_y, digits = nks+1)
mindf_y <- round((mdf_y - mdf), digits = nks+1)
modulo4 <- (mdf_y*10^(nks+1)) %% 4
if(modulo4 !=0){
  mdf_y <- mdf_y + (4-modulo4) * 10^(-(nks+1))
  mdf_y_break <- mdf_y -modulo4 * 10^(-(nks+1))
}else{
  mdf_y_break <- mdf_y -4 * 10^(-(nks+1))
}
  
modulo2 <- (mindf_y*10^(nks+1)) %% 2
if(modulo2 !=0){
  mindf_y <- mindf_y + (2-modulo2) * 10^(-(nks+1))  
}
mindf_y <- -mindf_y


dist_plot<-ggplot(data=disttab[disttab$bin<1001,], aes(x = bin, y=freq, color=dist,goup=dist))+
  geom_point()+
  geom_path()+
  scale_x_continuous(expand=c(0,0),limits = c(-50,1050))+
  scale_y_continuous(expand=c(0,0), limits=c(mindf_y, mdf_y), breaks = seq(0,mdf_y_break, by=(8 * 10^(-(nks+1))) ))+
  # geom_vline(xintercept = dist_cutoff, size=1, linetype="dashed")+
  labs(x="distance",y="relative frequency")+
  theme(panel.background = element_rect(fill=NA, color="white"), axis.title = element_text(size=22), axis.text.y = element_text(size=18), axis.text.x=element_text(size=18,vjust=0.5), legend.position = c(0.55,0.8), legend.key = element_rect(color="white", fill=NA), legend.title = element_blank(), legend.key.size = unit(1,"cm"), legend.text = element_text(size=16))+
  #scale_color_manual(values=c("darkorange","darkgreen"))+
  scale_color_manual(values=c("purple","brown2", "darkorange"),breaks = c("empirical", "within haplotype", "between haplotype"))+
  # scale_shape_manual(values=c(1,2, 3),breaks = c("empirical", "within haplotype", "between haplotype"))+
  geom_segment(x=0, y=mindf_y, xend=1000, yend=mindf_y, size=2, color="black")+
  geom_segment(x=-50, y=0, xend=-50, yend=mdf_y_break, size=2, color="black")+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  geom_vline(xintercept = dist_cutoff)+
  coord_cartesian()


#EM for geometric distributio


pdf(width = 16, height =8 , file = paste(outpath, ID,"_dist.pdf", sep = ""))
print(dist_plot)
dev.off() 

print(paste("distance cutoff: ", dist_cutoff, sep = ""))
if(verbose){
  print("writing results.....")
}

write.table(x = dist_cutoff,file =paste(outpath, ID,"_EMresults.dist", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
