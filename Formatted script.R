library(reshape2)
library(ggplot2)

#### Figure 1: examples ####


#### Values for the first 3 subfigures 
# One can play around with the values to see the effect on metrics
#        Fig.1  Fig.2  Fig.3          
FN =  c(     0,    0,    15)
TP =  c(   100,   85,    85)
FP =  c(   300,   15,     0)

# Prevalence for subfigures 4 5 & 6
prevalence = c(  0.25,  0.1,  0.01 )

# OPR & UP values for subfigures 4, 5 & 6
OPR = 0.55
UTP = 0.3

# Calculation of all the metrics
total = 10000 # (total number of possible values) - change this number to see the effect of sample size on metrics
obs  = TP+FN
obs  = c(obs, round( prevalence * total ))
FN   = c(FN,  round( prevalence * total * UTP ))
TP   = c(TP,  round( prevalence * total * (1-UTP) ))
FP   = c(FP,  round( TP[4:6] * (OPR/(1-OPR)) ))
pred = TP+FP
TN   = total-obs-FP


Jaccard  = round(   TP / (  TP+FN+FP)  ,2) 
Sorensen = round( 2*TP / (2*TP+FN+FP)  ,2)
OPR  = round( FP / pred ,2)
UTP  = round( FN / obs  ,2)
Sens = round( TP / (TP+FN) ,2)
Spe  = round( TN / (TN+FP) ,2)
TSS  = round( Sens + Spe - 1 ,2)
x = TP-pred*obs/total
Prevalence = round( (TP+FN) / total  ,2)

df <- data.frame(TP, TN, FP, FN, obs, pred, OPR, UTP,Jaccard, Sorensen, Sens, Spe, TSS, Prevalence )
print(df)


#### Graphical parameters 
width  = 1200 #pixels
height = 675  #pixels
row = 2
col = 3

titles=c(paste0('High overprediction (',
                FP[1]/TP[1] * 100, '%)'), 
         'Good fit, 15% over-prediction', 'Good fit, 15% under-prediction',
         paste('Prevalence = ', df$Prevalence[4],sep=""), 
         paste('Prevalence = ', df$Prevalence[5],sep=""), 
         paste('Prevalence = ' , df$Prevalence[6],sep=""))


#### Calculation of the relative size of circles for figures 
area=height*width
diametres.obs <-  2*sqrt(obs*area /(total*pi))/96 #inches
diametres.pred <- 2*sqrt(pred*area/(total*pi))/96 #inches
ratio = df$UTP

# Change the following lines 
par( mfrow = c(row,col), mar = c(0,0,0,0) )
#### Plots 
for (i in 1:6)
{
  if (i<=3){
    #     L='A'
    #     j=i
    x0=-1 
    y0=-1
    x1=diametres.pred[i]-diametres.obs[i]+x0}
  else{
    #     L='B'
    #     j=i-3
    x0=-2.5
    y0=0
    x1=ratio[i]*(diametres.pred[i]+diametres.obs[i])+x0}
  
  
  # Aire Observ?e
  symbols(  x0,y0, circles=diametres.obs[i],  ylim=c(-10,10), xlim=c(-10,10), inches = FALSE,  xaxt='n', yaxt='n', ann=FALSE, bg=gray(0.1, alpha = .5))
  #Aire Pr?dite
  symbols(  x1,y0, fg = grey(.3), circles=diametres.pred[i], ylim=c(-10,10), xlim=c(-10,10), inches = FALSE,  xaxt='n', yaxt='n', ann=FALSE, bg=gray(0.85, alpha = .5), add=T)
  
  
  # Titre et affichage des m?triques
  mtext(paste(letters[i], ". ", titles[i], sep=""), adj=0.03, side=1, line=-1.5, cex=1.15, font=1)
  # mtext(paste('TN =',   TN[i], '\nTP =',       TP[i], '\nFN =',            FN[i], '\nFP =',             FP[i], sep=" "), adj=0.02, side=3, line=-7, cex=1)
  mtext(paste('TSS = ', format(TSS[i], nsmall = 2),
              '\nJaccard = ',  format(Jaccard[i], nsmall = 2), 
              '\nSorensen = ', format(Sorensen[i], nsmall = 2), sep=""), adj=0.95, side=3, line=-6, cex=1)
  # mtext(paste('UTP = ', format(UTP[i], nsmall = 2), 
  #             '\nOPR = ',     format(OPR[i], nsmall = 2), 
  #             '\nSensitivity = ', format(Sens[i], nsmall = 2), 
  #             '\nSpecificity = ',   format(Spe[i], nsmall = 2), sep=""), adj=0.95, side=1, line=-2, cex=1)
}


#### Figure 2: simulations ####

metrics <- function(x, total)
{
  TP <- x[, "TP"]
  FP <- x[, "FP"]
  FN <- x[, "FN"]
  TN <- x[, "TN"]
  
  pred <- TP + FP
  obs <- TP + FN
  
  
  results <- data.frame(matrix(nc = 0, nr = nrow(x)))
  
  
  
  results$TP <- TP
  results$FP <- FP
  results$FN <- FN
  results$TN <- TN
  results$Jaccard <- TP / (TP + FN + FP)
  results$Sorensen = 2 * TP / (2 * TP + FN + FP)
  results$OPR  = FP / (TP + FP)
  results$UTP  = FN / (TP + FN)
  results$Sens = TP / (TP + FN)
  results$Spe  = TN / (TN + FP)
  results$TSS  = TP / (TP + FN) + TN / (TN + FP) - 1
  results$Kappa = (TP - pred * obs / total) / ((TP - pred * obs / total) + FP * (1/2) + FN * (1/2))
  results$prevalence = obs / total
  results$pred.prevalence = pred / total
  results$OPpc = FP / obs
  
  return(results)
  
}


# Figure 2
tns <- c(9900, 1900, 900, 300, 100, 33, 11, 5, 1)
all.res <- NULL
res <- data.frame(matrix(nr = 1, nc = 4, dimnames = list(1, c("TP", "FP", "FN", "TN"))))
for(TN in tns)
{
  res$TP <- 100
  res$FP <- 0
  res$FN <- 0
  res$TN <- TN
  all.res <- rbind(all.res, 
                   res)
  
  while(res$TP > 0)
  {
    if(res$TN > 0)
    {
      res$TP <- res$TP - 1
      res$FN <- res$FN + 1
      res$FP <- res$FP + 1
      res$TN <- res$TN - 1
      all.res <- rbind(all.res, 
                       res)
    } else
    {
      res$TP <- 0
    }
  }
}

results <- metrics(all.res, total = apply(all.res, 1, sum))

ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", "UTP", "OPpc"), 
            measure.vars = c("Jaccard", "TSS", "UTP", "Sens", "OPR", "Spe"))

ggr$prevalence <- round(ggr$prevalence, 2)
ggr$prev <- as.factor(ggr$prevalence)
ggr$OPR <- ggr$OPR * 100
ggr$UTP <- ggr$UTP * 100

levels(ggr$variable) <- c("Jaccard", "True Skill Statistic", "Unpredicted Presences", "Sensitivity", "OverPrediction Rate", "Specificity")

ggplot(ggr, aes(x = OPR, y = value, col = prev)) + 
  geom_line(size = .9) + 
  facet_wrap(~variable, scales = "free_y", nrow = 3) +
  xlab("Simultaneous increase in both over- and underprediction\n(% of actual presences)") + theme_bw() +
  guides(col=guide_legend(title = "Prevalence")) + ylab("Metric value") + 
  theme(legend.position = "top")


# Figure A1

tns <- c(9900, 1900, 900, 300, 100, 33, 11, 5, 1)
all.res <- NULL
res <- data.frame(matrix(nr = 1, nc = 4, dimnames = list(1, c("TP", "FP", "FN", "TN"))))
for(TN in tns)
{
  res$TP <- 100
  res$FP <- 0
  res$FN <- 0
  res$TN <- TN
  all.res <- rbind(all.res, 
                   res)
  
  while(res$TP > 0)
  {
    
    res$TP <- res$TP - 1
    res$FN <- res$FN + 1
    all.res <- rbind(all.res, 
                     res)
    
  }
}

results <- metrics(all.res, total = apply(all.res, 1, sum))
# results <- results[which(results$prevalence <= 0.5), ]

ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", "UTP", "OPpc"), 
            measure.vars = c("Jaccard", "TSS", "UTP", "Sens", "OPR", "Spe"))

ggr$prevalence <- round(ggr$prevalence, 2)
ggr$prev <- as.factor(ggr$prevalence)
ggr$OPR <- ggr$OPR * 100
ggr$UTP <- ggr$UTP * 100


levels(ggr$variable) <- c("Jaccard", "True Skill Statistic", "Unpredicted Presences", "Sensitivity", "OverPrediction Rate", "Specificity")

ggplot(ggr, aes(x = UTP, y = value, col = prev)) + 
  geom_line(size = .9) + 
  facet_wrap(~variable, nrow = 3) +
  xlab("Increase in underprediction\n(% of actual presences)") + theme_bw() +
  guides(col=guide_legend(title = "Prevalence")) + ylab("Metric value") + 
  theme(legend.position = "top")


# Figure A2
tns <- c(9900, 1900, 900, 300, 100, 33, 11, 5, 1)
all.res <- NULL
res <- data.frame(matrix(nr = 1, nc = 4, dimnames = list(1, c("TP", "FP", "FN", "TN"))))
for(TN in tns)
{
  res$TP <- 100
  res$FP <- 0
  res$FN <- 0
  res$TN <- TN
  all.res <- rbind(all.res, 
                   res)
  
  while(res$TN > 0 & (res$TP + res$FP) < 4 * (res$FN + res$TP))
  {
    
    res$FP <- res$FP + 1
    res$TN <- res$TN - 1
    all.res <- rbind(all.res, 
                     res)
    
  }
}

results <- metrics(all.res, total = apply(all.res, 1, sum))

ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", "UTP", "OPpc"), 
            measure.vars = c("Jaccard", "TSS", "UTP", "Sens", "OPR", "Spe"))

ggr$prevalence <- round(ggr$prevalence, 2)
ggr$prev <- as.factor(ggr$prevalence)
ggr$OPR <- ggr$OPR * 100
ggr$UTP <- ggr$UTP * 100
ggr$OPpc <- ggr$OPpc * 100

levels(ggr$variable) <- c("Jaccard", "True Skill Statistic", "Unpredicted Presences", "Sensitivity", "OverPrediction Rate", "Specificity")

ggplot(ggr, aes(x = OPpc, y = value, col = prev)) + 
  geom_line(size = .9) + 
  facet_wrap(~variable, nrow = 3) +
  xlab("Increase in overprediction\n(% of actual presences)") + theme_bw() +
  guides(col=guide_legend(title = "Prevalence")) + ylab("Metric value") + 
  theme(legend.position = "top")

