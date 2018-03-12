
#### Choisir la pr?valence pour les 3 derni?res figures ####
#        Fig.1  Fig.2  Fig.3          
FN =  c(     0,     0,    15); 
TP =  c(   100,   85,    85)
FP =  c(    300,   15,     0)

#Pr?valence :  Fig.4  Fig.5  Fig.6 
prevalence = c(  0.25,  0.1,   0.01 )

#OPR & UTP ? fixer :
OPR = 0.3
UTP = 0.3

total = 10000
obs  = TP+FN
obs  = c(obs, round( prevalence * total ))
FN   = c(FN,  round( prevalence * total * UTP ))
TP   = c(TP,  round( prevalence * total * (1-UTP) ))
FP   = c(FP,  round( TP[4:6] * (OPR/(1-OPR)) ))
pred = TP+FP
TN   = total-obs-FP





#### Calcul des m?triques ####
Jaccard  = round(   TP / (  TP+FN+FP)  ,2) 
Sorensen = round( 2*TP / (2*TP+FN+FP)  ,2)
OPR  = round( FP / pred ,2)
UTP  = round( FN / obs  ,2)
Sens = round( TP / (TP+FN) ,2)
Spe  = round( TN / (TN+FP) ,2)
TSS  = round( Sens + Spe - 1 ,2)
   x = TP-pred*obs/total
# TSS   = round( x / (x +FP*(obs/total) +FN*(1-obs/total))  ,2)
Kappa = round( x / (x +FP*(1/2)       +FN*(1/2)        )  ,2) 
Prevalence = round( (TP+FN) / total  ,2)

df <- data.frame(TP, TN, FP, FN, obs, pred, OPR, UTP,Jaccard, Sorensen, Sens, Spe, TSS, Kappa, Prevalence )
print(df)


#### Param?tres graphiques ####
width  = 1200 #pixels
height = 675  #pixels
row = 2
col = 3

titles=c(paste0(FP[1]/TP[1] * 100, '% overprediction,\n0% under-prediction,\nprevalence = 0.01'), 
         '15% over-prediction,\n0% under-prediction,\nprevalence = 0.01',
         '0% over-prediction,\n15% under-prediction,\nprevalence = 0.01',
         paste('30% over- & under-prediction,\nprevalence = ', df$Prevalence[4],sep=""), 
         paste('30% over- & under-prediction,\nprevalence = ', df$Prevalence[5],sep=""), 
         paste('30% over- & under-prediction,\nprevalence = ' , df$Prevalence[6],sep=""))


#### Param?tres d'affichage des cercles ####
area=height*width
diametres.obs <-  2*sqrt(obs*area /(total*pi))/96 #inches
diametres.pred <- 2*sqrt(pred*area/(total*pi))/96 #inches
ratio = df$UTP

# png( 'C:/Users/BorisMNHN/Google Drive/recherche/publis/21 - Metriques virtual species/Figure 1.png', width = width, height = height )  
png( './outputs/Figure 1.png', width = width * 4.2, height = height * 4.2, res = 300)  
par( mfrow = c(row,col), mar = c(0,0,0,0) )
#### Plots ####
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

if(i <= 3)
{
  adj <- .23
  lett <- ".\n\n"
} else
{
  adj <- .46
  lett <- ".\n"
}
mtext(paste(letters[i], lett, sep=""), adj=0.03, side=1, line=-1.5, cex=1.6, font=1)
mtext(titles[i], adj=.5, side=1, line=-1.5, cex=1.6, font=1)
# mtext(paste('TN =',   TN[i], '\nTP =',       TP[i], '\nFN =',            FN[i], '\nFP =',             FP[i], sep=" "), adj=0.02, side=3, line=-7, cex=1)
mtext(paste('TSS = ', format(TSS[i], nsmall = 2),
            '\nSensitivity = ', format(Sens[i], nsmall = 2),
            '\nSpecificity = ',   format(Spe[i], nsmall = 2), sep = ""), adj=0.05, side=3, line=-7, cex=1.6)
mtext(paste(#'\nJaccard = ',  format(Jaccard[i], nsmall = 2),
            '\nS\u00F8rensen = ',  format(Sorensen[i], nsmall = 2),
            '\nUPR = ', format(UTP[i], nsmall = 2), 
            '\nOPR = ', format(OPR[i], nsmall = 2), sep=""), adj=0.95, side=3, line=-7, cex=1.6)
# mtext(paste('UTP = ', format(UTP[i], nsmall = 2), 
#             '\nOPR = ',     format(OPR[i], nsmall = 2), 
#             '\nSensitivity = ', 
#             '\nSpecificity = ',   format(Spe[i], nsmall = 2), sep=""), adj=0.95, side=1, line=-2, cex=1)
}

dev.off()

