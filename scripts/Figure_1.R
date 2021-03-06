
#### Defining confusion matrices for figure 1 ####
#        Fig.1  Fig.2  Fig.3          
FN =  c(     0,    0,    15)
TP =  c(   100,   85,    85)
FP =  c(   300,   15,     0)

#Prevalence for  Fig.4  Fig.5  Fig.6 
prevalence = c(  0.25,  0.1,   0.01 )

# Equal amounts of overprediction and underprediction for fig 4/5/6
overpred = 0.3
underpred = 0.3

total = 10000
obs  = TP+FN
obs  = c(obs, round( prevalence * total ))
FN   = c(FN,  round( prevalence * total * underpred ))
TP   = c(TP,  round( prevalence * total * (1-underpred) ))
FP   = c(FP,  round( TP[4:6] * (overpred/(1-overpred)) ))
pred = TP+FP
TN   = total-obs-FP





#### Metric calculation ####
Jaccard  = round(   TP / (  TP+FN+FP)  ,2) 
Sorensen = round( 2*TP / (2*TP+FN+FP)  ,2)
OPR  = round( FP / pred ,2)
UTP  = round( FN / obs  ,2)
Sens = round( TP / (TP+FN) ,2)
Spe  = round( TN / (TN+FP) ,2)
TSS  = round( Sens + Spe - 1 ,2)
   x = TP-pred*obs/total
Kappa = round( x / (x +FP*(1/2)       +FN*(1/2)        )  ,2) 
Prevalence = round( (TP+FN) / total  ,2)

df <- data.frame(TP, TN, FP, FN, obs, pred, OPR, UTP,Jaccard, Sorensen, Sens, Spe, TSS, Kappa, Prevalence )
print(df)


#### Graphical parameters ####
width  = 1200 #pixels 1200 (png), inches(postscript)
height = 675  #pixels 675 (png), inches(postscript)
row = 2
col = 3

titles=c(paste0(OPR[1] * 100, '% overprediction,\n0% under-prediction,\nprevalence = 0.01'), 
         '15% over-prediction,\n0% under-prediction,\nprevalence = 0.01',
         '0% over-prediction,\n15% under-prediction,\nprevalence = 0.01',
         paste('30% over- & under-prediction,\nprevalence = ', df$Prevalence[4],sep=""), 
         paste('30% over- & under-prediction,\nprevalence = ', df$Prevalence[5],sep=""), 
         paste('30% over- & under-prediction,\nprevalence = ' , df$Prevalence[6],sep=""))


#### Circle parameters ####
area=height*width
diametres.obs <-  2*sqrt(obs*area /(total*pi))/96 #inches
diametres.pred <- 2*sqrt(pred*area/(total*pi))/96 #inches
ratio = df$UTP

# png( './outputs/Figure 1.png', width = width * 4.2, height = height * 4.2, res = 300)  
# cairo_ps("./outputs/figure1.ps", width = 12, height = 6.75)
cairo_pdf("./outputs/figure1.pdf", width = 12, height = 6.75)
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
  
  
# Observed range
symbols(  x0,y0, circles=diametres.obs[i],  ylim=c(-10,10), xlim=c(-10,10), inches = FALSE,  xaxt='n', yaxt='n', ann=FALSE, bg=gray(0.1, alpha = .5))
# Predicted range
symbols(  x1,y0, fg = grey(.3), circles=diametres.pred[i], ylim=c(-10,10), xlim=c(-10,10), inches = FALSE,  xaxt='n', yaxt='n', ann=FALSE, bg=gray(0.85, alpha = .5), add=T)



if(i <= 3)
{
  adj <- .23
  lett <- ".\n\n"
} else
{
  adj <- .46
  lett <- ".\n"
}
mtext(paste(letters[i], lett, sep=""), adj=0.03, side=1, line=-1.5, cex=1.1, font=1)
mtext(titles[i], adj=.5, side=1, line=-1.5, cex=1.1, font=1)

mtext(paste('TSS = ', format(TSS[i], nsmall = 2),
            '\nSensitivity = ', format(Sens[i], nsmall = 2),
            '\nSpecificity = ',   format(Spe[i], nsmall = 2), sep = ""), adj=0.05, side=3, line=-5.5, cex=1.1)
mtext(paste('\nS\u00F8rensen = ',  format(Sorensen[i], nsmall = 2),
            '\nUPR = ', format(UTP[i], nsmall = 2), 
            '\nOPR = ', format(OPR[i], nsmall = 2), sep=""), adj=0.95, side=3, line=-5.5, cex=1.1)
}

dev.off()

