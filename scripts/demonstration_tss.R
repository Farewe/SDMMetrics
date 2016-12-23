library(virtualspecies)
library(ggplot2)
library(reshape2)
load("./Raster 1000 PA/species4.0.01")
class(species4.0.01) <- "virtualspecies"


sp <- species4.0.01


prevalences <- c(.01, .02, .03, .04, .05, .075, .1, .15, .2, .25, .3, .4, .5, .6, .7, .8, .9, .99)
sp.list.thres <- list()
sp.list.prob <- list()

for (i in prevalences)
{
  sp.list.thres[[paste0("prev.", i)]] <- convertToPA(sp, PA.method = "threshold", species.prevalence = i)
  sp.list.prob[[paste0("prev.", i)]] <- convertToPA(sp, PA.method = "probability", species.prevalence = i, alpha = -.03)
}
seq.seuils <- seq(0, 1, length = 101)
results <- array(dim = c(length(seq.seuils), 6,
                         2,
                         length(prevalences)),
                 dimnames = list(seuil = seq.seuils,
                                 metric = c("TP", "TN", "FP", "FN", "jaccard", "tss"),
                                 method = c("threshold", "probability"),
                                 prevalence = prevalences))
results.best.thres <- array(dim = c(length(prevalences), 3, 2),
                            dimnames = list(prevalence = prevalences, metric = c("jaccard", "tss", "true"),
                                            method = c("threshold", "probability")))
for(prev in prevalences)
{
  for(method in dimnames(results)$method)
  {
    if(method == "threshold")
    {
      cur.sp <- sp.list.thres[[paste0("prev.", prev)]]
    } else
    {
      cur.sp <- sp.list.prob[[paste0("prev.", prev)]]
    }
    true <- getValues(cur.sp$pa.raster)
    for(t in seq.seuils)
    {
      pred <- reclassify(cur.sp$suitab.raster, rcl = matrix(nc = 3, nr = 2, c(0, t, 0, t, 1, 1), byrow = T))
      pred.dat <- getValues(pred)
      conf.mat <- table(pred.dat, true)
      if(prod(dim(conf.mat)) == 4)
      {
        TP <- conf.mat["1", "1"]
        TN <- conf.mat["0", "0"]
        FP <- conf.mat["1", "0"]
        FN <- conf.mat["0", "1"]
        jaccard <- TP / (TP + FP + FN)
        tss <- TP / (TP + FN) + TN / (TN + FP) - 1
        results[as.character(t), , method, as.character(prev)] <- c(TP, TN, FP, FN, jaccard, tss)
      } else
      {
        results[as.character(t), , method, as.character(prev)] <- rep(NA, 6)
      }
    }
    results.best.thres[as.character(prev), "jaccard", method] <- 
      as.numeric(names(which(results[, "jaccard", method, as.character(prev)] == 
                               max(results[, "jaccard", method, as.character(prev)], 
                                   na.rm = T))))
    results.best.thres[as.character(prev), "tss", method] <-
      as.numeric(names(which(results[, "tss", method, as.character(prev)] == 
                               max(results[, "tss", method, as.character(prev)], 
                                   na.rm = T))))
    if(method == "threshold")
    {
      results.best.thres[as.character(prev), "true", method] <- as.numeric(cur.sp$PA.conversion["cutoff"])
    } else
    {
      results.best.thres[as.character(prev), "true", method] <- as.numeric(cur.sp$PA.conversion["beta"])
    }
   
  }
  cat(paste0(Sys.time(), " -- prevalence ", prev, " done.\n"))
}

gg.res.thres <- melt(results.best.thres)
png("./seuils.png", w = 1000, h = 1000)
ggplot(gg.res.thres, aes(x = prevalence, y = value, col = metric)) +
  geom_line() + facet_wrap(~method)
dev.off()

result.stack <- stack()
layernames <- NULL
for (prev in dimnames(results.best.thres)[[1]])
{
  cur.sp <- sp.list.prob[[paste0("prev.", prev)]]
  layernames <- c(layernames,
                  paste0("prev.", prev, ".true"),
                  paste0("prev.", prev, ".jaccard"),
                  paste0("prev.", prev, ".tss"))
  result.stack <- addLayer(result.stack,
                           cur.sp$pa.raster,
                           reclassify(cur.sp$suitab.raster, rcl = matrix(nc = 3, nr = 2, 
                                                                         c(0, results.best.thres[prev, "jaccard", "probability"], 0, 
                                                                           results.best.thres[prev, "jaccard", "probability"], 1, 1), 
                                                                         byrow = T)),
                           reclassify(cur.sp$suitab.raster, rcl = matrix(nc = 3, nr = 2, 
                                                                         c(0, results.best.thres[prev, "tss", "probability"], 0, 
                                                                           results.best.thres[prev, "tss", "probability"], 1, 1), 
                                                                         byrow = T)))
}
names(result.stack) <- layernames

png("./maps.png", h = 2000, w = 900)
par(mfrow = c(6, 3))
prev.to.plot <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.7)
prev.to.plot <- as.vector(sapply(prev.to.plot, function(x) c(paste0("prev.", x, ".true"),
                                                   paste0("prev.", x, ".jaccard"),
                                                   paste0("prev.", x, ".tss"))))
for (i in prev.to.plot)
{
  plot(result.stack[[i]],
       main = paste0(i))
}
dev.off()

library(RColorBrewer)

result.stack
# TP = 1 - 2 * 1 = -1
# TN = 0 - 2 * 0 = 0
# FP = 1 - 2 * 0 = 1
# FN = 0 - 2 * 1 = -2
diff.tss.0.01 <- result.stack[["prev.0.01.tss"]] - 2 * result.stack[["prev.0.01.true"]]
plot(diff.tss.0.01, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)],
     legend = F, oma = c(4.1, 4.1, 2.1, 2.1))
diff.tss.0.1 <- result.stack[["prev.0.1.tss"]] - 2 * result.stack[["prev.0.1.true"]]
plot(diff.tss.0.1, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)])
diff.tss.0.5 <- result.stack[["prev.0.5.tss"]] - 2 * result.stack[["prev.0.5.true"]]
plot(diff.tss.0.5, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)])

diff.jaccard.0.01 <- result.stack[["prev.0.01.jaccard"]] - 2 * result.stack[["prev.0.01.true"]]
plot(diff.jaccard.0.01, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)])
diff.jaccard.0.1 <- result.stack[["prev.0.1.jaccard"]] - 2 * result.stack[["prev.0.1.true"]]
plot(diff.jaccard.0.1, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)])
diff.jaccard.0.5 <- result.stack[["prev.0.5.jaccard"]] - 2 * result.stack[["prev.0.5.true"]]
plot(diff.jaccard.0.5, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)])



source("./scripts/legend_large_boxes.R")

plot(species4.0.01$suitab.raster)
title("")
mtext("a. Species probability of occurrence", side = 3, cex = 1, line = -2, outer = T, at = .5)
plot(0, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
plot(0, type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")

# postscript("./figure_papier.eps", width = 6.5, height = 5,
#            fonts = "Helvetica", family = "Helvetica", horizontal = F)
png("./figure_papier.png", w = 1400, h = 800)
par(mfrow = c(2, 3), mar = c(4.1, 2.1, 3.1, 2.1), oma = c(4.1, 0.1, 0.1, 0.1))

plot(diff.tss.0.01, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)],
     legend = F, main = "a. prevalence = 0.01\nTSS-optimised threshold")
plot(diff.tss.0.1, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)],
     legend = F, main = "b. prevalence = 0.10\nTSS-optimised threshold")
plot(diff.tss.0.5, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)],
     legend = F, main = "c. prevalence = 0.50\nTSS-optimised threshold")

plot(diff.jaccard.0.01, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)],
     legend = F, main = "e. prevalence = 0.01\nJaccard-optimised threshold")
plot(diff.jaccard.0.1, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)],
     legend = F, main = "f. prevalence = 0.10\nJaccard-optimised threshold")
legend2(x = -110, y = 10, pch = 15,
        col =  c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)],
        legend = c("TP", "TN", "FP", "FN"),
        bty = "n", border = "white", cex = 2, y.intersp = 1,
        adj = c(0, .2), title = "", xpd = NA, horiz = T)
plot(diff.jaccard.0.5, breaks = c(-2, -1.9, -.9, .1, 1), col = c("lightgrey", brewer.pal(3, "Set1"))[c(3, 4, 1, 2)],
     legend = F, main = "g. prevalence = 0.50\nJaccard-optimised threshold")
dev.off()
