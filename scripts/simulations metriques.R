library(reshape2)
library(ggplot2)
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
  results$OPpc = FP / (TP + FN)
  
  return(results)
  
}



# Simulation 1: all possible parameterisations
prevalence <- seq(0, 1, length = 100)
max.size <- 10000
step <- 100


results <- NULL
for (i in prevalence)
{
  obs.pres <- round(i * max.size)
  obs.abs <- max.size - obs.pres
  
  
  TP <- seq(0, obs.pres, by = step)
  FP <- seq(0, obs.abs, by = step)
  testing.matrix <- expand.grid(TP, FP)
  colnames(testing.matrix) <- c("TP", "FP")
  
  testing.matrix$FN <- obs.pres - testing.matrix$TP
  testing.matrix$TN <- obs.abs - testing.matrix$FP 
  
  results <- rbind(results,
                   metrics(testing.matrix, total = max.size))
}


ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN"), 
            measure.vars = c("Jaccard", "Sorensen", "OPR", "UTP", "Sens", "Spe", "TSS", "Kappa"))


ggplot(ggr, aes(x = prevalence, y = value)) + geom_point(alpha = 1/10, aes(col = TN)) +
  facet_wrap(~variable) + stat_smooth()




# Simulation 2: 
# For varying prevalences,
# variations of under-prediction (0 to -100%) or overprediction (0 to +300%) and both.

prevalence <- seq(0, 1, length = 101)
max.size <- 10000


results <- NULL
results.overpred <- NULL
results.underpred <- NULL
for (i in prevalence)
{
  obs.pres <- round(i * max.size)
  obs.abs <- max.size - obs.pres
  
  # Underpred
  underpred <- data.frame(TP = seq(obs.pres, 0, length = 101),
                          FN = seq(0, obs.pres, length = 101),
                          FP = 0,
                          TN = obs.abs)
  
  # Overpred 
  overpred <- data.frame(TP = obs.pres,
                         FN = 0,
                         FP = seq(0, min(3 * obs.pres, obs.abs), length = 101),
                         TN = obs.abs - seq(0, min(3 * obs.pres, obs.abs), length = 101))
 
  # Under + Overpred
  TP <- seq(obs.pres, 0, length = 101)
  FP <- seq(0, min(3 * obs.pres, obs.abs), length = 101)
  
  testing.matrix <- expand.grid(TP, FP)
  colnames(testing.matrix) <- c("TP", "FP")
  
  testing.matrix$FN <- obs.pres - testing.matrix$TP
  testing.matrix$TN <- obs.abs - testing.matrix$FP 
  
  results.overpred <- rbind(results.overpred,
                            metrics(overpred, total = max.size))
  results.underpred <- rbind(results.underpred,
                            metrics(underpred, total = max.size)) 
  
  results <- rbind(results,
                   metrics(testing.matrix, total = max.size))
}


# Simulation 3: all possible parameterisations


a <- function(x) (100 - x * 100)/x


tns <- c(9900, 1900, 900, 300, 100, 33, 11, 5, 1)
for(TN in tns)
{
  res$TP <- 100
  res$FP <- 0
  res$FN <- 0
  
  while(res$TP > 0)
  {
    
  }
}
results <- NULL
for (i in prevalence)
{
  obs.pres <- round(i * max.size)
  obs.abs <- max.size - obs.pres
  
  
  TP <- seq(0, obs.pres, by = step)
  FP <- seq(0, obs.abs, by = step)
  testing.matrix <- expand.grid(TP, FP)
  colnames(testing.matrix) <- c("TP", "FP")
  
  testing.matrix$FN <- obs.pres - testing.matrix$TP
  testing.matrix$TN <- obs.abs - testing.matrix$FP 
  
  results <- rbind(results,
                   metrics(testing.matrix, total = max.size))
}


ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN"), 
            measure.vars = c("Jaccard", "Sorensen", "OPR", "UTP", "Sens", "Spe", "TSS", "Kappa"))


ggplot(ggr, aes(x = prevalence, y = value)) + geom_point(alpha = 1/10, aes(col = TN)) +
  facet_wrap(~variable) + stat_smooth() + coord_fixed(ratio = 1)



ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", "UTP"), 
            measure.vars = c("Jaccard", "Sorensen", "Sens", "Spe", "TSS", "Kappa"))
ggr1 <- ggr[ggr$prevalence %in% c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), ]

png("./outputs/Metriques vs OPR - covariation OPR et UTP.png", h = 600, w = 800)
ggplot(ggr1, aes(x = OPR, y = value, col = prevalence)) + 
  geom_point(alpha = 1/10) + 
  facet_wrap(~variable, scales = "free_y")  +
  coord_fixed(ratio = 1)
dev.off()

png("./outputs/Metriques vs UTP - covariation OPR et UTP.png", h = 600, w = 800)
ggplot(ggr1, aes(x = UTP, y = value, col = prevalence)) + 
  geom_point(alpha = 1/10) + 
  facet_wrap(~variable, scales = "free_y") +
  coord_fixed(ratio = 1) 
dev.off()

ggr.u <- melt(results.underpred, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", "UTP"), 
            measure.vars = c("Jaccard", "Sorensen", "Sens", "Spe", "TSS", "Kappa"))
ggr.u1 <- ggr.u[ggr.u$prevalence %in% c(0.01, 0.05, 0.10, 0.15, 0.2, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), ]
ggr.u1$prev <- as.factor(ggr.u1$prevalence)

png("./outputs/Metriques vs UTP - OPR = 0.png", h = 600, w = 800)
ggplot(ggr.u1, aes(x = UTP, y = value, col = prev)) + 
  geom_line() + 
  facet_wrap(~variable, scales = "free_y") + xlab("Underprediction (% of observed presences)") +
  coord_fixed(ratio = 1)
dev.off()
ggr.o <- melt(results.overpred, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", "UTP", "OPpc"), 
              measure.vars = c("Jaccard", "Sorensen", "Sens", "Spe", "TSS", "Kappa"))
ggr.o1 <- ggr.o[ggr.o$prevalence %in% c(0.01, 0.05, 0.10, 0.15, 0.2, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), ]
ggr.o1$prev <- as.factor(ggr.o1$prevalence)
png("./outputs/Metriques vs OPR - UTP = 0.png", h = 600, w = 800)
ggplot(ggr.o1, aes(x = OPR, y = value, col = prev)) + 
  geom_line() + 
  facet_wrap(~variable, scales = "free_y") +
  xlab("Overprediction (as % of predicted presences)") +
  coord_fixed(ratio = 1)
dev.off()

png("./outputs/Metriques vs OPpc - UTP = 0.png", h = 600, w = 800)
ggplot(ggr.o1, aes(x = OPpc, y = value, col = prev)) + 
  geom_line() + 
  facet_wrap(~variable, scales = "free_y") +
  xlab("Overprediction (as % of observed presences)") +
  coord_fixed(ratio = 1)
dev.off()

summary(results.overpred)

apply(results.overpred[, c("TP", "FP", "FN", "TN")], 1, sum)
