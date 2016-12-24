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
  
  
  results$sample.presences <- TP + FN
  results$sample.absences <- FP + TN
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
# 
# 
# ggplot(ggr, aes(x = prevalence, y = value)) + geom_point(alpha = 1/10, aes(col = TN)) +
#   facet_wrap(~variable) + stat_smooth()
# 
# 


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
# results <- results[which(results$prevalence <= 0.5), ]

ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", "UTP", "OPpc"), 
            measure.vars = c("Jaccard", "TSS", "UTP", "Sens", "OPR", "Spe"))

ggr$prevalence <- round(ggr$prevalence, 2)
ggr$prev <- as.factor(ggr$prevalence)
ggr$OPR <- ggr$OPR * 100
ggr$UTP <- ggr$UTP * 100


levels(ggr$variable) <- c("a. Jaccard", "b. True Skill Statistic", "c. Unpredicted Presences", 
                          "d. Sensitivity", "e. OverPrediction Rate", "f. Specificity")

# Modification pour masquer les couleurs sur les panneaux a c d et e
# Car il y a de l'overlap - à la demande du reviewer 3
ggr$value[which(ggr$variable %in% c("a. Jaccard", "c. Unpredicted Presences", 
                                    "d. Sensitivity", "e. OverPrediction Rate") &
                  ggr$prev != "0.01")] <- NA


tiff( './outputs/Figure 2.tiff', width = 450 * 4.2, height = 1000 * 4.2, res = 300)  
# png("./outputs/Figure 2.png", h = 800, w = 600)
ggplot(ggr, aes(x = OPR, y = value, col = prev)) + 
  geom_line(size = .9) + 
  facet_wrap(~ variable, nrow = 3) +
  xlab("Simultaneous increase in both over- and underprediction\n(% of actual presences)") + theme_bw() +
  guides(col=guide_legend(title = "Prevalence")) + ylab("Metric value") + 
  theme(legend.position = "top") +
  coord_fixed(ratio = 100)
dev.off()


# Simulation 4: overpred


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
# results <- results[which(results$prevalence <= 0.5), ]

ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", "UTP", "OPpc"), 
            measure.vars = c("Jaccard", "TSS", "UTP", "Sens", "OPR", "Spe"))

ggr$prevalence <- round(ggr$prevalence, 2)
ggr$prev <- as.factor(ggr$prevalence)
ggr$OPR <- ggr$OPR * 100
ggr$UTP <- ggr$UTP * 100
ggr$OPpc <- ggr$OPpc * 100

levels(ggr$variable) <- c("a. Jaccard", "b. True Skill Statistic", "c. Unpredicted Presences", 
                          "d. Sensitivity", "e. OverPrediction Rate", "f. Specificity")

# Modification pour masquer les couleurs sur les panneaux a c d et e
# Car il y a de l'overlap - à la demande du reviewer 3
ggr$value[which(ggr$variable %in% c("a. Jaccard", "c. Unpredicted Presences", 
                                  "d. Sensitivity", "e. OverPrediction Rate") &
                ggr$prev != "0.01")] <- NA


png("./outputs/Figure S2.2.png", h = 600, w = 1000)
ggplot(ggr, aes(x = OPpc, y = value, col = prev)) + 
  geom_line(size = .9) + 
  facet_wrap(~variable, nrow = 3) +
  xlab("Increase in overprediction\n(% of actual presences)") + theme_bw() +
  guides(col=guide_legend(title = "Prevalence")) + ylab("Metric value") + 
  theme(legend.position = "top") +
  coord_fixed(ratio = 100)
dev.off()

# Simulation 5: underpred


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


levels(ggr$variable) <- c("a. Jaccard", "b. True Skill Statistic", "c. Unpredicted Presences", 
                          "d. Sensitivity", "e. OverPrediction Rate", "f. Specificity")

# Modification pour masquer les couleurs sur les panneaux a c d et e
# Car il y a de l'overlap - à la demande du reviewer 3
ggr$value[which(ggr$variable %in% c("a. Jaccard", "b. True Skill Statistic", "c. Unpredicted Presences", 
                                    "d. Sensitivity", "e. OverPrediction Rate", "f. Specificity") &
                  ggr$prev != "0.01")] <- NA

png("./outputs/Figure S2.1.png", h = 800, w = 600)
ggplot(ggr, aes(x = UTP, y = value, col = prev)) + 
  geom_line(size = .9) + 
  facet_wrap(~variable, nrow = 3) +
  xlab("Increase in underprediction\n(% of actual presences)") + theme_bw() +
  guides(col=guide_legend(title = "Prevalence")) + ylab("Metric value") + 
  theme(legend.position = "top") +
  coord_fixed(ratio = 100)
dev.off()


# Simulation 6 : presence are fixed, increases in absences
absences <- seq(50, 9900, 10)
all.res <- NULL
res <- data.frame(matrix(nr = 1, nc = 4, dimnames = list(1, c("TP", "FP", "FN", "TN"))))
for(nbA in absences)
{
  res$TP <- 85
  res$FN <- 15
  res$FP <- floor(0.15 * nbA)
  res$TN <- ceiling(0.85 * nbA)
  all.res <- rbind(all.res, 
                   res)
}

results <- metrics(all.res, total = apply(all.res, 1, sum))
# results <- results[which(results$prevalence <= 0.5), ]

ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", 
                                 "UTP", "OPpc", "sample.presences", "sample.absences"), 
            measure.vars = c("Jaccard", "TSS", "UTP", "Sens", "OPR", "Spe"))

ggr$prevalence <- round(ggr$prevalence, 2)
ggr$prev <- as.factor(ggr$prevalence)
ggr$OPR <- ggr$OPR * 100
ggr$UTP <- ggr$UTP * 100


levels(ggr$variable) <- c("a. Jaccard", "b. True Skill Statistic", "c. Unpredicted Presences", 
                          "d. Sensitivity", "e. OverPrediction Rate", "f. Specificity")


png("./outputs/Figure S2.3.png", h = 800, w = 600)
ggplot(ggr, aes(x = sample.absences, y = value)) + 
  geom_line(size = .9) + 
  facet_wrap(~variable, nrow = 3) +
  xlab("Number of absences\n(Number of presences = 100)") + theme_bw() +
  guides(col=guide_legend(title = "Prevalence")) + ylab("Metric value") + 
  theme(legend.position = "top") #+
  # coord_fixed(ratio = 100)
dev.off()

# Simulation 7 : comparing metrics for a good and a bad model, with different sample sizes
presences <- rep(100, 3)
absences <- c(100, 1000, 10000)
true.prev <- c(.01, .10, .5)
area <- 100000
true.presences <- true.prev * area
true.absences <- area - (true.prev * area)
models <- list(M1 = list(TP = .85,
                         FP = .15,
                         FN = .15),
               M2 = list(TP = .60,
                         FP = .40,
                         FN = .40),
               M3 = list(TP = .80,
                         FP = 0,
                         FN = .20),
               M4 = list(TP = .80,
                         FP = .20,
                         FN = 0))

all.res <- NULL

sim <- 10

res <- data.frame(matrix(nr = 1, nc = 6, dimnames = list(1, c("model", "sp.prev", "TP", "FP", "FN", "TN"))))
for(i in 1:sim)
{
  for(nbP in presences)
  {
    for(nbA in absences)
    {
      for(md in 1:length(models))
      {
        for(sp.prev in 1:length(true.prev))
        {
          pred.pres <- c(rep(1, models[[md]]$TP * true.presences[sp.prev]),
                         rep(0, true.presences[sp.prev] - models[[md]]$TP * true.presences[sp.prev]))
          pred.abs <- c(rep(1, length(pred.pres) * models[[md]]$FP),
                        rep(0, true.absences[sp.prev] - length(pred.pres) * models[[md]]$FP))
          obs.pres <- sample(pred.pres, nbP)
          obs.abs <- sample(pred.abs, nbA)
          
          res$model <- md
          res$sp.prev <- true.prev[sp.prev]
          
          res$TP <- length(which(obs.pres == 1))
          res$FN <- length(which(obs.pres == 0))
          res$FP <- length(which(obs.abs == 1))
          res$TN <- length(which(obs.abs == 0))
          all.res <- rbind(all.res, 
                           res)
        }      
      }
    }
  }
  
}

results <- metrics(all.res, total = apply(all.res[, 3:6], 1, sum))
results$model <- all.res$model
results$sp.prev <- all.res$sp.prev
results$c.val <- results$sample.presences / ((results$sample.absences + results$sample.presences) * results$sp.prev)
results$Fcpb <- 2 * results$TP / (results$TP + results$FN + results$FP * results$c.val)
# results <- results[which(results$prevalence <= 0.5), ]

ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", 
                                 "UTP", "OPpc", "sample.presences", "sample.absences", "model", "sp.prev"), 
            measure.vars = c("Jaccard", "Fcpb", "TSS"))

ggr$prevalence <- round(ggr$prevalence, 2)
ggr$prev <- as.factor(ggr$prevalence)
ggr$OPR <- ggr$OPR * 100
ggr$UTP <- ggr$UTP * 100

ggr$sample.size <- ggr$sample.absences + ggr$sample.presences
ggr$model <- as.factor(ggr$model)
ggr$sample.absences <- as.factor(ggr$sample.absences)
levels(ggr$model) <- c("OPR = UP = 0.15", "OPR = UP = 0.40", "OPR = 0, UP = 0.20", "OPR = 0.20, UP = 0")

levels(ggr$variable) <- c("a. Jaccard", "b. Fcpb", "c. True Skill Statistic")


# png("./outputs/Figure S2.3.png", h = 800, w = 600)
ggplot(ggr, aes(x = sp.prev, y = value, col = model)) + 
  geom_point(aes(shape = sample.absences)) + 
  stat_smooth(se = T) +
  facet_wrap(~variable, nrow = 1) +
  # xlab("Number of absences\n(Number of presences = 100)") +
  theme_bw() +
  # guides(col=guide_legend(title = "Prevalence")) + ylab("Metric value") + 
  theme(legend.position = "top") #+
# coord_fixed(ratio = 100)
# dev.off()

# Clearly, the results indicate that none of the metrics can provide
# a value independent of prevalence.
# The problem is most likely linked to the sampling of absences which
# is dependent on species prevalence. For TSS, low prevalences mean
# high numbers of TN and thus high values of sensitivity --> high TSS.
# For Jaccard, low prevalences mean low likelihood of sampling 
# FP, and thus higher values. Conversely, high prevalences results
# in higher likelihood of sampling FPs, and thus lower values. 
