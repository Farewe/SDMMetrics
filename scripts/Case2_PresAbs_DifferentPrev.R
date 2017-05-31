library(reshape2)
library(ggplot2)
source("./scripts/functions/basic_metrics.R")


presences <- rep(500, 3)
absences <- c(500, 1000, 10000)
true.prev <- c(.01, .10, .5)
area <- 100000
true.presences <- true.prev * area
true.absences <- area - (true.prev * area)
models <- list(M1 = list(TP = .60, 
                         OPR = .40,
                         UP = .40), # HAS TO BE 1 - TP
               M2 = list(TP = .60,
                         OPR = 0,
                         UP = .40),
               M3 = list(TP = 1,
                         OPR = .40,
                         UP = 0))

all.res <- NULL

sim <- 100

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
          # pred.obs.pres : Predictions for the species distribution range only
          pred.obs.pres <- c(rep(1, models[[md]]$TP * true.presences[sp.prev]),
                             rep(0, true.presences[sp.prev] - models[[md]]$TP * true.presences[sp.prev]))
          # pred.obs.abs : Predictions for areas  where the species is absent 
          pred.obs.abs <- c(rep(1, models[[md]]$OPR * models[[md]]$TP * true.presences[sp.prev] / (1 - models[[md]]$OPR)),
                            rep(0, true.absences[sp.prev] - 
                                  models[[md]]$OPR * models[[md]]$TP * true.presences[sp.prev] / (1 - models[[md]]$OPR)))
          sample.pres <- sample(pred.obs.pres, nbP)
          sample.abs <- sample(pred.obs.abs, nbA)
          
          res$model <- md
          res$sp.prev <- true.prev[sp.prev]
          
          res$TP <- length(which(sample.pres == 1))
          res$FN <- length(which(sample.pres == 0))
          res$FP <- length(which(sample.abs == 1))
          res$TN <- length(which(sample.abs == 0))
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
results$x.val <- results$sample.presences * (1 - results$sp.prev) / (results$sample.absences * results$sp.prev)
results$OPRp <- results$FP * results$x.val / (results$TP + results$FP * results$x.val)
results$Jacp <- results$TP / (results$TP + results$FN + results$FP * results$x.val)
# results <- results[which(results$prevalence <= 0.5), ]

# Fcpb by Li & Guo
results$c.val <- results$sample.presences / (results$sp.prev * results$sample.absences)
results$Fcpb <- results$TP / (results$TP + results$FN + results$FP * results$c.val)


# Number of FP in the sample for each prevalence
mean(results$OPRp[which(results$sample.absences == 500 & results$model == 1 & results$sp.prev == 0.01)])
mean(results$OPRp[which(results$sample.absences == 500 & results$model == 1 & results$sp.prev == 0.1)])
mean(results$OPRp[which(results$sample.absences == 500 & results$model == 1 & results$sp.prev == 0.5)])

ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", 
                                 "UTP", "OPpc", "sample.presences", "sample.absences", "model", "sp.prev", "OPRp"), 
            measure.vars = c("TSS", "Jaccard", "Fcpb","Jacp"))

ggr$prevalence <- round(ggr$prevalence, 2)
ggr$prev <- as.factor(ggr$prevalence)
ggr$OPR <- ggr$OPR * 100
ggr$UTP <- ggr$UTP * 100

ggr$sample.size <- ggr$sample.absences + ggr$sample.presences
ggr$model <- as.factor(ggr$model)
ggr$sample.absences <- as.factor(ggr$sample.absences)
levels(ggr$model) <- c("OPR = UP = 0.40", "OPR = 0, UP = 0.40", "OPR = 0.40, UP = 0")
levels(ggr$sample.absences) <- c("500 absence points\nSample\nprevalence = 0.50", 
                                 "1000 absence points\nSample\nprevalence = 0.33", 
                                 "10 000 absence points\nSample\nprevalence = 0.05")
levels(ggr$variable) <- c("a. True Skill Statistic", "b. Jaccard", "c. Fcpb", "d. Prevalence calibrated\nJaccard")



png("./outputs/Figure presence-absence.png", h = 800, w = 960)
ggplot(ggr, aes(x = sp.prev, y = value, col = model)) +
  geom_point() +
  stat_smooth(se = T) +
  facet_grid(sample.absences~variable) +
  xlab("Species prevalence") +
  theme_bw() +
  guides(col=guide_legend(title = "Models")) + ylab("Metric value") +
  theme(legend.position = "top") +
  geom_hline(yintercept = c(6/14, 6/10), linetype = 2)
dev.off()


