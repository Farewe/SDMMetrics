library(reshape2)
library(ggplot2)
source("./scripts/functions/basic_metrics.R")


presences <- rep(500, 3)
absences <- c(500, 1000, 10000)
true.prev <- seq(.01, .6, by = .01)
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
results$Sorp <- 2 * results$TP / (2 * results$TP + results$FN + results$FP * results$x.val)


ggr <- melt(results, id.vars = c("prevalence", "pred.prevalence", "TP", "FP", "FN", "TN", "OPR", 
                                 "UTP", "OPpc", "sample.presences", "sample.absences", "model", "sp.prev", "OPRp"), 
            measure.vars = c("TSS", "Sorensen", "Sorp"))

ggr$prevalence <- round(ggr$prevalence, 2)
ggr$prev <- as.factor(ggr$prevalence)
ggr$OPR <- ggr$OPR * 100
ggr$UTP <- ggr$UTP * 100

ggr$sample.size <- ggr$sample.absences + ggr$sample.presences
ggr$model <- as.factor(ggr$model)
ggr$sample.absences <- as.factor(ggr$sample.absences)
levels(ggr$sample.absences) <- c("Sample\nprevalence = 0.50", 
                                 "Sample\nprevalence = 0.33", 
                                 "Sample\nprevalence = 0.05")
levels(ggr$variable) <- c("a. True Skill Statistic", "b. S\u00F8rensen", "c. Prevalence calibrated\nS\u00F8rensen")

levels(ggr$model) <- c("40% overprediction &\n40% underprediction", 
                       "40% underprediction", 
                       "40% overprediction")


png("./outputs/Figure 3 presence-absence differentprev.png", h = 800, w = 960)
ggplot(ggr, aes(x = sp.prev, y = value, col = model)) +
  geom_point(alpha = 1/10) +
  stat_smooth(se = T) +
  facet_grid(sample.absences~variable, switch = 'y') +
  xlab("Species prevalence") +
  theme_bw(base_size = 20) +
  guides(col=guide_legend(title = "Case studies\n")) + ylab("Metric value") +
  theme(legend.position = "right",
        legend.key = element_rect(color = "white"),
        legend.key.height	= unit(3, 'lines'),
        legend.key.width = unit(2, 'lines')) # +
  # geom_hline(yintercept = c(6/14, 6/10), linetype = 2)
dev.off()


