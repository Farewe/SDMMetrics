library(ggplot2)
library(raster)
source("i:/r/scripts/functions/summaryse.R")
source("i:/r/scripts/functions/multiplot.R")

results <- read.csv("./data/nsb_result.csv")

# Each 'realised' species is dependent on 3 columns: vs_size, vs_no & scenario
unique(cbind(results$vs_size, results$vs_no, results$scenario))

# For models, they replicated sampled occurrence points for each species 10 times
# (column replication_no)

# Calculating TSS values
results$TSS_A <- results$sensitivity_A + results$specificity_A - 1
results$TSS_Go <- results$sensitivity_Go + results$specificity_Go - 1
results$size_A <- results$A1 + results$A2


sp1 <- results[which(results$vs_size == "big" & 
                       results$vs_no == 1 &
                       results$scenario ==  "CB"), ]

fig3Atab <- data.frame(summarySE(data = results,
                                measurevar = "TSS_A",
                                groupvars = c("algorithm", "scenario")),
                      summarySE(data = results,
                                measurevar = "TSS_Go",
                                groupvars = c("algorithm", "scenario"))[, c("TSS_Go", "sd", "se", "ci")],
                      summarySE(data = results,
                                measurevar = "size_A",
                                groupvars = c("algorithm", "scenario"))$size_A,
                      summarySE(data = results,
                                measurevar = "A1",
                                groupvars = c("algorithm", "scenario"))$A1)
colnames(fig3Atab)[colnames(fig3Atab) %in% c("sd", "se", "ci", "sd.1", "se.1", "ci.1")] <- c("sd_A", "se_A", "ci_A", "sd_Go", "se_Go", "ci_Go")
colnames(fig3Atab)[((ncol(fig3Atab) - 1):ncol(fig3Atab))] <- c("size_A", "size_Go")

p1 <- ggplot(fig3Atab, aes(x = TSS_Go, y = TSS_A, col = algorithm)) + 
  geom_point(aes(shape = scenario, size = size_A))  +
  xlim(0.9, 1) + geom_errorbar(aes(ymin = TSS_A - ci_A, ymax = TSS_A + ci_A)) +
  geom_errorbarh(aes(xmin = TSS_Go - ci_Go, xmax = TSS_Go + ci_Go)) +
  ylim(0, 1) + theme_bw()

fig3Btab <- data.frame(summarySE(data = results,
                                 measurevar = "TSS_A",
                                 groupvars = c("algorithm", "vs_size")),
                       summarySE(data = results,
                                 measurevar = "TSS_Go",
                                 groupvars = c("algorithm", "vs_size"))[, c("TSS_Go", "sd", "se", "ci")],
                       summarySE(data = results,
                                 measurevar = "size_A",
                                 groupvars = c("algorithm", "vs_size"))$size_A,
                       summarySE(data = results,
                                 measurevar = "A1",
                                 groupvars = c("algorithm", "vs_size"))$A1)
colnames(fig3Btab)[colnames(fig3Btab) %in% c("sd", "se", "ci", "sd.1", "se.1", "ci.1")] <- c("sd_A", "se_A", "ci_A", "sd_Go", "se_Go", "ci_Go")
colnames(fig3Btab)[((ncol(fig3Btab) - 1):ncol(fig3Btab))] <- c("size_A", "size_Go")


p2 <- ggplot(fig3Btab, aes(x = TSS_Go, y = TSS_A, col = algorithm)) + 
  geom_point(aes(shape = vs_size, size = size_A))  +
  xlim(0.9, 1) + geom_errorbar(aes(ymin = TSS_A - ci_A, ymax = TSS_A + ci_A)) +
  geom_errorbarh(aes(xmin = TSS_Go - ci_Go, xmax = TSS_Go + ci_Go)) +
  ylim(0, 1) + theme_bw()


### Species raster
conds <- raster("./data/present.tiff.asc")
conds
sp2.big.hd <- raster("./data/big/2/present.tiff")

a1 <- results[which(results$vs_size == "big" & results$scenario == "HD"), ]

#### Conversion table ####
##     Corresp.       A          Go
# TP     a         A1 + A2       A1
# FP     b           A3       A2 + A3
# FN     c         A4 + A5       A4
# TN     d           A6       A5 + A6


results$sens_A <- (results$A1 + results$A2) / (results$A1 + results$A2 + results$A4 + results$A5)
results$spe_A <- (results$A6) / (results$A3 + results$A6)

# Values aren't strictly identical but very close...
a <- round(results$sensitivity_A, 3) - round(results $ sens_A, 3)
b <- round(results$specificity_A, 3) - round(results$spe_A, 3)
# % of different values
length(a[a != 0]) / length(a)
length(b[b != 0]) / length(a)
boxplot(a, b)


# Calculation of new indices
results$Jacc_A <- (results$A1 + results$A2) / (results$A1 + results$A2 + results$A3 + results$A4 + results$A5)
results$OPR_A <- (results$A3) / (results$A1 + results$A2 + results$A3)
results$UPR_A <- (results$A4 + results$A5) / (results$A1 + results$A2 + results$A4 + results$A5)

results$Jacc_Go <- (results$A1) / (results$A1 + results$A2 + results$A3 + results$A4)
results$OPR_Go <- (results$A2 + results$A3) / (results$A1 + results$A2 + results$A3)
results$UPR_Go <- (results$A4) / (results$A1 + results$A4)

results$Sor_A <- 2 * (results$A1 + results$A2) / (2 * (results$A1 + results$A2) + results$A3 + results$A4 + results$A5)
results$Sor_Go <- 2 * (results$A1) / (2 * results$A1 + results$A2 + results$A3 + results$A4)

fig.new.indA <- data.frame(summarySE(data = results,
                                    measurevar = "Jacc_A",
                                    groupvars = c("algorithm", "scenario")),
                          summarySE(data = results,
                                    measurevar = "Jacc_Go",
                                    groupvars = c("algorithm", "scenario"))[, c("Jacc_Go", "sd", "se", "ci")],
                          summarySE(data = results,
                                    measurevar = "size_A",
                                    groupvars = c("algorithm", "scenario"))$size_A,
                          summarySE(data = results,
                                    measurevar = "A1",
                                    groupvars = c("algorithm", "scenario"))$A1)
colnames(fig.new.indA)[colnames(fig.new.indA) %in% c("sd", "se", "ci", "sd.1", "se.1", "ci.1")] <- c("sd_A", "se_A", "ci_A", "sd_Go", "se_Go", "ci_Go")
colnames(fig.new.indA)[((ncol(fig.new.indA) - 1):ncol(fig.new.indA))] <- c("size_A", "size_Go")

p3 <- ggplot(fig.new.indA, aes(x = Jacc_Go, y = Jacc_A, col = algorithm)) + 
  geom_point(aes(shape = scenario, size = size_A)) +
  geom_errorbar(aes(ymin = Jacc_A - ci_A, ymax = Jacc_A + ci_A)) +
  geom_errorbarh(aes(xmin = Jacc_Go - ci_Go, xmax = Jacc_Go + ci_Go)) +
  ylim(0, 1) + xlim(0, 1) + theme_bw()



fig.new.indB <- data.frame(summarySE(data = results,
                                 measurevar = "Jacc_A",
                                 groupvars = c("algorithm", "vs_size")),
                       summarySE(data = results,
                                 measurevar = "Jacc_Go",
                                 groupvars = c("algorithm", "vs_size"))[, c("Jacc_Go", "sd", "se", "ci")],
                       summarySE(data = results,
                                 measurevar = "size_A",
                                 groupvars = c("algorithm", "vs_size"))$size_A,
                       summarySE(data = results,
                                 measurevar = "A1",
                                 groupvars = c("algorithm", "vs_size"))$A1)
colnames(fig.new.indB)[colnames(fig.new.indB) %in% c("sd", "se", "ci", "sd.1", "se.1", "ci.1")] <- c("sd_A", "se_A", "ci_A", "sd_Go", "se_Go", "ci_Go")
colnames(fig.new.indB)[((ncol(fig.new.indB) - 1):ncol(fig.new.indB))] <- c("size_A", "size_Go")

p4 <- ggplot(fig.new.indB, aes(x = Jacc_Go, y = Jacc_A, col = algorithm)) + 
  geom_point(aes(shape = vs_size, size = size_A)) +
  geom_errorbar(aes(ymin = Jacc_A - ci_A, ymax = Jacc_A + ci_A)) +
  geom_errorbarh(aes(xmin = Jacc_Go - ci_Go, xmax = Jacc_Go + ci_Go)) +
  xlim(0, 1) + ylim(0, 1) + theme_bw()

multiplot(p1, p3, p2, p4, cols = 2)

comp.tab <- data.frame(fig3Atab,
                       fig.new.indA)

p5 <- ggplot(comp.tab, aes(x = TSS_A, y = Jacc_A, col = algorithm)) +
  geom_point(aes(shape = scenario), size = 3)  +
  geom_errorbar(aes(ymin = Jacc_A - ci_A.1, ymax = Jacc_A + ci_A.1)) +
  geom_errorbarh(aes(xmin = TSS_A - ci_A, xmax = TSS_A + ci_A)) +
  xlim(0, 1) + ylim(0, 1)
  
results$total.area <- rowSums(results[, paste("A", 1:6, sep = "")])

results[, paste("A", 1:6, sep = "")] <- apply(results[, paste("A", 1:6, sep = "")] , 2, as.numeric)

results$Y <- rowSums(results[, c("A1", "A2", "A3")]) * rowSums(results[, c("A1", "A2", "A4", "A5")]) / results$total.area

reclass.Y <- function(x) {
  if (x < 0.001)
  { "#D8B365" } else
    if(x >= 0.001 & x <= 0.01)
    { "#F5F5F5" } else
      if(x > 0.1)
      { "#5AB4AC" }
}

results$ratio.TP.y <- results$Y / (results$A1 + results$A2)

results$class.Y <- NA
results$class.Y[which(results$ratio.TP.y < 0.001)] <- "#D8B365"
results$class.Y[which(results$ratio.TP.y >= 0.001 & results$ratio.TP.y <= 0.01)] <- "#F5F5F5"
results$class.Y[which(results$ratio.TP.y > 0.01)] <- "#5AB4AC"

results$kappa_A <- 2 * ((results$A1 + results$A2) * results$A6 - results$A3 * (results$A4 + results$A5)) /
  (results$A3^2 + results$A3 * results$A6 + (results$A4 + results$A5) * (results$A4 + results$A5 + results$A6) +
     (results$A1 + results$A2) * (results$A3 + results$A4 + results$A5 + 2 * results$A6))


ggplot(results, aes(x = kappa_A, y = Jacc_A)) +
  geom_point(aes(col = ratio.TP.y))

ggplot(results, aes(x = Kappa_A, y = Sor_A)) +
  geom_point(aes(col = class.Y))

ggplot(results, aes(x = Kappa_Go, y = Sor_Go)) +
  geom_point()

