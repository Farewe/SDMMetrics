library(raster)
library(virtualspecies)
load(file = "./Raster 1000 PA/species4.0.10")
load("./Raster 1000 PA/species4.0.01")
load("./Raster 1000 PA/species4.0.50")

class(species4.0.01) <- "virtualspecies"
class(species4.0.10) <- "virtualspecies"
class(species4.0.50) <- "virtualspecies"

# Need recheck méthode de calcul prévalence dans virtualspecies
species4.0.01 <- convertToPA(species4.0.01, species.prevalence = 0.015, alpha = -.025)
species4.0.10 <- convertToPA(species4.0.01, species.prevalence = 0.1, alpha = -.025)
species4.0.50 <- convertToPA(species4.0.01, species.prevalence = 0.5, alpha = -.025)

thres.selection <- array(dim = c(6, 12),
                         dimnames = list(c(rep(0.01, 2), rep(0.1, 2), rep(0.5, 2)),
                                         c("True.cutoff", "Cutoff", "TP", "TN", "FP", "FN", "Sensitivity", "Specificity", "TSS", "Jaccard", "OPR", "UP")))


for (cur.prev in c(0.01, 0.1, 0.5))
{
  cur.sp <- eval(parse(text = paste0("species4.", formatC(cur.prev, digits = 2, format = "f"))))
  
  # Penser à vérifier les arrondis dans virtualspecies... (cf. lignes suivantes)
  eval.points <- sampleOccurrences(cur.sp, n = 1201, type = "presence-absence", sample.prevalence = 200/1000)
  length(which(eval.points$sample.points$Real == 1))
  length(which(eval.points$sample.points$Real == 0))
  
  eval.points$sample.points$Suitability <- extract(cur.sp$suitab.raster, eval.points$sample.points[, c("x", "y")])
  
  thres.seq <- seq(0, 1, length = 101)
  results <- array(dim = c(length(thres.seq), 10),
                   dimnames = list(thres.seq, 
                                   c("TP", "TN", "FP", "FN", "Sensitivity", "Specificity", "TSS", "Jaccard", "OPR", "UP")))
  for (thresh in thres.seq)
  {
    tmp <- eval.points$sample.points
    tmp$prediction <- 0
    tmp$prediction[which(tmp$Suitability >= thresh)] <- 1
    TP <- length(which(tmp$Real == 1 & tmp$prediction == 1)) # TP
    TN <- length(which(tmp$Real == 0 & tmp$prediction == 0)) # TN
    FP <- length(which(tmp$Real == 0 & tmp$prediction == 1)) # FP
    FN <- length(which(tmp$Real == 1 & tmp$prediction == 0)) # FN
    results[which(rownames(results) == thresh), ] <-
      c(TP, 
        TN,
        FP,
        FN,
        TP / (TP + FN),
        TN / (TN + FP),
        TP / (TP + FN) + TN / (TN + FP) - 1,
        TP / (TP + FP + FN),
        FP / (TP + FP),
        FN / (TP + FN))
  }
  thres.selection[which(rownames(thres.selection) == cur.prev), ][1, ] <-
    c(as.numeric(cur.sp$PA.conversion["beta"]),
      as.numeric(rownames(results)[which(results[, "TSS"] == max(results[, "TSS"]))[1]]),
      results[which(results[, "TSS"] == max(results[, "TSS"]))[1], ])
  thres.selection[which(rownames(thres.selection) == cur.prev), ][2, ] <-
    c(as.numeric(cur.sp$PA.conversion["beta"]),
      as.numeric(rownames(results)[which(results[, "Jaccard"] == max(results[, "Jaccard"]))[1]]),
      results[which(results[, "Jaccard"] == max(results[, "Jaccard"]))[1], ])
}


thres.selection



