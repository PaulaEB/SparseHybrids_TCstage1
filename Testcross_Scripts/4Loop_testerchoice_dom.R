## Loop over dominance, tester levels, diagnostics to choose a tester with random sampling ##
# Here I just added some diagnostics proposed by DT to choose testers and ignore the phenotypic values simulation
# March 26, 2026, PaulaE
rm(list = ls())
# .rs.restartR()

library(AlphaSimR)
library(FieldSimR)
library(dplyr)
library(tibble)
# library(asreml)

foldertosave <- paste0("4Loop_results_", format(Sys.time(), "%d_%m_%H%M"))
dir.create(foldertosave, showWarnings = FALSE, recursive = TRUE)

# Testcross parameters
nparents <- 80 #Number of parents in both pools
nchr <- 10 #Number of chromosomes in maize
nqtl <- 500 #Lit causal QTL for a trait. So total is nchr x nqtl = 5000 QTL
nsnps <- 500 # Markers not necessarily in LD with QTL of with effects on the trait
nGenSplit <- c(2,50) #Heterotic pool split
# ncrosses <- 60

# Trait parameters
meanG <- 70
varG <- 20
h2 <- 0.3 # plot-level heritability, ratio of genetic to phenotypic variance, varG/(varG + varE)
meanDD<- c(1.2,0.5,0) #high, medium and low dominance according to Gonzalez-Dieguez 2025
# varDD <- 0.2

# Experimental design parameters: full factorial
nenvs <- 1
nblocks <- nreps <- 2
ngenos_ff <- (nparents/2)*(nparents/2)

# Experimental design parameters: testcross
nenvs <- 1
nblocks <- nreps <- 2

diagnostics_global <- list()
heterosis_ff <- list()
heterosis_tc <- list()
var_A_D <- list()
gca_sca_ff_tc_global <- list()
gca_sca_ff <- list()

# Loop over dominance levels #

t0 <- Sys.time()

for(s in nGenSplit){
for(d in meanDD){
  
  cat("Simulating founders with split of",s, "and meanDD of ", d, "\n")
  
  outdir <- file.path(foldertosave, paste0("split", s), paste0(d,"DD"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  founderpop <- runMacs(nInd=nparents, 
                        nChr=nchr, 
                        segSites= (nsnps+nqtl), #number of segregating sites per chromosome 1000 * 10= 10,000 sites
                        inbred=TRUE,
                        species="MAIZE",
                        split = s)
  
  SP <- SimParam$new(founderpop)
  SP$setTrackPed(TRUE)
  # SP$restrSegSites(minSnpPerChr = nsnps, overlap = F)
  SP$addSnpChip(nsnps)
  if (d == 0) {varDD <- 0} else {varDD <- 0.2}
  SP$addTraitAD(nQtlPerChr=nqtl, mean = meanG, var = varG, meanDD = d, varDD = varDD)
  
  #Simulation both pools only with split as divergence 
  
  poolA <- newPop(founderpop[1:(nparents/2)])
  poolB <- newPop(founderpop[((nparents/2)+1):(nparents)])
  
  ## Full factorial cross ##
  
  hybridpop <- hybridCross(females = poolA,
                           males = poolB, 
                           crossPlan = "testcross")
  
  ##Population structure and allele frequencies between pools##
  
  M_segsitesA <- pullSegSiteGeno(poolA)
  M_snpsA <- pullSnpGeno(poolA)
  M_qtlA <- pullQtlGeno(poolA)
  p_sitesA <- colMeans(M_segsitesA)/2
  p_snpsA<- colMeans(M_snpsA)/2
  p_qtlA <- colMeans(M_qtlA)/2
  
  M_segsitesB <- pullSegSiteGeno(poolB)
  M_snpsB <- pullSnpGeno(poolB)
  M_qtlB <- pullQtlGeno(poolB)
  p_sitesB <- colMeans(M_segsitesB)/2
  p_snpsB <- colMeans(M_snpsB)/2
  p_qtlB <- colMeans(M_qtlB)/2
  
  cat(" Overlap between QTL and SNPs in pool A: ", sum(colnames(M_qtlA) %in% colnames(M_snpsA)), "\n")
  
  delta_p_sites <- (p_sitesA - p_sitesB)
  delta_p_qtl <- (p_qtlA - p_qtlB)
  absdelta_p_sites <- abs(p_sitesA - p_sitesB)
  absdelta_p_qtl <- abs(p_qtlA - p_qtlB)
  
  png(file.path(outdir, paste0("1_split_",s,"DD_",d,"_hist_sitespools.png")), width = 2400, height = 1800, res = 600)
  hist(p_sitesA, breaks = 15, freq = FALSE,
       col = adjustcolor("purple2", 0.50), cex.main = 0.7,
       main = sprintf("All sites allfreq | split=%d | DD=%.2f\n absdelta=%.2f",
                      s, d, mean(absdelta_p_sites)))
  hist(p_sitesB, breaks = 15, add = TRUE, freq = FALSE,
       col = adjustcolor("grey20", 0.25))
  dev.off()
  
  png(file.path(outdir, paste0("1_split_",s,"DD_",d, "_hist_qtlpools.png")), width = 2400, height = 1800, res = 600)
  hist(p_qtlA, breaks = 15, freq = FALSE,
       col = adjustcolor("purple2", 0.50), cex.main = 0.7,
       main = sprintf("QTL allfreq | split=%d | DD=%.2f\n absdelta=%.2f",
                      s, d,mean(absdelta_p_qtl)))
  hist(p_qtlB, breaks = 15, add = TRUE, freq = FALSE,
       col = adjustcolor("grey20", 0.25))
  dev.off()
  
  # Allele frequencies correlation plots
  
  cat("Correlation between allele frequencies of all the sites", cor(p_sitesA, p_sitesB), "\n") 
  cor_allele_freqssites <- cor(p_sitesA, p_sitesB) 
  cor_allele_freqssnps <- cor(p_snpsA, p_snpsB) 
  cor_allele_freqsqtl <- cor(p_qtlA, p_qtlB) 
  
  png(file.path(outdir, paste0("1_split_",s,"_DD_", d, "_corallfreqssites.png")), width = 2400, height = 1800, res = 600)
  plot(p_sitesA, p_sitesB,
       xlab = "Allele frequencies Pool A",
       ylab = "Allele frequencies Pool B",
       main = paste0("All sites allele frequencies", "\n", "split ", s, " DD",d),
       pch = 1, col = "grey60")
  abline(0,1, lwd = 2)
  legend("topleft",
         legend = sprintf("correlation = %.3f", cor(p_sitesA, p_sitesB)),
         bty = "n")
  dev.off()

  png(file.path(outdir, paste0("1_split_",s,"_DD_", d, "_corallfreqssnps.png")), width = 2400, height = 1800, res = 600)
  plot(p_snpsA, p_snpsB,
       xlab = "Allele frequencies Pool A",
       ylab = "Allele frequencies Pool B",
       main = paste0("SNPs allele frequencies", "\n", "split ", s, " DD",d),
       pch = 1, col = "grey60")
  abline(0,1, lwd = 2)
  legend("topleft",
         legend = sprintf("correlation = %.3f", cor(p_snpsA, p_snpsB)),
         bty = "n")
  dev.off()
  
  png(file.path(outdir, paste0("1_split_",s,"_DD_", d, "_corallfreqsqtl.png")), width = 2400, height = 1800, res = 600)
  plot(p_qtlA, p_qtlB,
       xlab = "Allele frequencies Pool A",
       ylab = "Allele frequencies Pool B",
       main = paste0("QTL allele frequencies", "\n", "split ", s, " DD",d),
       pch = 1, col = "grey60")
  abline(0,1, lwd = 2)
  legend("topleft",
         legend = sprintf("correlation = %.3f", cor(p_qtlA, p_qtlB)),
         bty = "n")
  dev.off()
  
  #PCA plots between pools with SNPs
  
  G <- rbind(M_snpsA, M_snpsB)
  pools <- factor(c(rep("A", nrow(M_snpsA)), rep("B", nrow(M_snpsB))))
  pc <- prcomp(G, center = TRUE, scale. = FALSE)
  summary(pc)
  M_vaf1 <- 100*pc$sdev[1]^2/sum(pc$sdev^2)
  M_vaf2 <- 100*pc$sdev[2]^2/sum(pc$sdev^2)
  
  cent_A <- c(mean(pc$x[pools=="A", 1]), mean(pc$x[pools=="A", 2]))
  cent_B <- c(mean(pc$x[pools=="B", 1]), mean(pc$x[pools=="B", 2]))
  
  png(file.path(outdir, paste0("2_split_",s,"_DD_", d, "pca_pools.png")), width = 2800, height = 2200, res = 600)
  plot(pc$x[,1], pc$x[,2],
       col = c("A"="purple","B"="black")[pools], pch = 1, cex = 0.5,
       xlab = sprintf("PC1 (%.1f%%)", M_vaf1),
       ylab = sprintf("PC2 (%.1f%%)", M_vaf2), main = paste0("SNPs Between pools DD",d))
  points(cent_A[1], cent_A[2], pch = 4, cex = 1.1, lwd = 2, col = "purple4")
  points(cent_B[1], cent_B[2], pch = 4, cex = 1.1, lwd = 2, col = "black")
  
  legend("bottomleft", c("Pool A","Pool B"), col = c("purple","black"), pch = 19, bty = "n")
  dev.off()
  
  # PCA plot pools and hybrids with SNPs
  
  M_snpshybrids <- pullSnpGeno(hybridpop)
  G_all<-rbind(M_snpsA, M_snpsB, M_snpshybrids)
  grp <- factor(c(rep("A", nrow(M_snpsA)),
                  rep("B", nrow(M_snpsB)),
                  rep("H", nrow(M_snpshybrids))),
                levels = c("A","B","H"))
  
  pc_all <- prcomp(G_all, center = TRUE, scale. = FALSE)
  varExp <- pc_all$sdev^2 / sum(pc_all$sdev^2)
  
  colors <- c(H = "#FFC20A", A = "purple", B = "black")
  pchs <- c(H = 1, A = 19, B = 19)
  
  cent_H <- c(mean(pc_all$x[grp=="H", 1]), mean(pc_all$x[grp=="H", 2]))
  
  png(file.path(outdir, paste0("2_split_",s,"_DD_", d, "pca_all.png")), width = 2600, height = 2200, res = 600)
  plot(pc_all$x[grp=="H",1], pc_all$x[grp=="H",2],
       col = colors["H"], pch = pchs["H"], cex = 0.9, 
       xlab = sprintf("PC1 (%.1f%%)", 100*varExp[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100*varExp[2]),
       main = paste0("SNPs all pops split",s,"DD", d))
  
   points(pc_all$x[grp=="A",1], pc_all$x[grp=="A",2],
         col = colors["A"], pch = pchs["A"], cex = 1.0)
   points(pc_all$x[grp=="B",1], pc_all$x[grp=="B",2],
         col = colors["B"], pch = pchs["B"], cex = 1.0)
  
   points(cent_H[1], cent_H[2], pch = 4, cex = 1.5, lwd = 3, col = "red")
   
   legend("bottomleft",
         legend = c("Hybrids","Pool A","Pool B","Hybrid centroid"),
         col = c(colors["H"],colors["A"],colors["B"],"red"),
         pch = c(pchs["H"],pchs["A"],pchs["B"],4),
         pt.cex = c(0.9, 1.0, 1.0,0.5),
         bty = "n")
  dev.off()
  
  # Histograms of genetic values of hybrids and parents
  
  png(file.path(outdir, paste0("3_split_",s,"_DD_", d, "hist_gv.png")), width = 4000, height = 3800, res = 600)
  hist(hybridpop@gv[,1], breaks = 10, freq = FALSE,
       xlim = range(c(hybridpop@gv[,1], poolA@gv[,1], poolB@gv[,1])),
       ylim = c(0, max(density(hybridpop@gv[,1])$y, density(poolA@gv[,1])$y, density(poolB@gv[,1])$y) * 1.5),
       col = adjustcolor("#FFC20A", 0.30), border = NA, cex.main = 0.5,
       main = sprintf("GV split=%d | DD=%.2f\nHybrids mean=%.2f var=%.2f\nPoolA mean=%.2f var=%.2f\nPoolB mean=%.2f var=%.2f",
                      s, d,
                      mean(hybridpop@gv[,1]), popVar(hybridpop@gv[,1]),
                      mean(poolA@gv[,1]), popVar(poolA@gv[,1]),
                      mean(poolB@gv[,1]), popVar(poolB@gv[,1])))
  hist(poolA@gv[,1], breaks = 10, add = TRUE, freq = FALSE,
       col = adjustcolor("purple2", 0.25), border = NA)
   
  hist(poolB@gv[,1], breaks = 10, add = TRUE, freq = FALSE,
       col = adjustcolor("grey20", 0.25), border = NA)

  legend("topleft", bty = "n",
         legend = c("Hybrids", "Pool A", "Pool B"),
         fill = c("#FFC20A", "purple2", "grey20"),
         border = NA)

  abline(v = mean(hybridpop@gv[,1]), col = "#FFC20A",    lty = 1, lwd = 3)
  abline(v = mean(poolA@gv[,1]),     col = "purple2", lty = 1, lwd = 3)
  abline(v = mean(poolB@gv[,1]),     col = "grey20", lty = 1, lwd = 3)
  dev.off()
  
  ## Diagnostic functions for testers' choice
  # 
  # rep_score <- function(weights, T_idx) {
  #   pT <- colMeans(M_snpsB[T_idx, , drop = FALSE] / 2)
  #   pB_hat <- p_snpsB
  #   -sum(weights * (pT - pB_hat)^2)
  # }
  
  # redundancy <- function(T_idx){
  #   if(length(T_idx) < 2) return(0)
  #   G <- cor(t(XB[T_idx,,drop=FALSE]))
  #   # weighted correlation?
  #   mean(G[upper.tri(G)])
  # }
  
  # # Combine into an index (tune lambdas as you like)
  # tester_index <- function(T_idx, lam=c(alignment=1, rep=0.1, disc=0.1, red=0.1)){
  #   lam["alignment"]*gca_align(T_idx) + # not available in practice
  #     lam["rep"]*rep_score(weights = 1, T_idx) +
  #     lam["disc"]*log1p(disc_power(T_idx)) -
  #     lam["red"]*redundancy(T_idx)
  # }
   
  
  ## Heterosis and dominance levels and distance between populations with SNPs ##

  a = SP$traits[[1]]@addEff
  dom = SP$traits[[1]]@domEff
  int = SP$traits[[1]]@intercept
  q_qtlA <- 1 - p_qtlA
  q_qtlB <- 1 - p_qtlB
  
  y = p_qtlA - p_qtlB
  pbar = (p_qtlA + p_qtlB) / 2
  qbar = (q_qtlA + q_qtlB) / 2
  
  HF1 <- sum(dom* (y^2)) #Falconer's panmictic heterosis, according to deJong is panmictic
  BaseHF1 <- sum((2*pbar*qbar*dom) - ((1/2)*((y^2)*dom))) #deJong Lamkey 1999
  InbredHF1 <- sum((2*pbar*qbar*dom) + ((1/2)*((y^2)*dom))) #deJong Lamkey 1999
  MPH <- mean(hybridpop@gv[,1]) - ((mean(poolA@gv[,1]) + mean(poolB@gv[,1]))/2) #heterosis 
  BPH <- mean(hybridpop@gv[,1]) - max(mean(poolA@gv[,1]), mean(poolB@gv[,1]))
  
  # NeisDist <- mean((p_snpsA- p_snpsB)^2)
  # Fst <- NeisDist/(mean(p_snpsA*(1 - p_snpsB)) + mean(p_snpsB*(1 - p_snpsA))) 
  # NeisDist_qtl <- mean((p_qtlA - p_qtlB)^2)
  
  heterosis_ff[[as.character(s)]][[as.character(d)]] <- data.frame(split = s,
                                                                   meanDD = d, 
                                                                   cor_allele_freqssites = cor_allele_freqssites,
                                                                   cor_allele_freqssnps = cor_allele_freqssnps,
                                                                   cor_allele_freqsqtl = cor_allele_freqsqtl,
                                                                   delta_p_sites = mean(delta_p_sites),
                                                                   absdelta_p_sites = mean(absdelta_p_sites),
                                                                   delta_p_qtl = mean(delta_p_qtl),
                                                                   absdelta_p_qtl = mean(absdelta_p_qtl),
                                                                   HF1 = HF1, BaseHF1=BaseHF1, InbredHF1=InbredHF1,
                                                                   MPH = MPH, BPH = BPH)
  
  var_A_D[[as.character(s)]][[as.character(d)]] <- data.frame(split =s, 
                                                              meanDD = d, 
                                                              varA = varA(hybridpop), 
                                                              varD = varD(hybridpop), 
                                                              varAandD = varA(hybridpop) + varD(hybridpop),
                                                              varG = varG(hybridpop), 
                                                              varAovervarG = varA(hybridpop) / varG(hybridpop),
                                                              varDovervarG = varD(hybridpop) / varG(hybridpop))
  ## Genetic values of hybrids ##
  
  gv_ff_df <- data.frame(env = 1,
                         mother_poolA = factor(as.numeric(hybridpop@mother)), 
                         father_poolB = factor(as.numeric(hybridpop@father)), 
                         id = factor(as.numeric(hybridpop@id)),
                         rep = factor(rep(1:nreps, each = ngenos_ff)),
                         gv.TraitA = as.numeric(hybridpop@gv[,1]))
  
  gv_ff_df <- gv_ff_df[order(gv_ff_df$mother_poolA, gv_ff_df$father_poolB), ]
  
  # Create a matrix with the genetic values per combination of parents
  gv_ff_matrix <- with(gv_ff_df, tapply(gv.TraitA, list(mother_poolA, father_poolB), unique)) #matrix with the genetic values of the hybrids per combination of parents
  colnames(gv_ff_matrix) <- trimws(colnames(gv_ff_matrix)) #tidy column names, removing extra spaces
  pop_gv_mean_trait1 <- mean(gv_ff_matrix) #mean of the population in general, same as meanG = 70
  gv_ff_matrix <- gv_ff_matrix - pop_gv_mean_trait1 #genetic values without the intercept, so we can calculate the GCA effects directly
  
  # Calculate true GCA effects (all parent combos) for pool A trait 1
  gca_full_A_trait1 <- data.frame(mother_id = rownames(gv_ff_matrix),
                                  GCA = rowMeans(gv_ff_matrix))
  sigma_gA2 <- popVar(gca_full_A_trait1$GCA)

  png(file.path(outdir, paste0("4_split_",s,"_DD_", d, "histGCA.png")), width = 2600, height = 2000, res = 600)
  hist(gca_full_A_trait1$GCA, breaks = 20,cex.main = 0.8,
       main = sprintf("GCA effects Pool A | split=%d | DD=%.2f\nmean=%.2f var=%.2f",
                      s, d, mean(gca_full_A_trait1$GCA), sigma_gA2))
  abline(v = mean(gca_full_A_trait1$GCA), col = "purple2", lty = 1, lwd = 3)
  dev.off()

  gca_full_B_trait1 <- data.frame(father_id = colnames(gv_ff_matrix),
                                  GCA = colMeans(gv_ff_matrix))
  sigma_gB2 <- popVar(gca_full_B_trait1$GCA)
  
  sca_table_trait1 <- sweep(gv_ff_matrix, 1, gca_full_A_trait1$GCA, FUN = "-") 
  sca_table_trait1 <- sweep(sca_table_trait1, 2, gca_full_B_trait1$GCA, FUN = "-") #calculate the SCA effects as the deviation from the GCA effects of both parents
  
  sca_full_trait1 <- data.frame(mother_id = rep(rownames(gv_ff_matrix), each = ncol(gv_ff_matrix)),
                                father_id = as.numeric(rep(colnames(gv_ff_matrix), times = nrow(gv_ff_matrix))),
                                SCA = as.vector(sca_table_trait1))
  
  sigma_sca2 <- popVar(sca_full_trait1$SCA)

  png(file.path(outdir, paste0("4_split_",s,"_DD_", d, "histSCA.png")), width = 2600, height = 2000, res = 600)
  hist(sca_full_trait1$SCA, breaks = 20,cex.main = 0.8,
       main = sprintf("SCA effects | split=%d | DD=%.2f\nmean=%.2f var=%.2f",
                      s, d, mean(sca_full_trait1$SCA), sigma_sca2))
  abline(v = mean(sca_full_trait1$SCA), col = "#FFC20A", lty = 1, lwd = 3)
  dev.off()

  gca_sca_ff[[as.character(s)]][[as.character(d)]] <- data.frame(split = s,
                                                                   meanDD = d,
                                                                   mean_GCA_A = mean(gca_full_A_trait1$GCA),
                                                                   var_GCA_A = popVar(gca_full_A_trait1$GCA),
                                                                   mean_GCA_B = mean(gca_full_B_trait1$GCA),
                                                                   var_GCA_B = popVar(gca_full_B_trait1$GCA),
                                                                   mean_SCA   = mean(sca_full_trait1$SCA),
                                                                   var_SCA   = popVar(sca_full_trait1$SCA))
  
  png(file.path(outdir, paste0("4_split_",s,"_DD_", d, "barvarGCAsSCA.png")), width = 2100, height = 3200, res = 600)
  text(barplot(c(sigma_gA2,sigma_gB2,sigma_sca2), 
               names.arg = c("GCA_A", "GCA_B", "SCA"),
               col = c("purple2", "grey10", "#FFC20A"),
               main = paste0("Variance GCAs and SCA dd",d),
            ylab = "Variance", 
            ylim = c(0, max(c(sigma_gA2,sigma_gB2,sigma_sca2)) + 1),
            las = 1),
       c(sigma_gA2,sigma_gB2,sigma_sca2) + 0.5,
       labels = paste0("Var: ", round(c(sigma_gA2,sigma_gB2,sigma_sca2), 2),
                       "\nMean: ", round(c(mean(gca_full_A_trait1$GCA), mean(gca_full_B_trait1$GCA), mean(sca_full_trait1$SCA)), 6)),
       font = 1, cex = 0.5)
  dev.off()
  
  
 
  ## Loop over testers ##
  ntesters <- c(1,2,3,5)
  nsamples <- length(poolB@id)
  max_samples <- 5000
  top_n <- 10
  top_true <- as.numeric(gca_full_A_trait1$mother_id[order(gca_full_A_trait1$GCA, decreasing = TRUE)][1:top_n])
  
  diagnostics_testerset_list <- list()
  gca_sca_ff_tc_list <- list()
  het_gv_tc_list <- list()
  idx_testers_list <- list()
  
  ctr <- 0
  
    for(k in ntesters){
      niter <- if (k == 1) {
        3 * nsamples
      } else if (k == 2) {
        choose(nsamples, 2)
      }
      else {
        min(max_samples, ceiling(0.10 * choose(nsamples, k)))
      }
      
    for(i in 1:niter){
      ctr <- ctr + 1
      cat("N testers:",k, "it:",i,"\n")
      
      # Subset of testers from poolB
        if (k == 1) {
          testers_idx <- poolB@id[((i - 1) %% length(poolB@id)) + 1]
        } else {
          testers_idx <- sample(poolB@id, k, replace = FALSE)
        }
      
      tester_set <- paste(sort(testers_idx), collapse = ",")
      
      gca_acc_expected <- sqrt(sigma_gA2 / (sigma_gA2 + sigma_sca2 / k))
      sca_acc_expected <- sqrt(1 - 1 / k)
      var_gca_expected <- sigma_gA2 + sigma_sca2 / k
      var_sca_expected <- sigma_sca2 * (1 - 1/ k)
      
      # Correlation between the true GCA females from the full factorial cross and the GCA from the 3 testers
      gv_tc_matrix <- gv_ff_matrix[, as.character(testers_idx), drop = FALSE]
      mean(gv_tc_matrix)
      gca_testcross_A_trait1 <- data.frame(mother_id = rownames(gv_tc_matrix),
                                           GCA_tc = rowMeans(gv_tc_matrix))
      
      # gca_testcross_A_trait1 <- data.frame(mother_id = rownames(gv_ff_matrix),
      #                                      GCA_tc = rowMeans(cbind(gv_ff_matrix[, as.character(testers_idx)])))
      
      gca_testcross_B_trait1 <- data.frame(father_id = colnames(gv_tc_matrix),
                                      GCA_tc = NA)
      
      # gca_testcross_B_trait1 <- data.frame(father_id = as.character(testers_idx),
      #                                      GCA_tc = NA)
      
      sca_testcross_trait1 <- data.frame(mother_id = rep(rownames(gv_tc_matrix), each = ncol(gv_tc_matrix)),
                                    father_id = as.numeric(rep(colnames(gv_tc_matrix), times = nrow(gv_tc_matrix))),
                                    SCA_tc    = NA)
      
      # sca_testcross_trait1 <- data.frame(mother_id = rep(rownames(sca_table_trait1), times = length(testers_idx)),
      #                                    father_id = rep(testers_idx, each = nrow(sca_table_trait1)),
      #                                    SCA_tc    = NA)
       
      if (k > 1) {
      
      gca_testcross_B_trait1 <- data.frame(father_id = colnames(gv_tc_matrix),
                                             GCA_tc = colMeans(gv_tc_matrix))
          
      # gca_testcross_B_trait1 <- data.frame(father_id = as.character(testers_idx),
      #                                      GCA_tc    = colMeans(gv_ff_matrix[, as.character(testers_idx), drop = FALSE]))
      
      sca_tc_table_trait1 <- sweep(gv_tc_matrix, 1, gca_testcross_A_trait1$GCA_tc, FUN = "-") 
      sca_tc_table_trait1 <- sweep(sca_tc_table_trait1, 2, gca_testcross_B_trait1$GCA_tc, FUN = "-")
      
      sca_testcross_trait1 <- data.frame(mother_id = rep(rownames(sca_tc_table_trait1), times = length(testers_idx)),
                                           father_id = rep(testers_idx, each = nrow(sca_tc_table_trait1)),
                                           SCA_tc    = as.vector(sca_tc_table_trait1[, as.character(testers_idx)])) %>%
        arrange(father_id)
        
      }
      
      gv_ff_df_k_testers <- droplevels(gv_ff_df[gv_ff_df$father_poolB %in% testers_idx, ])
      
      MPH_gv_tc <- mean(gv_ff_df_k_testers$gv.TraitA) - (mean(poolA@gv[,1]) + mean(poolB@gv[,1]))/2
      BPH_gv_tc <- mean(gv_ff_df_k_testers$gv.TraitA) - max(mean(poolA@gv[,1]), mean(poolB@gv[,1]))
      
      idx_testers_list[[ctr]] <- data.frame(split = s,
                                            meanDD = d,
                                            ntesters = k,
                                            iter = i,
                                            tester_set = tester_set)
      
      gca_sca_ff_tc_list[[ctr]] <- data.frame(split = s,
                                              meanDD = d,
                                              ntesters = k,
                                              iter = i,
                                              tester_set = tester_set,
                                              gca_acc_expected = gca_acc_expected,
                                              var_gca_expected = var_gca_expected,
                                              sca_acc_expected = sca_acc_expected,
                                              var_sca_expected = var_sca_expected,
                                              mean_gcaA_full = mean(gca_full_A_trait1$GCA),
                                              var_gcaA_full = popVar(gca_full_A_trait1$GCA),
                                              mean_gcaB_full = mean(gca_full_B_trait1$GCA),
                                              var_gcaB_full = popVar(gca_full_B_trait1$GCA),
                                              mean_sca_full = mean(sca_full_trait1$SCA),
                                              var_sca_full = popVar(sca_full_trait1$SCA),
                                              mean_gcaA_tc   = mean(gca_testcross_A_trait1$GCA_tc),
                                              disc_power   = popVar(gca_testcross_A_trait1$GCA_tc),
                                              mean_gcaB_tc   = if (k > 1) mean(gca_testcross_B_trait1$GCA_tc) else NA,
                                              var_gcaB_tc = if (k > 1) popVar(gca_testcross_B_trait1$GCA_tc) else NA,
                                              mean_sca_tc   = if (k > 1) mean(sca_testcross_trait1$SCA_tc) else NA,
                                              var_sca_tc  = if (k > 1) popVar(sca_testcross_trait1$SCA_tc) else NA,
                                              gca_align = cor(gca_full_A_trait1$GCA, gca_testcross_A_trait1$GCA_tc),
                                              sca_align = if (k > 1) cor(as.vector(sca_table_trait1[, as.character(testers_idx), drop = FALSE]),as.vector(sca_tc_table_trait1)) else NA)
                                           
      het_gv_tc_list[[ctr]] <- data.frame(split = s,
                                          meanDD = d,
                                          ntesters = k,
                                          iter = i,
                                          tester_set = tester_set,
                                          MPH_tc = MPH_gv_tc,BPH_tc = BPH_gv_tc)
      
      
      #Diagnostics DT
      # Alignment
      gca_align <- cor(gca_full_A_trait1$GCA, gca_testcross_A_trait1$GCA_tc)
      
      # Discrimination among A using selected testers: Var of proxy g (bigger = more separation)
      disc_power <- popVar(gca_testcross_A_trait1$GCA_tc) 
      
      
      # Representativeness this would be to the mean within-pool:
      # weighted (not yet) squared distance between allele freqs in T and pool B (lower is better) 
      T_idx <- match(testers_idx, poolB@id) 
      representativeness <- sum((p_snpsB - colMeans(M_snpsB[T_idx, , drop = FALSE] / 2))^2)
      
      # Redundancy: average pairwise relatedness among chosen testers (lower is better)
      # Here: mean pairwise correlation across sites
      redundancy <- if (k < 2) {0} else {
        G <- cor(t(M_snpsB[T_idx, , drop = FALSE]))
        mean(G[upper.tri(G)])
      }
      
     diagnostics_testerset_list[[ctr]] <- data.frame(split = s,
                                            meanDD = d,
                                            ntesters = k,
                                            iter = i,
                                            tester_set = tester_set,
                                            gca_align = gca_align,
                                            disc_power = disc_power,
                                            representativeness = representativeness,
                                            redundancy = redundancy)
      } 
    }
  

  
  # Summaries/plots
  
  idx_testers_df <- dplyr::bind_rows(idx_testers_list)
  write.csv(idx_testers_df, file.path(outdir, paste0("idx_testers_df_s", s, "_DD", d, ".csv")), row.names = FALSE)
  gca_sca_ff_tc_df <- dplyr::bind_rows(gca_sca_ff_tc_list)
  gca_sca_ff_tc_df$ntesters <- factor(gca_sca_ff_tc_df$ntesters, levels = c(1,2,3,5))
  gca_sca_ff_tc_global[[as.character(s)]][[as.character(d)]] <- gca_sca_ff_tc_df
  write.csv(gca_sca_ff_tc_df, file.path(outdir, paste0("gca_sca_ff_tc_s", s, "_DD", d, ".csv")), row.names = FALSE)
  het_gv_tc_df <- dplyr::bind_rows(het_gv_tc_list)
  heterosis_tc[[as.character(s)]][[as.character(d)]] <- het_gv_tc_df
  diagnostics_testerset_df <- dplyr::bind_rows(diagnostics_testerset_list)
  diagnostics_global[[as.character(s)]][[as.character(d)]] <- diagnostics_testerset_df
  write.csv(diagnostics_testerset_df, file.path(outdir, paste0("diagnostics_testerchoice_s",s,"DD",d,".csv")), row.names = FALSE)
 
  png(file.path(outdir, paste0("5_split_",s,"_DD_", d, "GCA_SCA_expectedaccuracy.png")), width = 3100, height = 2300, res = 600)
  plot(gca_sca_ff_tc_df$ntesters,gca_sca_ff_tc_df$gca_acc_expected, type="b", pch = 16, cex=0.8, lwd = 2, col = "purple2",
       xlab = "Number of testers", ylab = "Expected accuracy", ylim = c(0, max(c(gca_sca_ff_tc_df$gca_acc_expected, gca_sca_ff_tc_df$sca_acc_expected), na.rm = TRUE) + 1),
                main = paste0("GCA expected accuracy | split=", s, " DD=", d),cex.main = 0.9,cex.lab = 0.8)
  lines(gca_sca_ff_tc_df$ntesters,gca_sca_ff_tc_df$sca_acc_expected, type="b", pch = 16, cex=0.8, lwd = 2, col = "#FFC20A")
  legend("topleft", legend = c(expression(GCA ==  sqrt(sigma[gA]^2 / (sigma[gA]^2 + sigma[sca]^2 / k))), expression(SCA == sqrt(1 - 1/k))),
         col = c("purple2", "#FFC20A"), pch = 16, bty = "n", cex = 0.5)
  dev.off()
  
  png(file.path(outdir, paste0("5_split_",s,"_DD_", d, "GCA_SCA_expectedvariance.png")), width = 3100, height = 2300, res = 600)
  plot(gca_sca_ff_tc_df$ntesters,gca_sca_ff_tc_df$var_gca_expected, type="b", pch = 16, cex=0.8, lwd = 2, col = "purple2",
       xlab = "Number of testers", ylab = "Expected variance", ylim = c(0, max(c(gca_sca_ff_tc_df$var_gca_expected, gca_sca_ff_tc_df$var_sca_expected), na.rm = TRUE) + 1),
       main = paste0("Expected variance GCA and SCA | split=", s, " DD=", d),cex.main = 0.9,cex.lab = 0.8)
  lines(gca_sca_ff_tc_df$ntesters,gca_sca_ff_tc_df$var_sca_expected, type="b", pch = 16, cex=0.8, lwd = 2, col = "#FFC20A")
  legend("topleft", legend = c(expression(GCA == sigma[gA]^2 + sigma[sca]^2 / k), expression(SCA == sigma[sca]^2 * (1 - 1/k))),
         col = c("purple2", "#FFC20A"), pch = 16, bty = "n", cex = 0.5)
  dev.off()
  
  png(file.path(outdir, paste0("6_split_",s,"_DD_", d, "GCA_alignment.png")), width = 2800, height = 3000, res = 600)
  boxplot(gca_align ~ ntesters, data = gca_sca_ff_tc_df, col = adjustcolor("#DA88A7", alpha.f = 0.65), border = "grey20",
          xlab = "Number of testers", ylab = "Correlation GCA full factorial and GCA testcross",
          ylimit = c(min(gca_sca_ff_tc_df$gca_align),1),
          main = paste0("GCA accuracy","\n","Divergence: ", s, " generations | Dominance: ", d), cex.main = 0.9, cex.lab = 0.8)
  points(x = as.numeric(gca_sca_ff_tc_df$ntesters), y = gca_sca_ff_tc_df$gca_acc_expected, col = "#DC3220", pch = 16, cex = 0.7)
  legend("topright", legend = "GCA expected accuracy testcross",
         col = "#DC3220", pch = 16, bty = "n", cex = 0.5)
  dev.off()
  
  png(file.path(outdir, paste0("6_split_",s,"_DD_", d, "GCA_discpower.png")), width = 2800, height = 3000, res = 600)
  boxplot(disc_power ~ ntesters, data = diagnostics_testerset_df,
          col = adjustcolor("#0C66C1", alpha.f = 0.50), border = "grey20",
          xlab = "Number of testers",
          ylab = "GCA variance of testcross",
          ylim = c(min(diagnostics_testerset_df$disc_power, na.rm = TRUE),
                   max(diagnostics_testerset_df$disc_power, na.rm = TRUE)),
          main = paste0("GCA variance Pool A","\n","Divergence: ", s, " generations | Dominance: ", d),
          cex.main = 0.9, cex.lab = 0.8)
  points(x = as.numeric(gca_sca_ff_tc_df$ntesters), y = gca_sca_ff_tc_df$var_gcaA_full, col = "#053B71", pch = 4, cex = 0.7)
  points(x = as.numeric(gca_sca_ff_tc_df$ntesters), y = gca_sca_ff_tc_df$var_gca_expected, col = "#053B71", pch = 16, cex = 0.7)
  legend("topright", legend = c("GCA expected variance testcross", "GCA variance full factorial"),
         col = c("#053B71", "#053B71"), pch = c(16, 4), bty = "n", cex = 0.5)
  dev.off()
  

} }

t1 <- Sys.time()  
runtime <- difftime(t1, t0, units = "mins")

summary_gca_sca_ff <- dplyr::bind_rows(lapply(gca_sca_ff, dplyr::bind_rows)) %>% dplyr::arrange(split, meanDD)
summary_var_A_D <- dplyr::bind_rows(lapply(var_A_D, dplyr::bind_rows)) %>% dplyr::arrange(split, meanDD)
summary_var_A_D <- summary_var_A_D %>% dplyr::rename(varA = Trait1,
                                 varD = Trait1.1,
                                 varAandD = Trait1.2,
                                 varG = Trait1.3,
                                 varAovervarG = Trait1.4,  
                                 varDovervarG = Trait1.5)
summary_heterosis_ff <- dplyr::bind_rows(lapply(heterosis_ff, dplyr::bind_rows)) %>% dplyr::arrange(split, meanDD)
summary_heterosis_tc <- dplyr::bind_rows(lapply(heterosis_tc, dplyr::bind_rows)) %>% dplyr::arrange(split, meanDD, ntesters, iter)
heterosis_ff_tc <- dplyr::left_join(summary_heterosis_ff, summary_heterosis_tc, by = c("split", "meanDD"))
summary_diagnostics <- dplyr::bind_rows(lapply(diagnostics_global, dplyr::bind_rows)) %>% dplyr::arrange(split, meanDD, ntesters, iter)
summary_gca_sca_global <- dplyr::bind_rows(lapply(gca_sca_ff_tc_global, dplyr::bind_rows)) %>% dplyr::arrange(split, meanDD, ntesters, iter)

write.csv(summary_gca_sca_global, file.path(foldertosave, "gca_sca_ff_tc_global.csv"), row.names = FALSE)
write.csv(summary_diagnostics, file.path(foldertosave, "diagnostics_global.csv"), row.names = FALSE)
write.csv(heterosis_ff_tc,file.path(foldertosave, "heterosis_ff_tc.csv"), row.names = FALSE)
write.csv(summary_gca_sca_ff, file.path(foldertosave, "gca_scaff.csv"), row.names = FALSE)
write.csv(summary_var_A_D, file.path(foldertosave, "var_ad_dom.csv"), row.names = FALSE)

for(ss in sort(unique(summary_var_A_D$split))){
  tmp <- subset(summary_var_A_D, split == ss)
  
  png(file.path(foldertosave, paste0("ratiovar_A_D_split_", ss, ".png")),
      width = 2100, height = 1500, res = 300)
  matplot(tmp$meanDD,
          as.matrix(tmp[, c("varAovervarG","varDovervarG")]),
          type = "b", pch = 16, lwd = 3, lty = 1,
          col = c("royalblue","orange3"),
          xlab = "meanDD", ylab = "Fraction",
          ylim = c(0,1),
          main = paste0("Genetic variance fractions in hybrids vs dominance | split=", ss))
  legend("topright", legend = c("varA/varG","varD/varG"),
         col = c("royalblue","orange3"),
         lwd = 3, pch = 16, bty = "n")
  dev.off()
}

# Faceted (side-by-side) matplot for split levels
splits_to_plot <- sort(unique(summary_var_A_D$split))

png(file.path(foldertosave, "ratiovar_A_D_global.png"),
    width = 9000, height = 5000, res = 600)

par(mfrow = c(1, length(splits_to_plot)),
    mar = c(4, 5.0, 1.2, 0.5) + 0.1,   # bottom, left, top, right
    oma = c(1, 0, 0.5, 1.5),
     yaxs = "i")

par(mgp = c(3.5, 1.2, 0))

for(j in seq_along(splits_to_plot)){
  ss <- splits_to_plot[j]
  tmp <- subset(summary_var_A_D, split == ss)
  
  matplot(tmp$meanDD,
          as.matrix(tmp[, c("varAovervarG","varDovervarG")]),
          type = "b", pch = 16, lwd = 3, lty = 1,
          col = c("royalblue","orange3"),
          xlab = "",
          ylab = if (j == 1) "Ratio" else "",
          cex.axis = 1.8, cex.lab = 2, cex.main = 2,
          xaxt = "n",
          ylim = c(-0.02, 1.02),
          main = paste0("Divergence ", ss, " generations"),
          bty = "n")   # removes top/right box lines
  
  axis(1, at = tmp$meanDD, labels = tmp$meanDD, cex.axis = 1.8)
  
  if(j == 2){
    legend("topright", legend = c("varA/varG","varD/varG"),
           col = c("royalblue","orange3"),
           lwd = 3, pch = 16, bty = "n", cex = 1.8)
  }
}

mtext("Dominance level", side = 1, outer = TRUE, line = 0.2, cex = 2)

dev.off()

for(ss in sort(unique(summary_gca_sca_ff$split))){
  tmp <- subset(summary_gca_sca_ff, split == ss)
  M <- t(as.matrix(tmp[, c("var_GCA_A","var_GCA_B","var_SCA")]))
  colnames(M) <- tmp$meanDD
  rownames(M) <- c("GCA_A","GCA_B","SCA")
  ratio <- tmp$var_SCA / (tmp$var_GCA_A + tmp$var_GCA_B)
  
  png(file.path(foldertosave, paste0("FF_barsGCAs_SCA_split", ss, ".png")),
      width = 2600, height = 1600, res = 300)
  par(mar = c(8, 4, 4, 2) + 0.1)
  
  text(x = barplot(M, beside = TRUE,
                   col = c("#0C66C1", "#C71D0B", "#AAA9BF"),
                   ylab = "Variance", xlab = "Dominance level",
                   # ylim = c(0, max(M, na.rm = TRUE) + 3),
                   ylim = c(0, 20),
                   main = paste0("Full factorial GCA and SCA variance components","\n","Divergence ", ss, " generations"),
                   las = 1),
       y = as.vector(M),
       labels = round(as.vector(M), 2), pos = 3, cex = 0.8)
  
  legend("topleft", legend = rownames(M),
         fill = c("#0C66C1","#C71D0B","#AAA9BF"), bty = "n")
  
  mtext(side = 1, line = 6,
        text = paste0("SCA to GCA variance ratio:", paste(round(ratio, 2), collapse = " | ")),
        cex = 0.9, font = 2)
  
  dev.off()
}

# Faceted (side-by-side) barplots for split 2 and 50
splits_to_plot <- sort(unique(summary_gca_sca_ff$split))

png(file.path(foldertosave, "FF_barsGCAs_SCA_global.png"), width = 9000, height = 5000, res = 600)

par(mfrow = c(1, length(splits_to_plot)),
    mar = c(4.2, 5.0, 1.8, 0.3) + 0.1,  # bottom, left, top, right (reduced)
    oma = c(1.0, 1.0, 0.6, 0),            # outer bottom/top reduced
    xaxs = "i", yaxs = "i")             # remove axis expansion padding
par(mgp = c(3.5, 1.2, 0))  # (title, labels, line)

for(j in seq_along(splits_to_plot)){
  ss <- splits_to_plot[j]
  tmp <- subset(summary_gca_sca_ff, split == ss)
  
  M <- t(as.matrix(tmp[, c("var_GCA_A","var_GCA_B","var_SCA")]))
  colnames(M) <- tmp$meanDD
  rownames(M) <- c("GCA_A","GCA_B","SCA")
  
  ratio <- tmp$var_SCA / (tmp$var_GCA_A + tmp$var_GCA_B)
  
  bp <- barplot(M, beside = TRUE,
                col = c("#0C66C1", "#C71D0B", "#AAA9BF"),
                xlab = "",                             # <- no per-panel x label
                ylim = c(0, 20),
                main = paste0("Divergence ", ss, " generations"),
                las = 1, cex.axis = 1.8, cex.lab = 2, cex.names = 1.8, cex.main = 2.2,
                ylab = if (j == 1) "Genetic variance components" else "",
                yaxt = if (j == 1) "s" else "n")
  
  text(x = bp, y = as.vector(M),
       labels = round(as.vector(M), 2), pos = 3, cex = 1.4)
  
  if(j == 1){
    legend("topleft", legend = rownames(M),
           fill = c("#0C66C1","#C71D0B","#AAA9BF"),
           bty = "n", cex = 1.8)
  }
  # 
  # mtext(side = 1, line = 6,
  #       text = paste0("SCA:(GCA_A + GCA_B) = ",
  #                     paste(round(ratio, 2), collapse = " | ")),
  #       cex = 0.9, font = 2)
}

# Shared X label centered across both panels
mtext("Dominance level", side = 1, outer = TRUE, line = 0.1, cex = 2)

dev.off()


summary_gca_sca_global$ntesters <- factor(summary_gca_sca_global$ntesters, levels = c(1,2,3,5))

align_sum <- summary_gca_sca_global %>%
  group_by(split, meanDD, ntesters) %>%
  summarise(align_med = median(gca_align, na.rm = TRUE),
            align_q25 = quantile(gca_align, 0.25, na.rm = TRUE),
            align_q75 = quantile(gca_align, 0.75, na.rm = TRUE),
            exp_acc   = median(gca_acc_expected, na.rm = TRUE),
            .groups = "drop") %>% arrange(split, meanDD, ntesters)

scen_cols <- c("split2_DD1.2"  = "#C51B7D",  
               "split2_DD0.5"  = "#F57C00",  
               "split2_DD0"    = "#FDB863", 
               "split50_DD1.2" = "#1B9E77", 
               "split50_DD0.5" = "#00A6D6", 
               "split50_DD0"   = "#7FD3E6")

scen_labs <- c("split2_DD1.2"  = "Div. 2 | Dom. mean 1.2",
               "split2_DD0.5"  = "Div. 2 | Dom. mean 0.5",
               "split2_DD0"    = "Div. 2 | Dom. mean 0",
               "split50_DD1.2" = "Div. 50 | Dom. mean 1.2",
               "split50_DD0.5" = "Div. 50 | Dom. mean 0.5",
               "split50_DD0"   = "Div. 50 | Dom. mean 0")

png(file.path(foldertosave, "GCA_accuracy_global.png"), width = 1600, height = 1500,res=450)
par(mar = c(3.2, 3, 0.5, 0.2) + 0.1)
par(mgp = c(1.5, 0.25, 0))  # (title, labels, line)

plot(NA, xlim = c(1,5), ylim = c(0.5,1), 
     xlab = "Number of testers", ylab = "Correlation GCA full factorial and GCA testcross", cex.lab = 0.7, cex.axis = 0.7)

for(ss in sort(unique(align_sum$split))){
  for(dd in sort(unique(align_sum$meanDD))){
    tmp <- subset(align_sum, split == ss & meanDD == dd)
    if(nrow(tmp) == 0) next
    scn <- paste0("split", ss, "_DD", dd)   # <-- define scenario key
    
    lines(as.numeric(as.character(tmp$ntesters)), tmp$align_med,
          col = scen_cols[scn], lty = 1, lwd = 1.5)
    
    arrows(x0 = as.numeric(as.character(tmp$ntesters)), y0 = tmp$align_q25,
           x1 = as.numeric(as.character(tmp$ntesters)), y1 = tmp$align_q75,
           angle = 90, code = 3, length = 0.02,
           col = adjustcolor(scen_cols[scn], 0.55), lwd = 1)
    
    # points(as.numeric(as.character(tmp$ntesters)), tmp$exp_acc,
    #        pch = ifelse(as.character(ss) == "2", 1, 16),
    #        cex = 0.5,
    #        col = split_cols[as.character(ss)],
    #        lwd = ifelse(as.character(ss) == "2", 1.2, NA))
  }
}

legend("bottomright", bty = "n", cex = 0.4,
       legend = unname(scen_labs[names(scen_cols)]),
       col = scen_cols, lty = 1, lwd = 2)

dev.off()

summary_diagnostics$ntesters <- factor(summary_diagnostics$ntesters, levels = c(1,2,3,5))

disc_sum <- summary_diagnostics %>%
  group_by(split, meanDD, ntesters) %>%
  summarise(
    disc_med = median(disc_power, na.rm = TRUE),
    disc_q25 = quantile(disc_power, 0.25, na.rm = TRUE),
    disc_q75 = quantile(disc_power, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(split, meanDD, ntesters)

gca_all <- dplyr::bind_rows(lapply( gca_sca_ff_tc_global, dplyr::bind_rows)) %>%
  dplyr::arrange(split, meanDD, ntesters, iter)

gca_all$ntesters <- factor(gca_all$ntesters, levels = c(1,2,3,5))

refs_sum <- gca_all %>%
  group_by(split, meanDD, ntesters) %>%
  summarise(
    var_gca_expected_med = median(var_gca_expected, na.rm = TRUE),
    var_gcaA_full_med    = median(var_gcaA_full, na.rm = TRUE),
    .groups = "drop"
  )

disc_sum <- left_join(disc_sum, refs_sum, by = c("split","meanDD","ntesters")) %>%
  arrange(split, meanDD, ntesters)

# --- Same as your previous working code: scenario colours + ALL lines continuous (no linetypes) ---

png(file.path(foldertosave, "Discpower_global.png"), width = 1600, height = 1500,res=450)
par(mar = c(3.2, 3, 0.5, 0.2) + 0.1)
par(mgp = c(1.5, 0.25, 0))  # (title, labels, line)

ymax <- max(disc_sum$disc_q75, disc_sum$var_gca_expected_med, disc_sum$var_gcaA_full_med, na.rm = TRUE)

plot(NA, xlim = c(1,5), ylim = c(4, ymax),
     xlab = "Number of testers",
     ylab = "Variance of GCA in Pool A", cex.lab = 0.7, cex.axis = 0.7)

for(ss in sort(unique(disc_sum$split))){
  for(dd in sort(unique(disc_sum$meanDD))){
    tmp <- subset(disc_sum, split == ss & meanDD == dd)
    if(nrow(tmp) == 0) next
    
    scn <- paste0("split", ss, "_DD", dd)  # scenario key
    
    # continuous lines, scenario colour
    lines(as.numeric(as.character(tmp$ntesters)), tmp$disc_med,
          col = scen_cols[scn],
          lty = 1, lwd = 1.5)
    
    arrows(x0 = as.numeric(as.character(tmp$ntesters)), y0 = tmp$disc_q25,
           x1 = as.numeric(as.character(tmp$ntesters)), y1 = tmp$disc_q75,
           angle = 90, code = 3, length = 0.02,
           col = adjustcolor(scen_cols[scn], 0.55), lwd = 1)
    
    # FULL FACTORIAL points only (same symbol logic you had)
    points(as.numeric(as.character(tmp$ntesters)), tmp$var_gcaA_full_med,
           pch = ifelse(as.character(dd) == "0", 4,                      # dominance 0 -> x
                        ifelse(as.character(ss) == "2", 1, 16)),          # split2 open, split50 filled
           cex = 0.5,
           col = scen_cols[scn],
           lwd = ifelse(as.character(dd) == "0", 1,
                        ifelse(as.character(ss) == "2", 1.2, 1.2)))
  }
}
## ---- Legend (single legend; scenarios as lines + Full factorial as filled dot) ----
legend("topright", bty = "n", cex = 0.4,
       legend = c("Full factorial var.",unname(scen_labs[names(scen_cols)])),
       col    = c("black",unname(scen_cols)),
       lty    = c(NA, rep(1, length(scen_cols))),
       lwd    = c(NA, rep(2, length(scen_cols))),
       pch    = c(16, rep(NA, length(scen_cols))),
       pt.cex = c(0.5,rep(NA, length(scen_cols))))

dev.off()


for(ss in sort(unique(heterosis_ff_tc$split))){
  
  # --- MPH ---
  png(file.path(foldertosave, paste0("MPHtc_split_", ss, ".png")),
      width = 3000, height = 1200, res = 300)
  
  op <- par(mfrow = c(1, length(sort(unique(heterosis_ff_tc$meanDD)))),
            mar = c(5, 5, 4, 1) + 0.1,
            oma = c(0, 0, 2, 0))
  
  for(dd in sort(unique(heterosis_ff_tc$meanDD))){
    
    tmp <- subset(heterosis_ff_tc, split == ss & meanDD == dd)
    
    # Ensure same x-axis order across panels
    tmp$ntesters <- factor(tmp$ntesters, levels = sort(unique(heterosis_ff_tc$ntesters)))
    
    MPH_FF <- unique(tmp$MPH)[1]
    
    boxplot(MPH_tc ~ ntesters, data = tmp,
            col = "grey85", border = "grey30",
            outline = FALSE,
            xlab = "Number of testers (k)",
            ylab = "MPH (testcross GV)",
            main = paste0("split=", ss, " | meanDD=", dd))
    
    abline(h = MPH_FF, col = "royalblue", lwd = 2, lty = 1)
    
    legend("topright", bty = "n", cex = 0.9,
           legend = c("Testcross sets", "Full factorial MPH"),
           fill = c("grey85", NA),
           border = c("grey30", NA),
           lty = c(NA, 2), lwd = c(NA, 3),
           col = c("grey30", "royalblue"))
  }
  
  mtext(paste0("Testcross Mid Parent Heterosis (split=", ss, ")"),
        outer = TRUE, cex = 1.1)
  
  par(op)
  dev.off()
  
}

## Facet MPH plots: split panels side-by-side, meanDD panels within each split
## Same dimensions + shared title style as your ratiovar_A_D_global plot

splits_to_plot <- sort(unique(heterosis_ff_tc$split))
dd_levels <- sort(unique(heterosis_ff_tc$meanDD))
k_levels <- sort(unique(heterosis_ff_tc$ntesters))

png(file.path(foldertosave, "MPHtc_global.png"),
    width = 9000, height = 5500, res = 600)

par(mfrow = c(length(splits_to_plot), length(dd_levels)),
    mar = c(5, 5.0, 2.2, 0.3) + 0.1,
    oma = c(2.5, 0, 2.5, 0))

par(mgp = c(3.5, 1.2, 0))

for(ss in splits_to_plot){
  for(dd in dd_levels){
    
    tmp <- subset(heterosis_ff_tc, split == ss & meanDD == dd)
    tmp$ntesters <- factor(tmp$ntesters, levels = k_levels)
    
    MPH_FF <- unique(tmp$MPH)[1]
    
    boxplot(MPH_tc ~ ntesters, data = tmp,
            col = "grey85", border = "grey30",
            outline = FALSE, ylim = c(min(tmp$MPH_tc, na.rm = TRUE) - 1.1, max(tmp$MPH_tc, na.rm = TRUE) + 1.1),
            xlab = "", cex.axis = 1.8, cex.lab = 1.8, cex.main = 1.8,
            ylab = if (dd == dd_levels[1]) "MPH (testcross gv)" else "",
            main = paste0("Divergence ", ss, " generations | Dom. level ", dd),
            bty = "n")   # remove top/right box lines
    
    abline(h = MPH_FF, col = "royalblue", lwd = 2, lty = 1)
    
    # Show legend once (top-left panel)
    if(ss == splits_to_plot[1] && dd == dd_levels[3]){
      legend("topright", bty = "n", cex = 1.2,
             legend = c("Testcross sets", "Full factorial MPH"),
             fill = c("grey85", NA),
             border = c("grey30", NA),
             lty = c(NA, 1), lwd = c(NA, 2),
             col = c("grey30", "royalblue"))
    }
  }
}

mtext("Number of testers", side = 1, outer = TRUE, line = 0, cex = 1.5)
mtext("Testcross Mid-Parent Heterosis (MPH)", outer = TRUE, line = 0.5, cex = 1.5)

dev.off()


cat("Total time for loop:", round(runtime, 3), "minutes\n",
    "Split levels:", length(nGenSplit), "\n",
    "Dominance levels:", length(meanDD), "\n",
    "Tester levels:", length(ntesters), "\n",
    "Iterations: k=1 ->", 5 * (nparents/2), "; k>1 ->", nsamples, "\n")
