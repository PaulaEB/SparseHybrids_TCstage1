## Loop over dominance and tester levels ##
# PaulaE
rm(list = ls())
# .rs.restartR()

library(AlphaSimR)
library(FieldSimR)
library(dplyr)
library(tibble)
library(asreml)

# Testcross parameters
nparents <- 80 #Number of parents in both pools
nchr <- 10 #Number of chromosomes in maize
nqtl <- 500 #Lit causal QTL for a trait. So total is nchr x nqtl = 5000 QTL
nsnps <- 500 # Markers not necessarily in LD with QTL of with effects on the trait
nGenSplit = 50 #Heterotic pool split
# ncrosses <- 60

# Trait parameters
meanG <- 70
varG <- 20
h2 <- 0.3 # plot-level heritability, ratio of genetic to phenotypic variance, varG/(varG + varE)
meanDD<- c(1.2, 0.9, 0.5, 0.2, 0) #high 0.9, medium 0.5 and low 0.2 dominance according to Gonzalez-Dieguez 2025
# varDD <- 0.2

# Experimental design parameters: full factorial
nenvs <- 1
nblocks <- nreps <- 2
ngenos_ff <- (nparents/2)*(nparents/2)

# Experimental design parameters: testcross
nenvs <- 1
nblocks <- nreps <- 2

heterosis_dd <- list()
vd_dd <- list()
gcaff_var_dd <- list()

# Loop over dominance levels #

t0 <- Sys.time()

for(d in meanDD){
  
  cat("Simulating founders with meanDD of ", d, "\n")
  
  dir.create(file.path("results", paste0(d,"DD")),
             showWarnings = FALSE, recursive = T)
  
  founderpop <- runMacs(nInd=nparents, 
                        nChr=nchr, 
                        segSites= (nsnps+nqtl), #number of segregating sites per chromosome 1000 * 10= 10,000 sites
                        inbred=TRUE,
                        species="MAIZE",
                        split = nGenSplit)
  
  SP <- SimParam$new(founderpop)
  SP$setTrackPed(TRUE)
  # SP$restrSegSites(minSnpPerChr = nsnps, overlap = F)
  SP$addSnpChip(nsnps)
  
  if (d == 0) {varDD <- 0} else {varDD <- 0.2}
  SP$addTraitAD(nQtlPerChr=nqtl, mean = meanG, var = varG, meanDD = d, varDD = varDD)
  
  poolA <- newPop(founderpop[1:(nparents/2)])
  poolB <- newPop(founderpop[((nparents/2)+1):(nparents)])
  
  ##Population structure and allele frequencies between pools##
  
  M_segsitesA <- pullSegSiteGeno(poolA)
  M_snpsA <- pullSnpGeno(poolA)
  M_qtlA <- pullQtlGeno(poolA)
  alleleFreqs_sitesA <- colMeans(M_segsitesA)/2
  alleleFreqs_snpsA <- colMeans(M_snpsA)/2
  alleleFreqs_qtlA <- colMeans(M_qtlA)/2
  
  M_segsitesB <- pullSegSiteGeno(poolB)
  M_snpsB <- pullSnpGeno(poolB)
  M_qtlB <- pullQtlGeno(poolB)
  alleleFreqs_sitesB <- colMeans(M_segsitesB)/2
  alleleFreqs_snpsB <- colMeans(M_snpsB)/2
  alleleFreqs_qtlB <- colMeans(M_qtlB)/2
  
  cat(" Overlap between QTL and SNPs in pool A: ", sum(colnames(M_qtlA) %in% colnames(M_snpsA)), "\n")
  overlap_A <- sum(colnames(M_qtlA) %in% colnames(M_snpsA))
  
  cat("Correlation between allele frequencies of all the sites", cor(alleleFreqs_sitesA, alleleFreqs_sitesB), "\n") 
  cor_allele_freqssites <- cor(alleleFreqs_sitesA, alleleFreqs_sitesB) 
  cor_allele_freqssnps <- cor(alleleFreqs_snpsA, alleleFreqs_snpsB) 
  cor_allele_freqsqtl <- cor(alleleFreqs_qtlA, alleleFreqs_qtlB)  
  
   
  png(file.path("results", paste0(d,"DD"),
                paste0("1dd_", d,"_corallfreqssites.png")), width = 2400, height = 1800, res = 600)
  plot(alleleFreqs_sitesA, alleleFreqs_sitesB,
       xlab = "Allele frequencies Pool A",
       ylab = "Allele frequencies Pool B",
       main = paste0("All sites allele frequencies", "\n", "between pools with split ", nGenSplit, "\n", "DD", d),
       pch = 1, col = "grey60")
  abline(0,1, lwd = 2)
  legend("topleft",
         legend = sprintf("correlation = %.3f", cor(alleleFreqs_sitesA, alleleFreqs_sitesB)),
         bty = "n")
  dev.off()
  
  png(file.path("results", paste0(d,"DD"),
                paste0("1dd_", d,"_corallfreqssnps.png")), width = 2400, height = 1800, res = 600)
  plot(alleleFreqs_snpsA, alleleFreqs_snpsB,
       xlab = "Allele frequencies Pool A",
       ylab = "Allele frequencies Pool B",
       main = paste0("SNPs allele frequencies", "\n", "between pools with split ", nGenSplit, "\n", "DD", d),
       pch = 1, col = "grey60")
  abline(0,1, lwd = 2)
  legend("topleft",
         legend = sprintf("correlation = %.3f", cor(alleleFreqs_snpsA, alleleFreqs_snpsB)),
         bty = "n")
  dev.off()
  
  png(file.path("results", paste0(d,"DD"),
                paste0("1dd_", d,"_corallfreqsqtl.png")), width = 2400, height = 1800, res = 600)
  plot(alleleFreqs_qtlA, alleleFreqs_qtlB,
       xlab = "Allele frequencies Pool A",
       ylab = "Allele frequencies Pool B",
       main = paste0("QTL allele frequencies", "\n", "between pools with split ", nGenSplit, "\n", "DD", d),
       pch = 1, col = "grey60")
  abline(0,1, lwd = 2)
  legend("topleft",
         legend = sprintf("correlation = %.3f", cor(alleleFreqs_qtlA, alleleFreqs_qtlB)),
         bty = "n")
  dev.off()

  G <- rbind(M_snpsA, M_snpsB)
  pools <- factor(c(rep("A", nrow(M_snpsA)), rep("B", nrow(M_snpsB))))
  pc <- prcomp(G, center = TRUE, scale. = FALSE)
  summary(pc)
  M_vaf1 <- 100*pc$sdev[1]^2/sum(pc$sdev^2)
  M_vaf2 <- 100*pc$sdev[2]^2/sum(pc$sdev^2)
  
  png(file.path("results", paste0(d,"DD"),
                paste0("2dd_", d,"pca_pools.png")), width = 2600, height = 2200, res = 600)
  plot(pc$x[,1], pc$x[,2],
       col = c("A"="purple","B"="black")[pools], pch = 19,
       xlab = sprintf("PC1 (%.1f%%)", M_vaf1),
       ylab = sprintf("PC2 (%.1f%%)", M_vaf2), main = paste0("SNPs Between pools DD",d))
  legend("topright", c("Pool A","Pool B"), col = c("purple","black"), pch = 19, bty = "n")
  dev.off()
  
  ## Full factorial cross ##
  
  hybridpop <- hybridCross(females = poolA,
                           males = poolB, 
                           crossPlan = "testcross")
  
  png(file.path("results", paste0(d,"DD"),
                paste0("3dd_", d,"hist_gv.png")), width = 2800, height = 2500, res = 600)
  hist(hybridpop@gv[,1], breaks = 10, freq = FALSE,
       xlim = range(c(hybridpop@gv[,1], poolA@gv[,1], poolB@gv[,1])),
       col = adjustcolor("#FFC20A", 0.30), border = NA,
       main = sprintf("GV Trait1 DD",d, "\n","Hmean%.2f var=%.2f\nPoolA mean=%.2f var=%.2f\nPoolB mean=%.2f var=%.2f",
                      mean(hybridpop@gv[,1]), popVar(hybridpop@gv[,1]),
                      mean(poolA@gv[,1]),     popVar(poolA@gv[,1]),
                      mean(poolB@gv[,1]),     popVar(poolB@gv[,1])), cex.main = 0.5)
    # Add pool A
  hist(poolA@gv[,1], breaks = 10, add = TRUE, freq = FALSE,
       col = adjustcolor("purple2", 0.25), border = NA)
    # Add pool B
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
  
  M_snpshybrids <- pullSnpGeno(hybridpop)
  G_all<-rbind(M_snpsA, M_snpsB, M_snpshybrids)
  grp <- factor(c(rep("H", nrow(M_snpshybrids)),
                  rep("A", nrow(M_snpsA)),
                  rep("B", nrow(M_snpsB))),
                levels = c("H","A","B"))
  
  pc_all <- prcomp(G_all, center = TRUE, scale. = FALSE)
  varExp <- pc_all$sdev^2 / sum(pc_all$sdev^2)
  
  colors <- c(H = "#FFC20A", A = "purple", B = "black")
  pchs <- c(H = 1, A = 19, B = 19)
  
  png(file.path("results", paste0(d,"DD"),
                paste0("2dd_", d,"pca_all.png")), width = 2600, height = 2200, res = 600)
  plot(pc_all$x[grp=="H",1], pc_all$x[grp=="H",2],
       col = colors["H"], pch = pchs["H"], cex = 0.9,
       xlab = sprintf("PC1 (%.1f%%)", 100*varExp[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100*varExp[2]),
       main = paste0("SNPs, all pops DD", d))
  
  points(pc_all$x[grp=="A",1], pc_all$x[grp=="A",2],
         col = colors["A"], pch = pchs["A"], cex = 1.2)
  points(pc_all$x[grp=="B",1], pc_all$x[grp=="B",2],
         col = colors["B"], pch = pchs["B"], cex = 1.2)
  
  legend("topright",
         legend = c("Hybrids","Pool A","Pool B"),
         col = colors[c("H","A","B")],
         pch = pchs[c("H","A","B")],
         pt.cex = c(0.9, 1.2, 1.2),
         bty = "n")
  dev.off()
  
  ## Heterosis and dominance levels ##

  addEff = SP$traits[[1]]@addEff
  domEff = SP$traits[[1]]@domEff
  
  mean_pool_allfreqsites <- mean(abs(alleleFreqs_sitesA - alleleFreqs_sitesB))
  mean_pool_allfreqqtl <- mean(abs(alleleFreqs_qtlA - alleleFreqs_qtlB))
  
  HF1 <- sum(domEff * (alleleFreqs_qtlA - alleleFreqs_qtlB)^2)
  MPH <- mean(hybridpop@gv[,1]) - (mean(poolA@gv[,1]) + mean(poolB@gv[,1]))/2
  BPH <- mean(hybridpop@gv[,1]) - max(mean(poolA@gv[,1]), mean(poolB@gv[,1]))
  
  heterosis_dd[[as.character(d)]] <- data.frame(meanDD = d, 
                                                cor_allele_freqssites = cor_allele_freqssites,
                                                cor_allele_freqsqtl = cor_allele_freqsqtl,
                                                mean_pool_allfreqqtl = mean_pool_allfreqqtl,
                                                mean_pool_allfreqsites = mean_pool_allfreqsites,
                                                HF1 = HF1, MPH = MPH, BPH = BPH)
  
  vd_dd[[as.character(d)]] <- data.frame(meanDD = d, 
                                         varA = varA(hybridpop), 
                                         varD = varD(hybridpop), 
                                         varG = varG(hybridpop), 
                                         varAovervarG = varA(hybridpop) / varG(hybridpop),
                                         varDovervarG = varD(hybridpop) / varG(hybridpop))

  ##Check with DT
  
  # a <- SP$traits[[1]]@addEff
  # dEff <- SP$traits[[1]]@domEff
  # 
  # Ai <- function(M_qtl) as.vector((M_qtl - 1) %*% a)
  # Di <- function(M_qtl) as.vector((M_qtl == 1) %*% dEff)
  # 
  # A_poolA <- Ai(pullQtlGeno(poolA))
  # D_poolA <- Di(pullQtlGeno(poolA))
  # 
  # A_poolB <- Ai(pullQtlGeno(poolB))
  # D_poolB <- Di(pullQtlGeno(poolB))
  # 
  # A_hyb <- Ai(pullQtlGeno(hybridpop))
  # D_hyb <- Di(pullQtlGeno(hybridpop))
  # 
  # MPH_A <- mean(A_hyb) - (mean(A_poolA) + mean(A_poolB))/2
  # MPH_D <- mean(D_hyb) - (mean(D_poolA) + mean(D_poolB))/2
  # 
  # MPH_G <- mean(hybridpop@gv[,1]) - (mean(poolA@gv[,1]) + mean(poolB@gv[,1]))/2
  
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
  pop_gv_mean_trait1 <- mean(gv_ff_matrix) #mean of the population in general, same as meanG = 70. Here as a matrix it takes correctly the mean of all the gvs centered 
  gv_ff_matrix <- gv_ff_matrix - pop_gv_mean_trait1 #genetic values without the intercept, so we can calculate the GCA effects directly
  
  # Calculate true GCA effects (all parent combos) for pool A trait 1
  gca_full_A_trait1 <- data.frame(mother_id = rownames(gv_ff_matrix),
                                  GCA = rowMeans(gv_ff_matrix))
 
  png(file.path("results", paste0(d,"DD"),
                paste0("4dd_", d,"histGCA.png")), width = 2600, height = 2000, res = 600)
  hist(gca_full_A_trait1$GCA, breaks = 20, 
       main = sprintf("GCA effects Pool A DD %.2f\nmean=%.2f var=%.2f",
                      d, mean(gca_full_A_trait1$GCA), popVar(gca_full_A_trait1$GCA)))
  abline(v = mean(gca_full_A_trait1$GCA), col = "purple2", lty = 1, lwd = 3) 
  dev.off() 
  
  gca_full_B_trait1 <- data.frame(father_id = colnames(gv_ff_matrix),
                                  GCA = colMeans(gv_ff_matrix))
  
  sca_table_trait1 <- sweep(gv_ff_matrix, 1, gca_full_A_trait1$GCA, FUN = "-") 
  sca_table_trait1 <- sweep(sca_table_trait1, 2, gca_full_B_trait1$GCA, FUN = "-") #calculate the SCA effects as the deviation from the GCA effects of both parents
    
  sca_full_trait1 <- data.frame(mother_id = rep(rownames(gv_ff_matrix), each = ncol(gv_ff_matrix)),
                              father_id = as.numeric(rep(colnames(gv_ff_matrix), times = nrow(gv_ff_matrix))),
                              SCA = as.vector(sca_table_trait1))
  
  png(file.path("results", paste0(d,"DD"),
                paste0("4dd_", d,"histSCA.png")), width = 2600, height = 2000, res = 600)
  hist(sca_full_trait1$SCA, breaks = 20,
       main = sprintf("SCA effects A DD %.2f\nmean=%.2f var=%.2f",
                      d, mean(sca_full_trait1$SCA), popVar(sca_full_trait1$SCA)))
  abline(v = mean(sca_full_trait1$SCA), col = "#FFC20A", lty = 1, lwd = 3)
  dev.off()
  
  gcaff_var_dd[[as.character(d)]] <- data.frame(meanDD = d,
                                                var_GCA_A = popVar(gca_full_A_trait1$GCA),
                                                var_GCA_B = popVar(gca_full_B_trait1$GCA),
                                                var_SCA   = popVar(sca_full_trait1$SCA))
  
  
 ## GRM ##
 
  Z <- sweep(M_snpsA, 2, 2*alleleFreqs_snpsA, "-") # the markers matrix from the SNPs (M_snpsA) in the mothers, apply by columns 2, 2pj substract so the mean is 0
  K <- 2 * sum(alleleFreqs_snpsA * (1 - alleleFreqs_snpsA), na.rm = TRUE)
  GRM_A <- (Z %*% t(Z)) #/ K
  dim(GRM_A)
  svd(GRM_A)$d
  GRM_A <- GRM_A + diag(1e-5, nrow(GRM_A)) #add a small value to the diagonal to make it positive definite
  GRM_A <- GRM_A/mean(diag(GRM_A)) #standardize the GRM so the mean of the diagonal is 1
  
  ## Loop over testers ##
  
  ntesters <- c(1,3,5,10,15)
  nsamples <- 120
  top_n <- 10
  top_true <- as.numeric(gca_full_A_trait1$mother_id[order(gca_full_A_trait1$GCA, decreasing = TRUE)][1:top_n])
  cor_gca_full_list <- cor_gca_tc_list <- meangca_list <- varE_list <- vargcas_scas_full_tc_list <- vargcaA_models_list <- ranking_intercept <- idx_list <- c()
  # cor_sca_tc_list
  ctr <- 0
  
  for(k in ntesters){
    niter <- if (k == 1) 2 * (length(poolB@id)) else nsamples
    
    for(i in 1:niter){
      ctr <- ctr + 1
      cat("ntesters:",k, "it:",i,"\n")
      
   
      # Subset of random testers from poolB
      if (k == 1) {
        testers_idx <- poolB@id[((i - 1) %% (length(poolB@id))) + 1]                  
      } else {
        testers_idx <- sample(poolB@id, k, replace = FALSE)   
      }
      
      tester_set <- paste(sort(testers_idx), collapse = ",")
      
      # Correlation between the true GCA females from the full factorial cross and the GCA from the 3 testers
      gca_testcross_A_trait1 <- data.frame(mother_id = rownames(gv_ff_matrix),
                                           GCA_tc = rowMeans(cbind(gv_ff_matrix[, testers_idx])))
      
      gca_testcross_B_trait1 <- data.frame(father_id = as.character(testers_idx),
                                           GCA_tc = NA)
      
      sca_testcross_trait1 <- data.frame(mother_id = rep(rownames(sca_table_trait1), times = length(testers_idx)),
                                         father_id = rep(testers_idx, each = nrow(sca_table_trait1)),
                                         SCA_tc    = NA)
      if (k > 1) {
        
      gca_testcross_B_trait1 <- data.frame(father_id = as.character(testers_idx),
                                           GCA_tc    = colMeans(gv_ff_matrix[, as.character(testers_idx), drop = FALSE]))
        
      sca_testcross_trait1 <- data.frame(mother_id = rep(rownames(sca_table_trait1), times = length(testers_idx)),
                                           father_id = rep(testers_idx, each = nrow(sca_table_trait1)),
                                           SCA_tc    = as.vector(sca_table_trait1[, as.character(testers_idx)])) %>%
        arrange(father_id)
        
      
      }
      
      # Environmental variance 
      gv_ff_df_k_testers <- droplevels(gv_ff_df[gv_ff_df$father_poolB %in% testers_idx, ])
      
      # (varE <- diag(varG(hybridpop))/h2 - diag(varG(hybridpop))) 
      
      ngenos_tc_k <- 40 * k
      nplots_tc_k <- ngenos_tc_k * nreps
      ncols_tc_k <- 10
      nrows_tc_k <- nplots_tc_k / ncols_tc_k
      
      varE <- popVar(gv_ff_df_k_testers$gv.TraitA)/h2 - popVar(gv_ff_df_k_testers$gv.TraitA)
      
      
      # var_pops <- data.frame(PoolA = c(mean(poolA@gv[,1]), varA(poolA), varD(poolA), varA(poolA)+varD(poolA), varG(poolA), varA(poolA)/varG(poolA), varD(poolA)/varG(poolA),varE),
      #                        PoolB = c(mean(poolB@gv[,1]), varA(poolB), varD(poolB), varA(poolB)+varD(poolB), varG(poolB), varA(poolB)/varG(poolB), varD(poolB)/varG(poolB),varE),
      #                        Hybrids = c(mean(hybridpop@gv[,1]), varA(hybridpop), varD(hybridpop), varA(hybridpop)+varD(hybridpop), varG(hybridpop), varA(hybridpop)/varG(hybridpop), varD(hybridpop)/varG(hybridpop),varE),
      #                        row.names = c("meangv", "varA", "varD","varA+D", "varG", "varA/varG", "varD/varG", "varE"))
      
      error_df <- FieldSimR::field_trial_error(varR = varE,
                                               ntraits = 1,
                                               nenvs = nenvs,
                                               nblocks = nblocks,
                                               ncols = ncols_tc_k,  
                                               nrows = nrows_tc_k,
                                               spatial.model = "AR1")
      
      gv_tc_df <- gv_ff_df_k_testers[, c(1, 4, 5, 6)]
      pheno_tc_df <- make_phenotypes(gv.df = gv_tc_df, error_df, randomise = T)
      pedigree <- gv_ff_df_k_testers %>% select(id, mother_poolA, father_poolB) %>% distinct(id, .keep_all = TRUE)
      pheno_tc_df <- pheno_tc_df %>% left_join(pedigree, by = "id")
      
      #Obtaining GCA from the phenotypic data without a model
      pheno_gca_A_nomodels <- pheno_tc_df %>% group_by(mother_poolA) %>% summarise(GCA_A_pheno = mean(y.Trait1)-mean(pheno_tc_df$y.Trait1)) %>% ungroup()
      
      #Obtaining SCA from the phenotypic data when more than one tester
      if (k > 1) {
      pheno_gca_B_nomodels <- pheno_tc_df %>% group_by(father_poolB) %>% summarise(GCA_B_pheno = mean(y.Trait1)-mean(pheno_tc_df$y.Trait1)) %>% ungroup()
      
      mu_pheno <- mean(pheno_tc_df$y.Trait1)
      
      pheno_sca_nomodels <- pheno_tc_df %>% group_by(mother_poolA, father_poolB) %>% 
        summarise(y_ij = mean(y.Trait1) - mu_pheno, .groups = "drop") %>% #mean per the 2 reps and without the mean to have the same units as gcas
        left_join(pheno_gca_A_nomodels, by = "mother_poolA") %>% left_join(pheno_gca_B_nomodels, by = "father_poolB") %>% #join the gcas per mom and dad
        mutate(SCA_pheno = y_ij - GCA_A_pheno - GCA_B_pheno) %>% #calculate the sca as the deviation from the gcas
        select(mother_poolA, father_poolB, SCA_pheno) %>% #remove the gcas
        arrange(father_poolB) #order the scas by mom and dad to have the same order as the sca_testcross_trait1
      }
      
      # Obtain the GCA using the three ASReml models
      
      model_testcross1 <- asreml(fixed = y.Trait1 ~ father_poolB - 1, 
                                 random = ~ mother_poolA + block, 
                                 residual = ~units,
                                 data = pheno_tc_df)
      
      model_testcross2 <- asreml(fixed = y.Trait1 ~ father_poolB - 1, 
                                 random = ~ mother_poolA +
                                   block, 
                                 residual = ~ar1(col):ar1(row),
                                 data = pheno_tc_df)
      
      model_testcross3 <- asreml(fixed = y.Trait1 ~ father_poolB - 1, 
                                 random = ~ vm(mother_poolA, GRM_A) +                                                      
                                   block, 
                                 residual = ~ar1(col):ar1(row),
                                 data = pheno_tc_df)
      
      m1 <- model_testcross1$coefficients$random %>% as.data.frame() %>% 
        rownames_to_column(var = "mother_id") %>% filter(grepl("^mother_poolA", mother_id)) %>% 
        mutate(mother_id = sub(".*_", "", mother_id), effect) 
      
      gca_vc_m1 <- as.numeric(summary(model_testcross1)$varcomp["mother_poolA", "component"])
      
      m2 <- model_testcross2$coefficients$random %>% as.data.frame() %>% 
        rownames_to_column(var = "mother_id") %>% filter(grepl("^mother_poolA", mother_id)) %>% 
        mutate(mother_id = sub(".*_", "", mother_id), effect)
      
      gca_vc_m2 <- as.numeric(summary(model_testcross2)$varcomp["mother_poolA", "component"])
      
      m3 <- model_testcross3$coefficients$random %>% as.data.frame() %>% 
        rownames_to_column(var = "mother_id") %>% filter(grepl("mother_poolA", mother_id)) %>% 
        mutate(mother_id = sub(".*_", "", mother_id), effect)
      gca_vc_m3 <- as.numeric(summary(model_testcross3)$varcomp["vm(mother_poolA, GRM_A)", "component"])
      
      top_tc_gv <- as.numeric(gca_testcross_A_trait1$mother_id[order(gca_testcross_A_trait1$GCA_tc, decreasing = TRUE)][1:top_n])
      top_tc_pheno_nomodel <- as.numeric(pheno_gca_A_nomodels$mother_poolA[order(pheno_gca_A_nomodels$GCA_A_pheno, decreasing = TRUE)][1:top_n])
      top_tc_gca_model1 <- as.numeric(m1$mother_id[order(m1$effect, decreasing = TRUE)][1:top_n])
      top_tc_gca_model2 <- as.numeric(m2$mother_id[order(m2$effect, decreasing = TRUE)][1:top_n])
      top_tc_gca_model3 <- as.numeric(m3$mother_id[order(m3$effect, decreasing = TRUE)][1:top_n])
      
      idx_list[[ctr]] <- data.frame(ntesters = k,
                                    iter = i,
                                    tester_set = tester_set)
      
      cor_gca_full_list[[ctr]] <- data.frame(ntesters = k,
                                         iter = i,
                                         tester_set = tester_set,
                                         gv_tc  = cor(gca_full_A_trait1$GCA, gca_testcross_A_trait1$GCA_tc),
                                         pheno  = cor(gca_full_A_trait1$GCA, pheno_gca_A_nomodels$GCA_A_pheno),
                                         model1 = cor(gca_full_A_trait1$GCA, m1$effect),
                                         model2 = cor(gca_full_A_trait1$GCA, m2$effect),
                                         model3 = cor(gca_full_A_trait1$GCA, m3$effect))
      
      cor_gca_tc_list[[ctr]] <- data.frame(ntesters = k,
                                       iter = i,
                                       tester_set = tester_set,
                                       gv_tc  = cor(gca_full_A_trait1$GCA, gca_testcross_A_trait1$GCA_tc),
                                       pheno  = cor(gca_testcross_A_trait1$GCA_tc, pheno_gca_A_nomodels$GCA_A_pheno),
                                       model1 = cor(gca_testcross_A_trait1$GCA_tc, m1$effect),
                                       model2 = cor(gca_testcross_A_trait1$GCA_tc, m2$effect),
                                       model3 = cor(gca_testcross_A_trait1$GCA_tc, m3$effect))
      ##Check with DT
      # cor_sca_tc_list[[ctr]] <- data.frame(ntesters = k,
      #                                  iter = i,
      #                                  tester_set = tester_set,
      #                                  pheno  = cor(sca_testcross_trait1$SCA_tc, pheno_sca_nomodels$SCA_pheno))
      
      meangca_list[[ctr]] <- data.frame(ntesters = k,
                                     iter = i,
                                     tester_set = tester_set,
                                     gv_full = mean(gca_full_A_trait1$GCA),
                                     gv_tc   = mean(gca_testcross_A_trait1$GCA_tc),
                                     pheno   = mean(pheno_gca_A_nomodels$GCA_A_pheno),
                                     model1  = mean(m1$effect),
                                     model2  = mean(m2$effect),
                                     model3  = mean(m3$effect))
      
      varE_list[[ctr]] <- data.frame(ntesters = k,
                                     iter = i,
                                     tester_set = tester_set,
                                     varE = varE)
      
      vargcas_scas_full_tc_list[[ctr]] <- data.frame(ntesters = k,
                                                     iter = i,
                                                     tester_set = tester_set,
                                                     gcaA_full = popVar(gca_full_A_trait1$GCA),
                                                     gcaB_full = popVar(gca_full_B_trait1$GCA),
                                                     sca_full = popVar(sca_full_trait1$SCA),
                                                     gcaA_tc   = popVar(gca_testcross_A_trait1$GCA_tc),
                                                     gcaB_tc = if (k > 1) popVar(gca_testcross_B_trait1$GCA_tc) else NA,
                                                     sca_tc  = if (k > 1) popVar(sca_testcross_trait1$SCA_tc) else NA)
      
      vargcaA_models_list[[ctr]] <- data.frame(ntesters = k,
                                               iter = i,
                                               tester_set = tester_set,
                                               gv_full = popVar(gca_full_A_trait1$GCA),
                                               gv_tc   = popVar(gca_testcross_A_trait1$GCA_tc),
                                               pheno   = popVar(pheno_gca_A_nomodels$GCA_A_pheno),
                                               model1  = gca_vc_m1,
                                               model2  = gca_vc_m2,
                                               model3  = gca_vc_m3)
      
      ranking_intercept[[ctr]] <- data.frame(ntesters = k,
                                             iter = i,
                                             tester_set = tester_set,
                                             method = c("gv_tc", "pheno", "model1", "model2", "model3"),
                                             overlap = c(gv_tc  = length(intersect(top_true, top_tc_gv)),
                                                         pheno  = length(intersect(top_true, top_tc_pheno_nomodel)),
                                                         model1 = length(intersect(top_true, top_tc_gca_model1)),
                                                         model2 = length(intersect(top_true, top_tc_gca_model2)),
                                                         model3 = length(intersect(top_true, top_tc_gca_model3))))
      
    
       } }
  # Summaries/plots
  
  cor_full_df <- dplyr::bind_rows(cor_gca_full_list)
  write.csv(cor_full_df, file.path("results", paste0(d,"DD"), paste0("cor_full_df_DD",d,".csv")), row.names = FALSE)
  cor_tc_df   <- dplyr::bind_rows(cor_gca_tc_list)
  write.csv(cor_tc_df, file.path("results", paste0(d,"DD"), paste0("cor_tc_df_DD",d,".csv")), row.names = FALSE)
  mean_df     <- dplyr::bind_rows(meangca_list)
  write.csv(mean_df, file.path("results", paste0(d,"DD"), paste0("mean_df_DD",d, ".csv")), row.names = FALSE)
  vargcas_scas_full_tc_list_df <- dplyr::bind_rows(vargcas_scas_full_tc_list)
  write.csv(vargcas_scas_full_tc_list_df, file.path("results", paste0(d,"DD"), paste0("vargcas_scas_full_tc_list_df_DD",d,".csv")), row.names = FALSE)
  vargcaA_models_df      <- dplyr::bind_rows(vargcaA_models_list)
  write.csv(vargcaA_models_df, file.path("results", paste0(d,"DD"), paste0("vargca_df_DD", d,".csv")), row.names = FALSE)
  rank_df     <- dplyr::bind_rows(ranking_intercept)
  write.csv(rank_df, file.path("results", paste0(d,"DD"), paste0("rank_df_DD",d,".csv")), row.names = FALSE)
  idx_df      <- dplyr::bind_rows(idx_list)
  write.csv(idx_df, file.path("results", paste0(d,"DD"), paste0("idx_df_DD",d,".csv")), row.names = FALSE)
  # write.csv(var_pops, file.path("results", paste0(d,"DD"), paste0("var_pops_DD",d,".csv")), row.names = TRUE)
  
  cor_full_mean <- aggregate(cbind(gv_tc, pheno, model1, model2, model3) ~ ntesters,
                             data = cor_full_df, FUN = mean, na.rm = TRUE)
  
  png(file.path("results", paste0(d,"DD"),
                paste0("5dd_",d,"_cortoFF.png")), width = 2100, height = 1500, res = 300)
  plot(cor_full_mean$ntesters, cor_full_mean$gv_tc, col = "green", ylim = c(0,1), 
       type="l", lwd=3,xaxt="n",
       main = paste0("Cor FF with TC n testers", "\n","dominance",d))
  axis(1, at = cor_full_mean$ntesters)
  lines(cor_full_mean$ntesters, cor_full_mean$pheno, pch = 16, col = "royalblue")
  lines(cor_full_mean$ntesters, cor_full_mean$model1, lty=3, lwd=3,pch = 16, col = "pink")
  lines(cor_full_mean$ntesters, cor_full_mean$model2, pch = 16, col = "magenta")
  lines(cor_full_mean$ntesters, cor_full_mean$model3, lty=3, lwd=3,pch = 16, col = "aquamarine")
  legend("bottomright", legend = c("gv_tc","pheno","model1","model2","model3"),
         col = c("green","royalblue","pink","magenta","aquamarine"), lwd =1, bty = "n")
  dev.off()
  
  cor_tc_mean <- aggregate(cbind(gv_tc, pheno, model1, model2, model3) ~ ntesters,
                           data = cor_tc_df, FUN = mean, na.rm = TRUE)
  
  png(file.path("results", paste0(d,"DD"),
                paste0("5dd_",d,"_cortoTC.png")), width = 2100, height = 1500, res = 300)
  plot(cor_tc_mean$ntesters, cor_tc_mean$gv_tc, col = "green", ylim = c(0,1), type="l",lwd=3, xaxt="n",
       main = paste0("Cor TC with TC n testers", "\n","dominance",d))
  axis(1, at = cor_tc_df$ntesters)
  lines(cor_tc_mean$ntesters, cor_tc_mean$pheno, pch = 16, col = "royalblue")
  lines(cor_tc_mean$ntesters, cor_tc_mean$model1, lty=3, lwd=3,pch = 16, col = "pink")
  lines(cor_tc_mean$ntesters, cor_tc_mean$model2, pch = 16, col = "magenta")
  lines(cor_tc_mean$ntesters, cor_tc_mean$model3, lty=3, lwd=3,pch = 16, col = "aquamarine")
  legend("bottomright", legend = c("gv_tc","pheno","model1","model2","model3"),
         col = c("green","royalblue","pink","magenta","aquamarine"), lwd =1, bty = "n")
  dev.off()
  
  if (any(cor_tc_df$ntesters == 1)) {
    
    cor_1tc_mean <- aggregate(cbind(gv_tc, pheno, model1, model2, model3) ~ tester_set,
                             data = cor_tc_df[cor_tc_df$ntesters == 1, ], FUN = mean, na.rm = TRUE) |> arrange(desc(gv_tc)) 
    cor_1tc_mean$tester_set <- as.numeric(cor_1tc_mean$tester_set) 
    
    png(file.path("results", paste0(d,"DD"),
                  paste0("5dd_",d,"_cortoFF_1tester.png")), width = 2100, height = 1500, res = 300)
    plot(seq_along(cor_1tc_mean$tester_set), cor_1tc_mean$gv_tc, col = "green", ylim = c(0,1), type="l",lwd=3, xaxt="n",
         main = paste0("Cor TC with TC 1 testers", "\n","dominance",d))
    axis(1,at = seq_along(cor_1tc_mean$tester_set), labels = cor_1tc_mean$tester_set, las = 2)
    lines(seq_along(cor_1tc_mean$tester_set), cor_1tc_mean$pheno, pch = 16, col = "royalblue")
    lines(seq_along(cor_1tc_mean$tester_set), cor_1tc_mean$model1, lty=3, lwd=3,pch = 16, col = "pink")
    lines(seq_along(cor_1tc_mean$tester_set), cor_1tc_mean$model2, pch = 16, col = "magenta")
    lines(seq_along(cor_1tc_mean$tester_set), cor_1tc_mean$model3, lty=3, lwd=3,pch = 16, col = "aquamarine")
    legend("bottomright", legend = c("gv_tc","pheno","model1","model2","model3"),
           col = c("green","royalblue","pink","magenta","aquamarine"), lwd =1, bty = "n")
    dev.off()
  } 
  
  vargca_mean <- aggregate(cbind(gv_full,gv_tc, pheno, model1, model2, model3) ~ ntesters,
                           data = vargcaA_models_df, FUN = mean, na.rm = TRUE)
  
  png(file.path("results", paste0(d,"DD"),
                paste0("6dd_",d,"_varGCA_A.png")), width = 2100, height = 1500, res = 300)
  plot(vargca_mean$ntesters, vargca_mean$gv_full, pch=10,col = "green", 
       type="l", lwd=2, xaxt="n", main = paste0("Var of GCA pool A models n testers", "\n","dominance",d),
       ylim= range(vargca_mean[, c("gv_full","gv_tc","pheno","model1","model2","model3")]))
  axis(1, at = vargca_mean$ntesters)
  lines(vargca_mean$ntesters, vargca_mean$gv_tc, pch = 16, col = "forestgreen")
  lines(vargca_mean$ntesters, vargca_mean$pheno, pch = 16, col = "royalblue")
  lines(vargca_mean$ntesters, vargca_mean$model1, lty=3, lwd=3,pch = 16, col = "pink")
  lines(vargca_mean$ntesters, vargca_mean$model2, pch = 16, col = "magenta")
  lines(vargca_mean$ntesters, vargca_mean$model3, lty=3, lwd=3,pch = 16, col = "aquamarine")
  # points(vargca_mean$ntesters, vargca_mean$gv_tc,  col="blue",    pch=16)
  # points(vargca_mean$ntesters, vargca_mean$pheno,  col="black",   pch=16)
  # points(vargca_mean$ntesters, vargca_mean$model1, col="pink",    pch=16)
  # points(vargca_mean$ntesters, vargca_mean$model2, col="magenta", pch=16)
  # points(vargca_mean$ntesters, vargca_mean$model3, col="purple4", pch=16)
  legend("topright", legend = c("gv_full","gv_tc","pheno","model1","model2","model3"),
         col = c("green","forestgreen","royalblue","pink","magenta","aquamarine"), lwd =1, bty = "n")
  dev.off()
  
  vargcas_scas_full_tc_mean <- aggregate(cbind(gcaA_full, gcaB_full, sca_full, gcaA_tc, gcaB_tc, sca_tc) ~ ntesters,
                                         data = vargcas_scas_full_tc_list_df,FUN = mean, na.rm = TRUE)
  
  png(file.path("results", paste0(d,"DD"), paste0("6dd_",d,"var_GCA_SCA.png")),  width = 2100, height = 1500, res = 300)
  
  plot(vargcas_scas_full_tc_mean$ntesters, vargcas_scas_full_tc_mean$gcaA_full, pch=10,col = "purple", 
       type="l", lwd=2, xaxt="n", main = paste0("Var of GCA pool A models n testers", "\n","dominance",d),
       ylim= range(vargcas_scas_full_tc_mean[, c("gcaA_full", "gcaB_full", "sca_full",  "gcaA_tc",  "gcaB_tc", "sca_tc")]))
  axis(1, at = vargcas_scas_full_tc_mean$ntesters)
  lines(vargcas_scas_full_tc_mean$ntesters, vargcas_scas_full_tc_mean$gcaA_tc, lty=2, lwd=3,pch = 16, col = "purple")
  lines(vargcas_scas_full_tc_mean$ntesters, vargcas_scas_full_tc_mean$gcaB_full, pch = 16, lwd=2, col = "black")
  lines(vargcas_scas_full_tc_mean$ntesters, vargcas_scas_full_tc_mean$gcaB_tc, lty=2, lwd=3, pch = 16, col = "black")
  lines(vargcas_scas_full_tc_mean$ntesters, vargcas_scas_full_tc_mean$sca_full, pch = 16, lwd=2, col = "#FFC20A")
  lines(vargcas_scas_full_tc_mean$ntesters, vargcas_scas_full_tc_mean$sca_tc, lty=2, lwd=3,pch = 16, col = "#FFC20A")
  # points(vargca_mean$ntesters, vargca_mean$gv_tc,  col="blue",    pch=16)
  # points(vargca_mean$ntesters, vargca_mean$pheno,  col="black",   pch=16)
  # points(vargca_mean$ntesters, vargca_mean$model1, col="pink",    pch=16)
  # points(vargca_mean$ntesters, vargca_mean$model2, col="magenta", pch=16)
  # points(vargca_mean$ntesters, vargca_mean$model3, col="purple4", pch=16)
  legend("topright", legend = c("FF_GCA_A","FF_GCA_B","FF_SCA","TC_GCA_A","TC_GCA_B","TC_SCA"),
         col = c("purple","black","#FFC20A","purple", "black","#FFC20A"),
         lty = c(1,1,1,2,2,2), lwd = 2, bty = "n", cex = 0.5)
  dev.off()
  png(file.path("results", paste0(d,"DD"),
                paste0("7dd_",d,"_overlaptesters.png")), width = 2100, height = 1500, res = 300)
  boxplot(overlap ~ ntesters, data = subset(rank_df, method == "pheno"), main= paste0("Overlap FF-TC of top ",top_n,"\n", "DD",d))
  dev.off()
  
  png(file.path("results", paste0(d,"DD"),
                paste0("8dd_",d,"_overlapmodels.png")), width = 2100, height = 1500, res = 300)
  boxplot(overlap ~ factor(method, levels = c("gv_tc","pheno","model1","model2","model3")), data = rank_df, main= paste0("Overlap FF-TC of top ",top_n,"\n", "DD",d))
  dev.off()
  
}

t1 <- Sys.time()  
runtime <- difftime(t1, t0, units = "mins")
cat("Total time for loop:", round(runtime, 3), "minutes\n",
  "Dominance levels:", length(meanDD), "\n",
  "Tester levels:", length(ntesters), "\n",
  "Iterations: k=1 ->", 5 * (nparents/2), "; k>1 ->", nsamples, "\n")

het_df <- dplyr::bind_rows(heterosis_dd) |> dplyr::arrange(meanDD)
dom_df <- dplyr::bind_rows(vd_dd) |> dplyr::arrange(meanDD)

colnames(dom_df)[1:6] <- c("meanDD", "varA", "varD", "varG", "varAovervarG", "varDovervarG")
write.csv(het_df, file.path("results", "heterosis_metrics_DD.csv"), row.names = FALSE)
write.csv(dom_df, file.path("results", "dominance_variance_DD.csv"), row.names = FALSE)

png(file.path("results", "heterosis.png"), width = 2100, height = 1500, res = 300)
matplot(het_df$meanDD,
        as.matrix(het_df[, c("HF1","MPH","BPH")]),
        type = "b", pch = 16, lwd = 3, lty = 1,
        col = c("grey40","royalblue","forestgreen"),
        xlab = "meanDD", ylab = "Value",
        main = "Heterosis metrics vs dominance")

legend("topleft", legend = c("HF1","MPH","BPH"),
       col = c("grey40","royalblue","forestgreen"),
       lwd = 3, pch = 16, bty = "n")
dev.off()

png(file.path("results", "dominancevar.png"), width = 2100, height = 1500, res = 300)
matplot(dom_df$meanDD,
        as.matrix(dom_df[, c("varAovervarG","varDovervarG")]),
        type = "b", pch = 16, lwd = 3, lty = 1,
        col = c("purple4","orange3"),
        xlab = "meanDD", ylab = "Fraction",
        ylim = c(0,1),
        main = "Genetic variance fractions in hybrids vs dominance")

legend("right", legend = c("varA/varG","varD/varG"),
       col = c("purple4","orange3"),
       lwd = 3, pch = 16, bty = "n")
dev.off()

gcaff_var_dd <- dplyr::bind_rows(gcaff_var_dd) %>% arrange(meanDD)
write.csv(gcaff_var_dd, file.path("results", "gca_metrics_DD.csv"), row.names = FALSE)
M <- t(as.matrix(gcaff_var_dd[, c("var_GCA_A","var_GCA_B","var_SCA")]))
colnames(M) <- paste0("DD=", gcaff_var_dd$meanDD)
rownames(M) <- c("GCA_A","GCA_B","SCA")
varff_scagca_ratio <- gcaff_var_dd$var_SCA / (gcaff_var_dd$var_GCA_A + gcaff_var_dd$var_GCA_B)

png(file.path("results", "FF_GCAs_SCA.png"), width = 2600, height = 1600, res = 300)
par(mar = c(8, 4, 4, 2) + 0.1)
text(x = barplot(M, beside = TRUE,
                 col = c("purple2", "grey10", "#FFC20A"),
                 ylab = "Variance", ylim = c(0, max(M, na.rm = TRUE) + 3),
                 main = "Full factorial GCAs", las = 2),
     y = as.vector(M),
     labels = round(as.vector(M), 2), pos = 3, cex = 0.8)

legend("topleft", legend = rownames(M), fill = c("purple2","grey10","#FFC20A"), bty = "n")
mtext(
  side = 1, line = 6,
  text = paste0("ratio varSCA/varGCAs: ", paste(round(varff_scagca_ratio, 2), collapse = " | ")),
  cex = 0.9, font = 2
)
dev.off()


runtime <- difftime(t1, t0, units = "mins")
cat("Total time for loop:", round(runtime, 3), "minutes\n",
    "Dominance levels:", length(meanDD), "\n",
    "Tester levels:", length(ntesters), "\n",
    "Iterations: k=1 ->", 5 * (nparents/2), "; k>1 ->", nsamples, "\n")

