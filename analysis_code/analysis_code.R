library(data.table)
library(dplyr)
library(forestploter)
library(geni.plots)
library(ggplot2)
library(ggrepel)
library(ggVolcano)
library(gwasglue)
library(grid)
library(ieugwasr) 
library(ldscr)
library(MungeSumstats)
library(paletteer)
library(plinkbinr)
library(psych)
library(purrr)
library(RColorBrewer)
library(reshape2)
library(readr)
library(tidyr)
library(TwoSampleMR)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(biomaRt)
library(boot)
library(caret)
library(readxl)
library(pheatmap)
library(pROC)
library(stringr)
library(survival)
library(survminer)
library(tidyr)
library(timeROC)
library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)
library(pROC)
library(Seurat)
library(utils)

setwd("D:/Desktop/ICT_HCC/")
#The key step code involved in the study

##### pQTL_type change #####
gtf_file <- "A:/Biology data/Plasma protein/gencode.v44lift37.annotation.gtf.gz"
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

genes <- genes(txdb)
gene_positions <- data.frame(seqnames = seqnames(genes),
                             start = start(genes),
                             end = end(genes),
                             ensgid = genes$gene_id)

gene_positions$ensgid <- sapply(strsplit(as.character(gene_positions$ensgid),
                                         "\\."), "[", 1)
gene_positions$symbol <- mapIds(org.Hs.eg.db, keys = gene_positions$ensgid,
                                column = "SYMBOL", keytype = "ENSEMBL")

annotate_and_save_pqtl_type <- function(file_path, out_file_name1, out_file_name2) {
  exp <- fread(file_path)
  exp <- as.data.frame(exp)
  names(exp)[1] <- "chromosome";names(exp)[2] <- "position"#
  exp$chromosome <- gsub("chr", "", exp$chromosome)
  exp$chromosome <- paste0("chr", exp$chromosome)
  
  cleaned_name <- strsplit(tools::file_path_sans_ext(basename(file_path)), "_")[[1]][2]
  target_gene <- subset(gene_positions, symbol == cleaned_name)
  
  if (nrow(target_gene) == 0) {
    empty_file_path <- paste(out_file_name1, "/", cleaned_name, ".txt", sep = "")
    file.create(empty_file_path)
    cat(paste("Empty file created:", empty_file_path), "\n")
    return()
  }
  
  exp_position <- exp$position
  exp_chromosome <- exp$chromosome
  target_gene_seqnames <- target_gene$seqnames
  matching_rows_start <- target_gene$start
  matching_rows_end <- target_gene$end
  exp$pqtl_type <- "trans-pqtl"
  
  for (i in 1:nrow(exp)) {
    current_chromosome <- exp_chromosome[i]
    matching_indices <- which(target_gene_seqnames == current_chromosome)
    
    if (length(matching_indices) > 0) {
      position <- exp_position[i]
      within_range <- (position >= matching_rows_start[matching_indices] - 1000000) & 
        (position <= matching_rows_end[matching_indices] + 1000000)
      
      if (any(within_range)) {
        exp$pqtl_type[i] <- "cis-pqtl"
      }
    }
  }
  
  exp_cis <- subset(exp, pqtl_type == "cis-pqtl")
  exp_trans <- subset(exp, pqtl_type == "trans-pqtl")
  exp_cis <- subset(exp_cis, select = -pqtl_type)
  exp_trans <- subset(exp_trans, select = -pqtl_type)
  
  new_file_path1 <- paste(out_file_name1, "/", tools::file_path_sans_ext(basename(file_path)), "_cis.gz", sep = "")
  new_file_path2 <- paste(out_file_name2, "/", tools::file_path_sans_ext(basename(file_path)), "_trans.gz", sep = "")
  
  write.csv(exp_cis, gzfile(new_file_path1), row.names = FALSE)
  write.csv(exp_trans, gzfile(new_file_path2), row.names = FALSE)
  print(cleaned_name)
}

a_files <- list.files("D:/Desktop/pwas/", 
                      pattern = "\\.gz$", full.names = TRUE)
out_file_name1 <- "E:/ICT/cis-pQTL"
out_file_name2 <- "E:/ICT/trans_pQTL"

for (file_path in a_files) {
  annotate_and_save_pqtl_type(file_path, out_file_name1, out_file_name2)
}

##### GWAS local clump #####
a <- read.table("A:/Biology data/Plasma protein/round_all_chrom_merged_lite_INFO.snpstats",header = T)
dup_index <- duplicated(a[c("chromosome", "position")])
dups <- a[dup_index, c("chromosome", "position")]
a <- a[!duplicated(a[c("chromosome", "position")]), ]

clump_local_data <- function(file_path, threshold, out_file_name){
  tryCatch({
    file_name <- basename(file_path)
    file_name_without_ext <- tools::file_path_sans_ext(file_name)
    file_name_without_ext <- strsplit(file_name_without_ext, "_qc_1e-4")[[1]][1]
    new_file_path <- paste(out_file_name, file_name_without_ext, "_", threshold, "clump.csv", sep = "")
    
    temp_dat <- read.table(file_path,header = T)
    temp_dat <- merge(temp_dat, a, by = c("chromosome", "position"))
    temp_dat$N <- 3301
    temp_dat <- temp_dat[, c(9, 10, 11, 13, 6, 7, 8, 19)]
    
    names(temp_dat) <- c("SNP","other_allele", "effect_allele", "eaf","beta","se","pval", "samplesize")
    temp_dat <- subset(temp_dat, pval < as.numeric(threshold))
    temp_dat <- format_data(dat = temp_dat,type = "exposure")

    if (T) {
      temp_dat$id <- temp_dat$id.exposure
      temp_dat$rsid <- temp_dat$SNP
      temp_dat$pval <- temp_dat$pval.exposure
      
      temp_dat <- ld_clump(temp_dat,
                           plink_bin = get_plink_exe(),
                           bfile = "A:/Biology data/Plasma protein/g1000_eur/g1000_eur",
                           clump_kb = 10000, clump_r2 = 0.001)
      
      temp_dat$rsid <- NULL
      temp_dat$pval <- NULL
      temp_dat$exposure <- file_name_without_ext
      write.csv(temp_dat,new_file_path,row.names = F)
      print(file_name_without_ext)
    }
  }, error = function(e) {
    message(paste("Error processing file:", file_path, "- Skipping to next file. Error message:", e$message))
  })
}

a_files <- list.files("A:\\Biology data\\Plasma protein\\INTERVAL\\INTERVAL_EUR_QC\\", 
                      pattern = "\\.txt$", full.names = TRUE)
for(file_path in a_files) {
  clump_local_data(file_path, "5e-8",
                   "A:\\Biology data\\Plasma protein\\INTERVAL\\INTERVAL_EUR_QC_clump_5e-8\\")
}

##### MR anlysis       #####
circulating_MR <- function(folder_path, save_path, exp_dat) {
  file_names <- list.files(folder_path, pattern = "\\.gz$", full.names = TRUE)
  for (file_name in file_names) {
    c <- fread(file_name)
    c <- as.data.frame(c)
    exp_dat <- subset(exp_dat,mr_keep.exposure == TRUE)
    out <- format_data(c,
                       type = "outcome",
                       snps = exp_dat$SNP,
                       header = TRUE,
                       phenotype_col = "Phenotype",
                       snp_col = "rsids",
                       beta_col = "beta",
                       se_col = "sebeta",
                       eaf_col = "af_alt",
                       effect_allele_col = "alt",
                       other_allele_col = "ref",
                       pval_col = "pval")
    
    mydata <- harmonise_data(
      exposure_dat=exp_dat,
      outcome_dat=out,
      action=2)
    mydata <- subset(mydata, pval.outcome >= 5e-8)
    mydata <- subset(mydata, mr_keep == TRUE)
    
    if (nrow(mydata) == 0){next
    }
    
    res_mi <- generate_odds_ratios(mr(mydata))
    res_mi$outcome <- file_name
    res_het <- mr_heterogeneity(mydata)
    res_ple <- mr_pleiotropy_test(mydata)
    
    mydata$R2_2 <- 2*mydata$eaf.exposure*(1-mydata$eaf.exposure)*(mydata$beta.exposure^2)
    exposure_list <- split(mydata, mydata$exposure)
    
    data_F <- data.frame(exposure = character(),F_statistic_b_R2_2 = numeric(),stringsAsFactors = FALSE ) 
    
    for (i in seq_along(exposure_list)) {  
      df <- exposure_list[[i]]  
      F_statistic_b_R2_2 <- min(sum(df$R2_2) * (df$samplesize.exposure - nrow(df) - 1) / (nrow(df) * (1 - sum(df$R2_2))))  
      
      data_F <- rbind(data_F, data.frame(  
        exposure = names(exposure_list)[i],  
        F_statistic_b_R2_2 = F_statistic_b_R2_2,  
        stringsAsFactors = FALSE  
      ))  
    }  
    
    if (all(is.na(res_het))) {
      merged_df <- cbind(res_mi, matrix(NA, nrow = nrow(res_mi), ncol = 3))
    } else {
      merged_df <- left_join(res_mi, res_het, by = c("id.exposure", "id.outcome", "method"))
      merged_df <- merged_df[, -c((ncol(merged_df) - 3):(ncol(merged_df) - 4))]
    }
    if (all(is.na(res_het))) {
      merged_df <- cbind(merged_df, matrix(NA, nrow = nrow(merged_df), ncol = 3))
    } else {
      merged_df <- left_join(merged_df, res_ple, by = c("id.exposure", "id.outcome"))
      merged_df <- merged_df[, -c((ncol(merged_df) - 3):(ncol(merged_df) - 4))]
    }
    colname1 <- c(colnames(res_mi),"Q","Q_df","Q_pval", "egger_intercept","egger_intercept_se"," egger_intercept_pval")
    colnames(merged_df) <- colname1
    merged_df <- left_join(merged_df, data_F, by = c("exposure"))
    
    file_name <- gsub(".gz","",strsplit(file_name, "/")[[1]][5])
    save_name1 <- paste0(save_path,file_name,"_MR.RData")
    save_name2 <- paste0(save_path,file_name, "_MR.csv")
    save(mydata, file = save_name1)
    write.csv(merged_df, save_name2, row.names = FALSE)
    gc()
  }
}

pQTLs_Sum <- read_csv("pQTLs_Sum.csv")
Finngen_MR("A:/Biology data/FinnGen R10 summary data/data/","phewas4/",df)
##### spearman_cor     #####
df <- read_excel("D:/Desktop/disposal data.xlsx", sheet = "Sheet2")
df <- as.data.frame(df)

calculate_spearman_bootstrap <- function(df, var1, var2) {
  df_clean <- na.omit(df[, c(var1, var2)])
  cor_test_result <- cor.test(df_clean[[var1]], df_clean[[var2]], method = "spearman")
  spearman_cor <- cor_test_result$estimate
  p_value <- cor_test_result$p.value
  bootstrap_spearman <- function(data, indices) {
    boot_data <- data[indices, ]
    cor(boot_data[[var1]], boot_data[[var2]], method = "spearman")
  }
  
  n_boot <- 1000
  boot_results <- boot(data = df_clean, statistic = bootstrap_spearman, R = n_boot)
  boot_ci <- boot.ci(boot_results, type = "perc") 
  cat(paste(var2, "rho:", spearman_cor, "pval:", p_value, "95_ci:", boot_ci$percent[4], "-", boot_ci$percent[5], "\n"))
  cat("\n")
}

variables <- c(names(df)[3],names(df)[9:33])
lapply(variables, function(i) calculate_spearman_bootstrap(df, "NCAN;Neurocan core protein", i))
##### KM plot          #####
df <- read_excel("disposal data0107.xlsx", sheet = "Sheet3")
df <- as.data.frame(df)
df <- df[, c(1, 2, 30, 31, 32, 33, 34, 36, 58, 59)]
df <- subset(df, `Any_liver_disease_time` > 0)
df$`Time_of_death` <- df$`Time_of_death` / 365.25
df$`Any_liver_disease_time` <- df$`Any_liver_disease_time` / 365.25

quartiles <- quantile(df$`NCAN;Neurocan core protein`, probs = seq(0, 1, 1/4), na.rm = TRUE)
df$NCAN_quartile <- cut(df$`NCAN;Neurocan core protein`,
                        breaks = quartiles, 
                        include.lowest = TRUE, 
                        labels = FALSE)

surv_obj <- Surv(df$Any_liver_disease_time, df$Any_liver_disease)
surv_obj2 <- Surv(df$Time_of_death, df$`Cause of liver death: ICD10 | Instance 0`)

fit <- survfit(surv_obj ~ NCAN_quartile, data = df)
fit2 <- survfit(surv_obj2 ~ NCAN_quartile, data = df)

deg_point_fill <- brewer.pal(length(unique(df$NCAN_quartile)), "RdYlBu")

ggsurvplot(fit, data = df, 
           pval = TRUE,
           risk.table = TRUE,
           palette = deg_point_fill,
           xlab = "Time (years)", 
           ylab = "Incidence probability (%)",
           ylim = c(0.9, 1),
           xlim = c(0, 15),
           break.time.by = 3,
           legend = "bottom"
)
ggsurvplot(fit2, data = df, 
           pval = TRUE,
           risk.table = TRUE,
           palette = deg_point_fill,
           xlab = "Time (years)", 
           ylab = "Survival probability (%)",
           ylim = c(0.99, 1),
           xlim = c(0, 15),
           break.time.by = 3,
           legend = "bottom"
)
##### enrichment plot  #####
x <- c('AGER', 'ATF6B', 'AKR1C3', 'APOC1', 'APOE', 'ASL',
       'BRD2', 'CCL19', 'F3', 'GCKR', 'GPN1', 'HLA-DRA',
       'IL18R1', 'LTA', 'MASP1', 'NCAN', 'NRBP1', 'NUDT9', 'SHMT1')

test = bitr(x, fromType="SYMBOL", toType="ENTREZID",OrgDb="org.Hs.eg.db")

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)

ego_ALL <- enrichGO(gene = test$ENTREZID, 
                    universe = names(geneList),
                    OrgDb = org.Hs.eg.db,
                    #keytype = 'ENSEMBL',
                    ont = "ALL",
                    pAdjustMethod = "bonferroni", 
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.05,
                    readable = TRUE) 
write.csv(ego_ALL,'Go_enrich.csv', row.names = F)

ego_20 <- ego_ALL[1:20,]
ego_20$pathway <- factor(ego_20$Description, levels = rev(ego_20$Description))

mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11), 
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), 
                 legend.title = element_text(size = 13), 
                 legend.text = element_text(size = 11))

ego_20$pathway <- sapply(ego_20$pathway, function(x) {
  str_wrap(x, width = 50)
})

ggplot(data = ego_20, 
       aes(x = Count, y = pathway, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.8) + 
  scale_fill_distiller(palette = 'Blues', direction = 1) +
  labs(x = "Number of Gene",
       title = "Gene Ontology enrichment analysis") +
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 1))

##### scRNA-seq plot   #####

datasets <- c("GSE136103", "GSE162616", "GSE174748", "GSE212047", "GSE202379")

for (dataset in datasets) {
  zip_path <- paste0("./", dataset, "/filtered_feature_bc_matrix.zip")
  if (file.exists(zip_path)) {
    unzip(zipfile = zip_path, exdir = paste0("./", dataset, "/"))
  }
}
seurat_list <- list()
for (dataset in datasets) {
  data_dir <- paste0("./", dataset, "/filtered_feature_bc_matrix/")
  counts <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = dataset,
    min.cells = 3,
    min.features = 200
  )
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_list[[dataset]] <- seurat_obj
}

merged_seurat <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[2:5],
  add.cell.ids = datasets
)

pdf("combined_qc.pdf", width = 10, height = 6)
VlnPlot(merged_seurat, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, 
        group.by = "orig.ident",
        pt.size = 0.1) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

merged_seurat <- subset(
  merged_seurat,
  subset = 
    nFeature_RNA > 200 & 
    nFeature_RNA < 7000 &
    nCount_RNA < 20000 &
    percent.mt < 20
)

scedata <- merged_seurat
scedata <- NormalizeData(scedata, 
                         normalization.method = "LogNormalize", 
                         scale.factor = 10000)
scedata <- FindVariableFeatures(scedata, 
                                selection.method = "vst",
                                nfeatures = 2000)
scedata <- ScaleData(scedata,verbose = T) 
scedata <- RunPCA(scedata, npcs = 50, verbose = T)

pdf("pca.pdf",width = 6,height = 6)
VizDimLoadings(scedata, dims = 1:2, reduction = "pca")
DimPlot(scedata, reduction = "pca")
DimHeatmap(scedata, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()

scedata <- JackStraw(scedata, num.replicate = 100,dims = 50)
scedata <- ScoreJackStraw(scedata, dims = 1:50)
pdf(paste0("JACKStrawplots.pdf"),width = 6,height = 5)
JackStrawPlot(object = scedata, dims = 1:50)
dev.off()


pdf(paste0("PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(scedata,ndims = 50)
dev.off()

dim.use = 1:30
scedata <- FindNeighbors(scedata, dims = dim.use, reduction = "pca")
scedata <- FindClusters(scedata, resolution = 0.5, cluster.name = "unintegrated_clusters")
scedata <- RunUMAP(scedata, dims = dim.use, 
                   reduction = "pca", 
                   reduction.name = "umap.unintegrated")
Idents(scedata) = scedata$orig.ident
DimPlot(scedata, reduction = "umap.harmony")


#Harmony
scedata <- IntegrateLayers(
  object = scedata, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = T
)
#####Harmony
scedata <- FindNeighbors(scedata, reduction = "harmony", dims = dim.use)
scedata <- FindClusters(scedata, resolution = 0.2)

scedata <- RunUMAP(scedata, reduction = "harmony", 
                   dims = dim.use, 
                   reduction.name = "umap.harmony")

scedata <- JoinLayers(scedata)
scedata
save(scedata,file = "scedata.rdata")

Idents(scedata) = scedata$RNA_snn_res.0.2


markers <- c("IGKC","IGLC2","IGHG4","IGHG3","IGLC3",
             "IGHG1","JCHAIN","IGHG2","BICC1",
             "ANXA4","PKHD1","TM4SF4","KCNQ1OT1",
             "KRT18","KRT8","CD24","AKAP12","TM4SF1","IFI27",
             "DNASE1L3","SPARCL1","IFITM3","FCN3","MGP",
             "ALB","APOA2","APOC3","SERPINA1","CYP3A5",
             "APOA1","CYP2B6","APOC1","HLA-DRA","CD74","FTL","HLA-DPB1",
             "HLA-DPA1","HLA-DRB1","C1QB","FTH1","TUBA1B","STMN1",
             "HIST1H4C","TUBB","GAPDH","HMGN2","H2AFZ","HMGB1","ACTA2",
             "TAGLN","MYL9","CALD1","TPM2","ADIRF","DSTN","CCL5",
             "CXCR4","CD69","CCL4","NKG7","B2M","KLRB1","TMSB4X"
)




color_gradient <- scale_color_gradientn(
  colors = c("#000000", "blue", "#ffac00", "#ff0000"),
  values = rescale(c(0, 0.25, 0.5, 1))
)

dotplot = DotPlot(scedata, features = markers,dot.min = 0.01) +
  coord_flip() + theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c("lightgrey", "firebrick3")) + 
  labs(x = NULL, y = NULL) + guides(size = guide_legend(order = 3))+
  color_gradient +
  geom_vline(xintercept = c(7.5,15.5,23.5,31.5,39.5,47.5,54.5), linetype = "dashed", color = "grey40")  
pdf(("marker_named.pdf"), width = 6, height = 8)
print(dotplot)
dev.off()



new.cluster.ids <- c("0"="Hep", 
                     "1"="Hep", 
                     "2"="Hep", 
                     "3"="Endo", 
                     "4"="Hep", 
                     "5"="Prol", 
                     "6"="Myel",
                     "7"="Prol",
                     "8"="Hep",
                     "9"="T",
                     "10"="Hep",
                     "11"="Stell",
                     "12"="Chol",
                     "13"="Endo",
                     "14"="Myel",
                     "15"="Prol")


scedata <- RenameIdents(scedata, new.cluster.ids)                        
scedata$celltype <- scedata@active.ident


pvn <- VlnPlot(scedata, features = c("NCAN","NUDT9", "ADH1B", "APOC1", "APOH", "GCKR"),pt.size = 0)
pdf(paste0("vln.pdf"),width = 10,height = 6)
print(pvn)
dev.off()

pdf(("dimpolt.pdf"), width = 6, height = 5)
DimPlot(scedata,reduction = "umap.harmony", label = T,repel = T,label.size = 3)
dev.off()
##### timeROC 5fold UK #####
df <- read_excel("disposal data0107.xlsx", sheet = "Sheet3")
df <- as.data.frame(df)
df <- df[,c(1,2,30,31,32,33,34,36,39,40,58,59)]
names(df)
df <- subset(df, `Any_liver_disease_time` > 0)
df$`Time_of_death` <- df$`Time_of_death` / 365.25
df$`Time of HCC diagnosis | Instance 0` <- df$`Time of HCC diagnosis | Instance 0` / 365.25
df$`Any_liver_disease_time` <- df$`Any_liver_disease_time` / 365.25 ; max(df$`Time of HCC diagnosis | Instance 0`)
time_points <- seq(1, 16)

set.seed(123) 

cross_validate_roc <- function(df, time_points) {
  folds <- createFolds(df$`Type of HCC: ICD10 | Instance 0`, k = 5, list = TRUE)
  
  roc_results <- list()
  
  for (fold in folds) {
    train_data <- df[-fold, ]
    test_data <- df[fold, ]
    
    cox_model_NCAN <- coxph(Surv(`Any_liver_disease_time`, `Any_liver_disease`) ~ `NCAN;Neurocan core protein`, data=train_data)
    cox_model_NCAN_FIB4 <- coxph(Surv(`Any_liver_disease_time`, `Any_liver_disease`) ~ `NCAN;Neurocan core protein` + `FIB-4`, data=train_data)
    cox_model_NCAN_Liverrisk <- coxph(Surv(`Any_liver_disease_time`, `Any_liver_disease`) ~ `NCAN;Neurocan core protein` + `liverRisk score`, data=train_data)
    cox_model_NCAN_aMAP <- coxph(Surv(`Any_liver_disease_time`, `Any_liver_disease`) ~ `NCAN;Neurocan core protein` + `aMAP`, data=train_data)
    cox_model_NCAN_mPAGEB <- coxph(Surv(`Any_liver_disease_time`, `Any_liver_disease`) ~ `NCAN;Neurocan core protein` + `mPAGE-B`, data=train_data)
    
    test_data$NCAN <- predict(cox_model_NCAN, newdata=test_data, type="risk")
    test_data$NCAN_FIB4 <- predict(cox_model_NCAN_FIB4, newdata=test_data, type="risk")
    test_data$NCAN_Liverrisk <- predict(cox_model_NCAN_Liverrisk, newdata=test_data, type="risk")
    test_data$NCAN_aMAP <- predict(cox_model_NCAN_aMAP, newdata=test_data, type="risk")
    test_data$NCAN_mPAGEB <- predict(cox_model_NCAN_mPAGEB, newdata=test_data, type="risk")
    
    roc_list <- list()
    
    for (marker in c("NCAN", "FIB-4", "NCAN_FIB4", "liverRisk score", "NCAN_Liverrisk", "aMAP", "NCAN_aMAP", "mPAGE-B", "NCAN_mPAGEB")) {
      roc_res <- timeROC(T=test_data$`Any_liver_disease_time`,
                         delta=test_data$`Any_liver_disease`,
                         marker=test_data[[marker]],
                         cause=1,
                         times=time_points,
                         iid=FALSE)
      roc_list[[marker]] <- roc_res$AUC
    }
    
    roc_results[[length(roc_results) + 1]] <- roc_list
  }
  return(roc_results)
}

auroc_values <- cross_validate_roc(df, time_points)

final_df <- data.frame(Time = integer(), fold = integer(), auc = numeric(), Model = character())

for (fold_idx in 1:length(auroc_values)) {
  fold_data <- auroc_values[[fold_idx]]
  for (model_idx in 1:length(fold_data)) {
    model_data <- fold_data[[model_idx]]
    temp_df <- data.frame(
      Time = 1:length(model_data),  
      fold = fold_idx,              
      auc = model_data,             
      Model = names(auroc_values[[fold_idx]])[model_idx]
    )
    final_df <- rbind(final_df, temp_df)
  }
}

df <- final_df

df_stats <- df %>%
  group_by(Time, Model) %>%
  summarise(
    mean_auc = mean(auc),
    sd_auc = sd(auc),
    min_auc = min(auc),
    max_auc = max(auc),
    .groups = 'drop'
  )

model_order <- c("NCAN", "FIB-4", "NCAN_FIB4", "liverRisk score", 
                 "NCAN_Liverrisk", "aMAP", "NCAN_aMAP", "mPAGE-B", "NCAN_mPAGEB")
df$Model <- factor(df$Model, levels = model_order)
df_stats$Model <- factor(df_stats$Model, levels = model_order)

df <- df %>%
  mutate(Model = case_when(
    Model == "NCAN_FIB4" ~ "NCAN_FIB-4",
    Model == "liverRisk score" ~ "LiverRisk score",
    Model == "NCAN_Liverrisk" ~ "NCAN_LiverRisk score",
    Model == "NCAN_mPAGEB" ~ "NCAN_mPAGE-B",
    TRUE ~ Model
  ))
df_stats <- df_stats %>%
  mutate(Model = case_when(
    Model == "NCAN_FIB4" ~ "NCAN_FIB-4",
    Model == "liverRisk score" ~ "LiverRisk score",
    Model == "NCAN_Liverrisk" ~ "NCAN_LiverRisk score",
    Model == "NCAN_mPAGEB" ~ "NCAN_mPAGE-B",
    TRUE ~ Model
  ))

model_order <- c("NCAN", "FIB-4", "NCAN_FIB-4", "LiverRisk score", 
                 "NCAN_LiverRisk score", "aMAP", "NCAN_aMAP", "mPAGE-B", "NCAN_mPAGE-B")
df$Model <- factor(df$Model, levels = model_order)
df_stats$Model <- factor(df_stats$Model, levels = model_order)

colors <- brewer.pal(9, "RdYlBu")

ggplot() +
  geom_point(data = df, aes(x = factor(Time), y = auc, color = Model), 
             size = 2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0, dodge.width = 1)) +
  geom_errorbar(data = df_stats, 
                aes(x = factor(Time), ymin = min_auc, ymax = max_auc, color = Model), 
                width = 0.4, position = position_dodge(1)) +
  geom_rect(data = df_stats,
            aes(xmin = as.numeric(factor(Time)) - 0.35, 
                xmax = as.numeric(factor(Time)) + 0.35,
                ymin = mean_auc - sd_auc,
                ymax = mean_auc + sd_auc,
                fill = Model),
            alpha = 0.4, position = position_dodge(1)) +
  geom_point(data = df_stats, 
             aes(x = factor(Time), y = mean_auc, fill = Model), 
             shape = 21, size = 3.5, 
             position = position_dodge(1)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) + 
  theme_minimal() +
  labs(x = "Time-point (Years)", y = "AUC", color = "Model", fill = "Model") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white")
  )

write.csv(df_stats,'5 fold_auc_ald.csv', row.names = F)

##### ROC Chinese      #####
set.seed(125)

data <- read_excel("5f_roc.xlsx") 
data <- subset(data, !is.na(data$FIB)) # HCC
label <- data$FIB
data <- data[, c("NCAN", "FIB4", "LiverRisk_score", "aMAP", "mPAGEB")]

model_names <- c("NCAN", "FIB4", "NCAN_FIB4",
                 "LiverRisk_score", "NCAN_LiverRisk_score")

# model_names <- c("NCAN", "aMAP", "NCAN_aMAP",
#                  "mPAGEB", "NCAN_mPAGEB")

folds <- createFolds(label, k = 5, list = TRUE, returnTrain = TRUE)
tprs <- list()
aucs <- list()
model_metrics <- list()
model_predictions <- list()
mean_fpr <- seq(0, 1, length.out = 100)

calculate_roc <- function(model, X_train, Y_train, X_test, Y_test) {
  model <- glm(Y_train ~ ., data = X_train, family = binomial())
  y_pred <- predict(model, newdata = X_test, type = "response")
  
  roc_curve <- roc(Y_test, y_pred, direction = "<")
  fpr <- rev(1 - roc_curve$specificities)
  tpr <- rev(roc_curve$sensitivities)
  
  tpr_interp <- approx(fpr, tpr, xout = mean_fpr, ties = mean)$y
  tpr_interp[1] <- 0
  
  auc_value <- auc(roc_curve)
  return(list(tpr_interp = tpr_interp, auc_value = auc_value, y_pred = y_pred, Y_test = Y_test))
}

calculate_metrics <- function(Y_true, Y_pred) {
  threshold <- 0.5
  predicted_labels <- ifelse(Y_pred > threshold, 1, 0)
  
  tp <- sum((Y_true == 1) & (predicted_labels == 1))
  tn <- sum((Y_true == 0) & (predicted_labels == 0))
  fp <- sum((Y_true == 0) & (predicted_labels == 1))
  fn <- sum((Y_true == 1) & (predicted_labels == 0))
  
  sensitivity <- tp / (tp + fn) 
  specificity <- tn / (tn + fp) 
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- tp / (tp + fp) 
  f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  youden_index <- sensitivity + specificity - 1
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  
  return(list(
    sensitivity = sensitivity,
    specificity = specificity,
    accuracy = accuracy,
    precision = precision,
    f1_score = f1_score,
    youden_index = youden_index,
    ppv = ppv,
    npv = npv
  ))
}

calculate_Yoden_metrics <- function(Y_true, Y_pred) {
  best_threshold <- 0
  best_youden_index <- -1
  
  for (threshold in seq(0, 1, by = 0.01)) {
    predicted_labels <- ifelse(Y_pred > threshold, 1, 0)
    
    tp <- sum((Y_true == 1) & (predicted_labels == 1))
    tn <- sum((Y_true == 0) & (predicted_labels == 0))
    fp <- sum((Y_true == 0) & (predicted_labels == 1))
    fn <- sum((Y_true == 1) & (predicted_labels == 0))
    
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    youden_index <- sensitivity + specificity - 1
    
    if (youden_index > best_youden_index) {
      best_youden_index <- youden_index
      best_threshold <- threshold
    }
  }
  
  predicted_labels <- ifelse(Y_pred > best_threshold, 1, 0)
  
  tp <- sum((Y_true == 1) & (predicted_labels == 1))
  tn <- sum((Y_true == 0) & (predicted_labels == 0))
  fp <- sum((Y_true == 0) & (predicted_labels == 1))
  fn <- sum((Y_true == 1) & (predicted_labels == 0))
  
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- tp / (tp + fp)
  f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  
  roc_curve <- roc(Y_true, Y_pred)
  auc_value <- auc(roc_curve)
  
  return(list(
    sensitivity = sensitivity,
    specificity = specificity,
    accuracy = accuracy,
    precision = precision,
    f1_score = f1_score,
    youden_index = best_youden_index,
    ppv = ppv,
    npv = npv,
    auc = auc_value,
    best_cut_off = best_threshold
  ))
}

for (model_name in model_names) {
  selected_features <- switch(model_name,
                              "NCAN" = c("NCAN"),
                              "FIB4" = c("FIB4"),
                              "NCAN_FIB4" = c("NCAN", "FIB4"),
                              "LiverRisk_score" = c("LiverRisk_score"),
                              "NCAN_LiverRisk_score" = c("NCAN", "LiverRisk_score"),
                              
                              "aMAP" = c("aMAP"),
                              "NCAN_aMAP" = c("NCAN", "aMAP"),
                              "mPAGEB" = c("mPAGEB"),
                              "NCAN_mPAGEB" = c("NCAN", "mPAGEB")
  )
  
  model_tprs <- list()
  model_aucs <- c()
  all_predicted_probs <- c()
  all_true_labels <- c()
  model_predictions[[model_name]] <- list(y_pred = c(), Y_test = c())
  
  for (i in 1:length(folds)) {
    train_index <- folds[[i]]
    test_index <- setdiff(1:length(label), train_index)
    
    X_train <- data[train_index, selected_features]
    Y_train <- label[train_index]
    X_test <- data[test_index, selected_features]
    Y_test <- label[test_index]
    
    roc_results <- calculate_roc(NULL, X_train, Y_train, X_test, Y_test)
    model_aucs <- c(model_aucs, roc_results$auc_value)
    model_tprs[[i]] <- roc_results$tpr_interp
    
    all_predicted_probs <- c(all_predicted_probs, roc_results$y_pred)
    all_true_labels <- c(all_true_labels, Y_test)
    
    model_predictions[[model_name]]$y_pred <- c(model_predictions[[model_name]]$y_pred, roc_results$y_pred)
    model_predictions[[model_name]]$Y_test <- c(model_predictions[[model_name]]$Y_test, Y_test)
  }
  
  mean_auc <- mean(model_aucs)
  std_auc <- sd(model_aucs)
  
  mean_tpr <- rowMeans(do.call(cbind, model_tprs))
  std_tpr <- apply(do.call(cbind, model_tprs), 1, sd)
  
  aucs[[model_name]] <- mean_auc
  tprs[[model_name]] <- list(mean_tpr = mean_tpr, std_tpr = std_tpr)
  
  cat(paste0(model_name, "AUC: ", round(mean_auc, 4), "\n"))
  metrics <- calculate_Yoden_metrics(all_true_labels, all_predicted_probs)

  model_metrics[[model_name]] <- c(mean_auc, metrics$sensitivity, metrics$specificity, metrics$accuracy, 
                                   metrics$precision, metrics$f1_score, metrics$youden_index, metrics$ppv, metrics$npv)
  model_metrics_df <- do.call(rbind, lapply(model_metrics, function(x) as.data.frame(t(x))))
  colnames(model_metrics_df) <- c("AUC", "Sensitivity", "Specificity", "Accuracy", "Precision", "F1_score", "Youden_index", "PPV", "NPV")
  model_metrics_df <- cbind(Model = rownames(model_metrics_df), model_metrics_df)
  row.names(model_metrics_df) <- NULL
}

print(model_metrics_df)

plot_data <- data.frame(
  mean_fpr = rep(mean_fpr, length(model_names)),
  mean_tpr = unlist(lapply(tprs, function(x) x$mean_tpr)),
  model_name = rep(model_names, each = length(mean_fpr))
)

model_colors <- c("#D7191CFF", "#FDAE61FF", "#ABD9E9FF", "#2C7BB6FF", "#313695")

plot_data$model_name <- factor(plot_data$model_name, levels = model_names)
p <- ggplot(plot_data, aes(x = mean_fpr, y = mean_tpr, color = model_name, group = model_name)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = model_colors) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "False Positive Rate", y = "True Positive Rate", title = "HCC ROC") +
  theme_minimal() +
  theme(
    legend.position = "bottomright", 
    text = element_text(size = 12), 
    plot.title = element_text(hjust = 0.5), 
    plot.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key = element_blank()
  )

for (i in 1:length(model_names)) {
  p <- p + annotate("text", x = 0.5, y = 0.35 - (i * 0.06), 
                    label = paste0(model_names[i], ": AUC = ", round(aucs[[model_names[i]]], 3)), 
                    color = model_colors[i], hjust = 0, alpha = 1)
}

p