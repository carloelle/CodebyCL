#initial tissue classification for ST data. same as report (if you want, look at zenodo link)

GSM7058757_C2$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',GSM7058757_C2$Stromal_C2L_percentage>GSM7058757_C2$Tumoral_C2L_percentage)))
GSM7058756_C1$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',GSM7058756_C1$Stromal_C2L_percentage>GSM7058756_C1$Tumoral_C2L_percentage)))
GSM7058758_C3$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',GSM7058758_C3$Stromal_C2L_percentage>GSM7058758_C3$Tumoral_C2L_percentage)))
GSM7058759_C4$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',GSM7058759_C4$Stromal_C2L_percentage>GSM7058759_C4$Tumoral_C2L_percentage)))

V573_1$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V573_1$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V573_1$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V573_1$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V573_1$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V573_1$CMS1_KoreanDeconv_percentageKorean+
                                                                    V573_1$CMS2_KoreanDeconv_percentageKorean+
                                                                    V573_1$CMS3_KoreanDeconv_percentageKorean+
                                                                    V573_1$CMS4_KoreanDeconv_percentageKorean))))
V573_2$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V573_2$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V573_2$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V573_2$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V573_2$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V573_2$CMS1_KoreanDeconv_percentageKorean+
                                                                    V573_2$CMS2_KoreanDeconv_percentageKorean+
                                                                    V573_2$CMS3_KoreanDeconv_percentageKorean+
                                                                    V573_2$CMS4_KoreanDeconv_percentageKorean))))
V371_1$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V371_1$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V371_1$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V371_1$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V371_1$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V371_1$CMS1_KoreanDeconv_percentageKorean+
                                                                    V371_1$CMS2_KoreanDeconv_percentageKorean+
                                                                    V371_1$CMS3_KoreanDeconv_percentageKorean+
                                                                    V371_1$CMS4_KoreanDeconv_percentageKorean))))
V371_2$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V371_2$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V371_2$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V371_2$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V371_2$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V371_2$CMS1_KoreanDeconv_percentageKorean+
                                                                    V371_2$CMS2_KoreanDeconv_percentageKorean+
                                                                    V371_2$CMS3_KoreanDeconv_percentageKorean+
                                                                    V371_2$CMS4_KoreanDeconv_percentageKorean))))
V838_1$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V838_1$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V838_1$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V838_1$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V838_1$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V838_1$CMS1_KoreanDeconv_percentageKorean+
                                                                    V838_1$CMS2_KoreanDeconv_percentageKorean+
                                                                    V838_1$CMS3_KoreanDeconv_percentageKorean+
                                                                    V838_1$CMS4_KoreanDeconv_percentageKorean))))
V838_2$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V838_2$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V838_2$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V838_2$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V838_2$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V838_2$CMS1_KoreanDeconv_percentageKorean+
                                                                    V838_2$CMS2_KoreanDeconv_percentageKorean+
                                                                    V838_2$CMS3_KoreanDeconv_percentageKorean+
                                                                    V838_2$CMS4_KoreanDeconv_percentageKorean))))
V763_1$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V763_1$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V763_1$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V763_1$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V763_1$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V763_1$CMS1_KoreanDeconv_percentageKorean+
                                                                    V763_1$CMS2_KoreanDeconv_percentageKorean+
                                                                    V763_1$CMS3_KoreanDeconv_percentageKorean+
                                                                    V763_1$CMS4_KoreanDeconv_percentageKorean))))
V763_2$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V763_2$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V763_2$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V763_2$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V763_2$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V763_2$CMS1_KoreanDeconv_percentageKorean+
                                                                    V763_2$CMS2_KoreanDeconv_percentageKorean+
                                                                    V763_2$CMS3_KoreanDeconv_percentageKorean+
                                                                    V763_2$CMS4_KoreanDeconv_percentageKorean))))
V688_1$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V688_1$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V688_1$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V688_1$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V688_1$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V688_1$CMS1_KoreanDeconv_percentageKorean+
                                                                    V688_1$CMS2_KoreanDeconv_percentageKorean+
                                                                    V688_1$CMS3_KoreanDeconv_percentageKorean+
                                                                    V688_1$CMS4_KoreanDeconv_percentageKorean))))
V688_2$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V688_2$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V688_2$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V688_2$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V688_2$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V688_2$CMS1_KoreanDeconv_percentageKorean+
                                                                    V688_2$CMS2_KoreanDeconv_percentageKorean+
                                                                    V688_2$CMS3_KoreanDeconv_percentageKorean+
                                                                    V688_2$CMS4_KoreanDeconv_percentageKorean))))
V015_1$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V015_1$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V015_1$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V015_1$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V015_1$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V015_1$CMS1_KoreanDeconv_percentageKorean+
                                                                    V015_1$CMS2_KoreanDeconv_percentageKorean+
                                                                    V015_1$CMS3_KoreanDeconv_percentageKorean+
                                                                    V015_1$CMS4_KoreanDeconv_percentageKorean))))
V015_2$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V015_2$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V015_2$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V015_2$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V015_2$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V015_2$CMS1_KoreanDeconv_percentageKorean+
                                                                    V015_2$CMS2_KoreanDeconv_percentageKorean+
                                                                    V015_2$CMS3_KoreanDeconv_percentageKorean+
                                                                    V015_2$CMS4_KoreanDeconv_percentageKorean))))
V797_1$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V797_1$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V797_1$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V797_1$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V797_1$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V797_1$CMS1_KoreanDeconv_percentageKorean+
                                                                    V797_1$CMS2_KoreanDeconv_percentageKorean+
                                                                    V797_1$CMS3_KoreanDeconv_percentageKorean+
                                                                    V797_1$CMS4_KoreanDeconv_percentageKorean))))
V797_2$Tissue_Classification=factor(gsub('FALSE','Tumoral',gsub('TRUE','Stromal',(V797_2$Stromal.1_KoreanDeconv_percentageKorean+
                                                                                    V797_2$Stromal.2_KoreanDeconv_percentageKorean+
                                                                                    V797_2$Stromal.3_KoreanDeconv_percentageKorean+
                                                                                    V797_2$Myofibroblasts_KoreanDeconv_percentageKorean)
                                                                >(V797_2$CMS1_KoreanDeconv_percentageKorean+
                                                                    V797_2$CMS2_KoreanDeconv_percentageKorean+
                                                                    V797_2$CMS3_KoreanDeconv_percentageKorean+
                                                                    V797_2$CMS4_KoreanDeconv_percentageKorean))))
Idents(GSM7058756_C1)='Tissue_Classification'
Idents(GSM7058757_C2)='Tissue_Classification'
Idents(GSM7058758_C3)='Tissue_Classification'
Idents(GSM7058759_C4)='Tissue_Classification'
Idents(V573_1)='Tissue_Classification'
Idents(V573_2)='Tissue_Classification'
Idents(V371_1)='Tissue_Classification'
Idents(V371_2)='Tissue_Classification'
Idents(V838_1)='Tissue_Classification'
Idents(V838_2)='Tissue_Classification'
Idents(V763_1)='Tissue_Classification'
Idents(V763_2)='Tissue_Classification'
Idents(V688_1)='Tissue_Classification'
Idents(V688_2)='Tissue_Classification'
Idents(V015_1)='Tissue_Classification'
Idents(V015_2)='Tissue_Classification'
Idents(V797_1)='Tissue_Classification'
Idents(V797_2)='Tissue_Classification'

V015_1<-NormalizeData(V015_1)
V015_2<-NormalizeData(V015_2)
V371_1<-NormalizeData(V371_1)
V371_2<-NormalizeData(V371_2)
V573_1<-NormalizeData(V573_1)
V573_2<-NormalizeData(V573_2)
V688_1<-NormalizeData(V688_1)
V688_2<-NormalizeData(V688_2)
V763_1<-NormalizeData(V763_1)
V763_2<-NormalizeData(V763_2)
V797_1<-NormalizeData(V797_1)
V797_2<-NormalizeData(V797_2)
V838_1<-NormalizeData(V838_1)
V838_2<-NormalizeData(V838_2)

rownames(V015_1@assays$Spatial@features)->rownames(V015_1@assays$Spatial@layers$counts)
rownames(V015_1@assays$Spatial@features)->rownames(V015_1@assays$Spatial@layers$data)
rownames(V015_2@assays$Spatial@features)->rownames(V015_2@assays$Spatial@layers$counts)
rownames(V015_2@assays$Spatial@features)->rownames(V015_2@assays$Spatial@layers$data)
rownames(V371_1@assays$Spatial@features)->rownames(V371_1@assays$Spatial@layers$counts)
rownames(V371_1@assays$Spatial@features)->rownames(V371_1@assays$Spatial@layers$data)
rownames(V371_2@assays$Spatial@features)->rownames(V371_2@assays$Spatial@layers$counts)
rownames(V371_2@assays$Spatial@features)->rownames(V371_2@assays$Spatial@layers$data)
rownames(V573_1@assays$Spatial@features)->rownames(V573_1@assays$Spatial@layers$counts)
rownames(V573_1@assays$Spatial@features)->rownames(V573_1@assays$Spatial@layers$data)
rownames(V573_2@assays$Spatial@features)->rownames(V573_2@assays$Spatial@layers$counts)
rownames(V573_2@assays$Spatial@features)->rownames(V573_2@assays$Spatial@layers$data)
rownames(V688_1@assays$Spatial@features)->rownames(V688_1@assays$Spatial@layers$counts)
rownames(V688_1@assays$Spatial@features)->rownames(V688_1@assays$Spatial@layers$data)
rownames(V688_2@assays$Spatial@features)->rownames(V688_2@assays$Spatial@layers$counts)
rownames(V688_2@assays$Spatial@features)->rownames(V688_2@assays$Spatial@layers$data)
rownames(V763_1@assays$Spatial@features)->rownames(V763_1@assays$Spatial@layers$counts)
rownames(V763_1@assays$Spatial@features)->rownames(V763_1@assays$Spatial@layers$data)
rownames(V763_2@assays$Spatial@features)->rownames(V763_2@assays$Spatial@layers$counts)
rownames(V763_2@assays$Spatial@features)->rownames(V763_2@assays$Spatial@layers$data)
rownames(V797_1@assays$Spatial@features)->rownames(V797_1@assays$Spatial@layers$counts)
rownames(V797_1@assays$Spatial@features)->rownames(V797_1@assays$Spatial@layers$data)
rownames(V797_2@assays$Spatial@features)->rownames(V797_2@assays$Spatial@layers$counts)
rownames(V797_2@assays$Spatial@features)->rownames(V797_2@assays$Spatial@layers$data)
rownames(V838_1@assays$Spatial@features)->rownames(V838_1@assays$Spatial@layers$counts)
rownames(V838_1@assays$Spatial@features)->rownames(V838_1@assays$Spatial@layers$data)
rownames(V838_2@assays$Spatial@features)->rownames(V838_2@assays$Spatial@layers$counts)
rownames(V838_2@assays$Spatial@features)->rownames(V838_2@assays$Spatial@layers$data)

sign$mrCAF_Cords<-unique(c(sign$iCAF_Cords,sign$mCAF_Cords,sign$MF1))
sign$mrCAF_sign<-unique(c(sign$iCAF,sign$mCAF,sign$MF1))
sign$iCAF_Cords=NULL
sign$mCAF_Cords=NULL
sign$dCAF_Cords=NULL
sign$apCAF_Cords=NULL
sign$pericyte_Cords=NULL
sign$vCAFrCAF_Cords=NULL
sign$tCAF_Cords=NULL
sign$iCAF=NULL
sign$MF2=NULL
sign$MF3=NULL
sign$MF4=NULL
sign$MF1=NULL
sign$cCAF=NULL
sign$dCAF=NULL
sign$vCAF=NULL
sign$mCAF=NULL
sign$myCAF=NULL
CRIS<-read.csv('/home/carlo/Desktop/transferdata/CRIS/CRIS_sign.csv',header = T)
CRIS<-data.frame(Genes=CRIS$X.1[3:dim(CRIS)[1]],CRIS=CRIS$X.2[3:dim(CRIS)[1]])
sign$CRIS_A<-CRIS%>%filter(CRIS=='CRIS-A')%>%select(Genes)%>%unlist()%>%unname()
sign$CRIS_B<-CRIS%>%filter(CRIS=='CRIS-B')%>%select(Genes)%>%unlist()%>%unname()
sign$CRIS_C<-CRIS%>%filter(CRIS=='CRIS-C')%>%select(Genes)%>%unlist()%>%unname()
sign$CRIS_D<-CRIS%>%filter(CRIS=='CRIS-D')%>%select(Genes)%>%unlist()%>%unname()
sign$CRIS_E<-CRIS%>%filter(CRIS=='CRIS-E')%>%select(Genes)%>%unlist()%>%unname()


library(reticulate)
use_python('/home/carlo/.local/share/r-miniconda/envs/giotto_env/bin/python', required = TRUE)

library(circlize)
library(data.table)
library(limma)

#two version for process_PAGEanalysis, values of PAGE does not change, p-values yes 
#v.2 will correct p-values corrected for empirical bayesian adjustment --> we aim at correcting extreme p-values that might screwup our analysis
#v.1 will just compute p-values based on a permutation test
process_PAGEanalysis2 <- function(seurat_object, signatures_all, only_fibro = TRUE) {
  if (only_fibro) {
    signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]
  }
  raw_exprs <- seurat_object@assays$Spatial@layers$counts
  rownames(raw_exprs) <- rownames(seurat_object@assays$Spatial@features@.Data)
  colnames(raw_exprs) <- rownames(seurat_object@meta.data)
  spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[, 1:2])
  colnames(spatial_locs) <- c("x", "y")
  
  myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs, spatial_locs = spatial_locs)
  myGiottoObj <- normalizeGiotto(gobject = myGiottoObj)
  
  # Create signature matrix for initial PAGE analysis
  all_signatures <- names(signatures_all)
  signature_matrix_complete <- makeSignMatrixPAGE(
    sign_names = all_signatures,
    sign_list = lapply(all_signatures, function(sign) signatures_all[[sign]])
  )
  
  # Run initial PAGE enrichment analysis on all signatures
  myGiottoObj_initial <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 5,
    output_enrichment = c("original")
  )
  
  # Run PAGE enrichment analysis with p-values
  myGiottoObj_pval <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 5,
    p_value = TRUE,
    reverse_log_scale = FALSE,
    return_gobject = FALSE
  )
  
  # Extract p-value results
  reshaped_pval <- myGiottoObj_pval$DT %>%
    select(cell_ID, cell_type, pval) %>%
    pivot_wider(names_from = cell_type, values_from = pval, values_fill = 0)
  
  # Apply empirical Bayes adjustment to p-values
  adjust_pvalues <- function(pval_matrix) {
    apply(pval_matrix, 2, function(p) {
      non_na_pvals <- na.omit(p)
      
      # Add a small constant to avoid zero p-values
      non_na_pvals[non_na_pvals == 0] <- 1e-10
      
      # Fit the linear model without log transformation
      if (length(non_na_pvals) > 1) {
        fit <- eBayes(lmFit(non_na_pvals))
        adjusted <- topTable(fit, number = Inf, sort.by = "none")$P.Value
      } else {
        adjusted <- non_na_pvals
      }
      
      # Expand adjusted p-values back to the original length, keeping NA values
      adjusted_full <- rep(NA, length(p))
      adjusted_full[!is.na(p)] <- adjusted
      return(adjusted_full)
    })
  }
  
  pval_matrix <- as.matrix(reshaped_pval[, -1])
  adjusted_pvals <- adjust_pvalues(pval_matrix)
  
  adjusted_pvals_df <- as.data.frame(adjusted_pvals)
  colnames(adjusted_pvals_df) <- paste0(colnames(adjusted_pvals_df), "_pvalPAGE")
  
  # Combine adjusted p-values with metadata
  seurat_object@meta.data <- cbind(seurat_object@meta.data, adjusted_pvals_df)
  
  # Extract z-score results from initial PAGE analysis
  zscore_df <- as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
  zscore_df <- as.data.frame(zscore_df)
  colnames(zscore_df) <- paste0(colnames(zscore_df), "PAGE")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_df)
  
  return(seurat_object)
}
process_PAGEanalysis1 <- function(seurat_object, signatures_all,only_fibro=T) {
  if(only_fibro==T){signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]}
  raw_exprs <- seurat_object@assays$Spatial@layers$counts
  rownames(raw_exprs) <- rownames(seurat_object@assays$Spatial@features@.Data)
  colnames(raw_exprs) <- rownames(seurat_object@meta.data)
  spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[,1:2])
  colnames(spatial_locs) <- c("x", "y")
  
  myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs, spatial_locs = spatial_locs)
  myGiottoObj <- normalizeGiotto(gobject = myGiottoObj)
  
  # Create signature matrix for initial PAGE analysis
  all_signatures <- names(signatures_all)
  signature_matrix_complete <- makeSignMatrixPAGE(
    sign_names = all_signatures,
    sign_list = lapply(all_signatures, function(sign) signatures_all[[sign]])
  )
  
  # Run initial PAGE enrichment analysis on all signatures
  myGiottoObj_initial <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 5,
    output_enrichment = c("original")
  )
  
  # Run PAGE enrichment analysis with p-values
  myGiottoObj_pval <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 5,
    p_value = TRUE,
    reverse_log_scale = FALSE,
    return_gobject = FALSE
  )
  
  # Extract p-value results
  reshaped_pval <- myGiottoObj_pval$DT %>%
    select(cell_ID, cell_type, pval) %>%  
    pivot_wider(names_from = cell_type, values_from = pval, values_fill = 0) 
  
  pzscore_df<-as.data.frame(reshaped_pval)
  colnames(pzscore_df)=paste0(colnames(pzscore_df), "_pvalPAGE")
  seurat_object@meta.data<-cbind(seurat_object@meta.data, pzscore_df)
  
  zscore_df<-as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
  zscore_df <- as.data.frame(zscore_df)
  colnames(zscore_df) <- paste0(colnames(zscore_df), "PAGE")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_df)
  
  return(seurat_object)
}

library(Giotto)
library(tidyr)
V573_1 <- process_PAGEanalysis1(seurat_object = V573_1, signatures_all = sign)
V573_2 <- process_PAGEanalysis1(seurat_object = V573_2, signatures_all = sign)
V371_1 <- process_PAGEanalysis1(seurat_object = V371_1, signatures_all = sign)
V371_2 <- process_PAGEanalysis1(seurat_object = V371_2, signatures_all = sign)
V838_1 <- process_PAGEanalysis1(seurat_object = V838_1, signatures_all = sign)
V838_2 <- process_PAGEanalysis1(seurat_object = V838_2, signatures_all = sign)
V763_1 <- process_PAGEanalysis1(seurat_object = V763_1, signatures_all = sign)
V763_2 <- process_PAGEanalysis1(seurat_object = V763_2, signatures_all = sign)
V688_1 <- process_PAGEanalysis1(seurat_object = V688_1, signatures_all = sign)
V688_2 <- process_PAGEanalysis1(seurat_object = V688_2, signatures_all = sign)
V015_1 <- process_PAGEanalysis1(seurat_object = V015_1, signatures_all = sign)
V015_2 <- process_PAGEanalysis1(seurat_object = V015_2, signatures_all = sign)
V797_1 <- process_PAGEanalysis1(seurat_object = V797_1, signatures_all = sign)
V797_2 <- process_PAGEanalysis1(seurat_object = V797_2, signatures_all = sign)

GSM7058756_C1<-process_PAGEanalysis1(seurat_object = GSM7058756_C1, signatures_all = sign)
GSM7058757_C2<-process_PAGEanalysis1(seurat_object = GSM7058757_C2, signatures_all = sign)
GSM7058758_C3<-process_PAGEanalysis1(seurat_object = GSM7058758_C3, signatures_all = sign)
GSM7058759_C4<-process_PAGEanalysis1(seurat_object = GSM7058759_C4, signatures_all = sign)


# Updated function to calculate both geometric mean (GM) and standard deviation for each signature, btw it is written variance but it calculates sd 
calculate_corrected_scores_with_variance <- function(seurat_obj) {
  for (name in names(sign)) {
    genes <- sign[[name]]
    valid_genes <- genes[genes %in% rownames(seurat_obj@assays$Spatial@layers$data)]
    
    if (length(valid_genes) > 0) {
      # Extract the expression data for the valid genes
      expr_data <- seurat_obj@assays$Spatial@layers$data[valid_genes, ]
      
      # Calculate geometric mean for the valid genes
      mean.exp <- exp(colMeans(expr_data, na.rm = TRUE))
      
      # Calculate sd of the expression data for valid genes
      sd.exp <- apply(expr_data, 2, sd, na.rm = TRUE)
      
      # Store GM and sd in the metadata for each spot
      score_column_gm <- paste0(name, ".score_GM1")
      score_column_sd <- paste0(name, ".score_sd_GM1")
      
      seurat_obj@meta.data[[score_column_gm]] <- mean.exp
      seurat_obj@meta.data[[score_column_sd]] <- sd.exp
    } else {
      warning(paste("No valid genes found for signature:", name))
    }
  }
  
  return(seurat_obj)
}

V573_1<-calculate_corrected_scores_with_variance(V573_1)
V573_2<-calculate_corrected_scores_with_variance(V573_2)
V371_1<-calculate_corrected_scores_with_variance(V371_1)
V371_2<-calculate_corrected_scores_with_variance(V371_2)
V838_1<-calculate_corrected_scores_with_variance(V838_1)
V838_2<-calculate_corrected_scores_with_variance(V838_2)
V763_1<-calculate_corrected_scores_with_variance(V763_1)
V763_2<-calculate_corrected_scores_with_variance(V763_2)
V688_1<-calculate_corrected_scores_with_variance(V688_1)
V688_2<-calculate_corrected_scores_with_variance(V688_2)
V015_1<-calculate_corrected_scores_with_variance(V015_1)
V015_2<-calculate_corrected_scores_with_variance(V015_2)
V797_1<-calculate_corrected_scores_with_variance(V797_1)
V797_2<-calculate_corrected_scores_with_variance(V797_2)
GSM7058756_C1<-calculate_corrected_scores_with_variance(GSM7058756_C1)
GSM7058757_C2<-calculate_corrected_scores_with_variance(GSM7058757_C2)
GSM7058758_C3<-calculate_corrected_scores_with_variance(GSM7058758_C3)
GSM7058759_C4<-calculate_corrected_scores_with_variance(GSM7058759_C4)

library(tidyr)
setwd('/home/carlo/Desktop/Deltas')

process_seurat_object<-function(seurat_object,obj_name){
  
  metadata <- seurat_object@meta.data
  
  # Load necessary libraries
  library(Seurat)
  library(dplyr)
  
  sign_names <- names(sign)
  
  # Define the prefixes
  valid_categories <- sort(paste0(sign_names, "PAGE"))
  valid_pvals <- sort(paste0(sign_names, "_pvalPAGE"))
  valid_gm <- sort(paste0(sign_names, ".score_GM1"))
  valid_gm_var <- sort(paste0(sign_names, ".score_sd_GM1"))

  # Define the names of PAGE, pval, GM, and variance columns
  # These are the columns to use for each CAF subtype
  #valid_categories <-sort(c("EMRPAGE","MF2PAGE","MF4PAGE","cCAFPAGE","mrCAFPAGE","myCAFPAGE","vCAFPAGE"))  
  #valid_pvals <- sort(c("mrCAF_pvalPAGE","myCAF_pvalPAGE","vCAF_pvalPAGE","MF2_pvalPAGE","cCAF_pvalPAGE","MF4_pvalPAGE","EMR_pvalPAGE"))
  #valid_gm <- sort(c("EMR.score_GM1", "MF2.score_GM1","MF4.score_GM1", "cCAF.score_GM1", "mrCAF.score_GM1", "myCAF.score_GM1", "vCAF.score_GM1"))
  #valid_gm_var <- sort(c("EMR.score_sd_GM1", "MF2.score_sd_GM1","MF4.score_sd_GM1","cCAF.score_sd_GM1", "mrCAF.score_sd_GM1", "myCAF.score_sd_GM1", "vCAF.score_sd_GM1"))


  # Filter the metadata for only valid categories
  metadata_filtered <- metadata %>%
  select(all_of(valid_categories), all_of(valid_pvals), all_of(valid_gm), all_of(valid_gm_var))

  # Normalize PAGE values using Min-Max normalization (for valid categories)
  normalize_min_max <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

  for (cat in valid_categories) {
  metadata_filtered[[paste0(cat, "_norm")]] <- normalize_min_max(as.numeric(metadata_filtered[[cat]]))  # Ensure numeric values
  }

  # Define the function to calculate the weighted score for a given subtype
  calculate_weighted_score_per_subtype <- function(row, category, pval_col, gm_col, gm_var_col) {
    # Check if the required columns are in the row
    if (!(pval_col %in% names(row) && gm_col %in% names(row) && gm_var_col %in% names(row) && paste0(category, "_norm") %in% names(row))) {
      return(NA_real_)
    }
    
    # Calculate the log10 of the p-value for the current subtype
    log_pval <- if (!is.na(as.numeric(row[[pval_col]]))) {
      pval <- as.numeric(row[[pval_col]])
      # Add small constant if p-value is 0
      if (pval == 0) {
        pval <- 1e-07
      }
      -log(pval)  # Use log for the final result
    } else {
      NA_real_  # If p-value is missing, return NA
    }
  
  # Normalize PAGE value
  normalized_page <- if (!is.na(row[[paste0(category, "_norm")]])) {
    as.numeric(row[[paste0(category, "_norm")]])
  } else {
    NA_real_
  }
  
  # Get the GM and GM sd values
  gm_score <- if (!is.na(row[[gm_col]])) {
    as.numeric(row[[gm_col]])
  } else {
    NA_real_
  }
  gm_var <- if (!is.na(row[[gm_var_col]])) {
    as.numeric(row[[gm_var_col]])
  } else {
    NA_real_
  }
  
  # Compute the weighted score for the current subtype
  # Compute the weighted score
  if (!is.na(log_pval) && !is.na(normalized_page) && !is.na(gm_score) && !is.na(gm_var)) {
    weighted_score <- log(log_pval * normalized_page) + log(gm_score / (gm_var + 0.000001))
  } else if (is.na(log_pval)) {
    weighted_score <- log(gm_score / (gm_var + 0.000001))
  } else {
    weighted_score <- NA_real_
  }
}
 
  
  # Apply the function rowwise to calculate weighted scores for each subtype and store in a new column
  metadata_filtered <- metadata_filtered %>%
  rowwise() %>%
  mutate(
    weighted_score_EMR = calculate_weighted_score_per_subtype(cur_data(), "EMRPAGE", "EMR_pvalPAGE", "EMR.score_GM1", "EMR.score_sd_GM1"),
    weighted_score_CRIS_A = calculate_weighted_score_per_subtype(cur_data(), "CRIS_APAGE", "CRIS_A_pvalPAGE", "CRIS_A.score_GM1", "CRIS_A.score_sd_GM1"),
    weighted_score_CRIS_B = calculate_weighted_score_per_subtype(cur_data(), "CRIS_BPAGE", "CRIS_B_pvalPAGE", "CRIS_B.score_GM1", "CRIS_B.score_sd_GM1"),
    weighted_score_CRIS_C = calculate_weighted_score_per_subtype(cur_data(), "CRIS_CPAGE", "CRIS_C_pvalPAGE", "CRIS_C.score_GM1", "CRIS_C.score_sd_GM1"),
    weighted_score_CRIS_D = calculate_weighted_score_per_subtype(cur_data(), "CRIS_DPAGE", "CRIS_D_pvalPAGE", "CRIS_D.score_GM1", "CRIS_D.score_sd_GM1"),
    weighted_score_CRIS_E = calculate_weighted_score_per_subtype(cur_data(), "CRIS_EPAGE", "CRIS_E_pvalPAGE", "CRIS_E.score_GM1", "CRIS_E.score_sd_GM1"),
    weighted_score_mrCAF_Cords = calculate_weighted_score_per_subtype(cur_data(), "mrCAF_CordsPAGE", "mrCAF_Cords_pvalPAGE", "mrCAF_Cords.score_GM1", "mrCAF_Cords_sd_GM1"),
    weighted_score_mrCAF = calculate_weighted_score_per_subtype(cur_data(), "mrCAFPAGE", "mrCAF_pvalPAGE", "mrCAF.score_GM1", "mrCAF.score_sd_GM1"),
    weighted_score_mrCAF_sign = calculate_weighted_score_per_subtype(cur_data(), "mrCAF_signPAGE", "mrCAF_sign_pvalPAGE", "mrCAF_sign.score_GM1", "mrCAF_sign.score_sd_GM1")
    )

  # Function to compute Delta for each spot
  compute_delta <- function(row) {
  sign_names<-names(sign)
  weighted_scores <- as.numeric(row[paste0("weighted_score_", sign_names)])
  #weighted_scores <- as.numeric(row[c("weighted_score_EMR", "weighted_score_MF2",
  #                                    "weighted_score_MF4", "weighted_score_cCAF", 
  #                                    "weighted_score_mrCAF", "weighted_score_myCAF", 
  #                                    "weighted_score_vCAF")])
  
  # Rank the weighted scores
  ranked_scores <- sort(weighted_scores, decreasing = TRUE)
  
  # Calculate the delta (difference between the top two scores)
  delta <- ranked_scores[1] - ranked_scores[2]
  
  return(delta)
}

  # Apply the delta calculation function to each row
  metadata_filtered <- metadata_filtered %>%
  rowwise() %>%
  mutate(delta = compute_delta(cur_data()))

  # Plot the distribution of delta values
  library(ggplot2)
  pdf_filename <- paste0("DeltaValues_", obj_name, ".pdf")
  
  hist_data <- hist(sort(metadata_filtered$delta), breaks = 50, plot = FALSE)
  max_freq_index <- which.max(hist_data$counts)
  inflection_point <- 0.3 #change it
  p <- ggplot(metadata_filtered, aes(x = delta)) +
      geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black", alpha = 0.7) +
      ggtitle("Distribution of Delta Values") +
      xlab("Delta (Difference between Top Two Scores)") +
      ylab("Frequency") +
      theme_minimal() +
      scale_x_continuous(breaks = seq(0, max(metadata_filtered$delta, na.rm = TRUE), by = 0.1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.title.x = element_text(margin = margin(t = 10)),
            axis.title.y = element_text(margin = margin(r = 10))) +
      geom_vline(xintercept = inflection_point, linetype = "dashed", color = "red", size = 1.2)
    
    ggsave(pdf_filename, plot = p, height = 5, width = 12)
  
  
  # Function to calculate Delta and assign a unique classification for each spot
  classify_spot <- function(row, inflection_point) {
    sign_names<-names(sign)
    weighted_scores <- as.numeric(row[paste0("weighted_score_", sign_names)])
    # Extract weighted scores for all subtypes
    #weighted_scores <- as.numeric(row[c("weighted_score_EMR", "weighted_score_MF2",
    #                                     "weighted_score_MF4", "weighted_score_cCAF", 
    #                                     "weighted_score_mrCAF", "weighted_score_myCAF", 
    #                                     "weighted_score_vCAF")])
    
    # Get the names of the subtypes for easy reference
    #subtype_names <- c("EMR", "MF2", "MF4", "cCAF", "mrCAF", "myCAF","vCAF")
    subtype_names<-names(sign)
    
    # Rank the weighted scores
    ranked_scores <- sort(weighted_scores, decreasing = TRUE)
    
    # Calculate the delta (difference between the top two scores)
    delta <- ranked_scores[1] - ranked_scores[2]
    
    # If delta is greater than the threshold, assign the spot to the top scoring subtype
    if (delta > inflection_point) {
      assigned_subtype <- subtype_names[which.max(weighted_scores)]
    } else {
      assigned_subtype <- "negative"  # Classify as negative if delta is too small
    }
    
    return(assigned_subtype)
  }
  
  # Apply the classification function to each row (spot)
  metadata_filtered <- metadata_filtered %>%
    rowwise() %>%
    mutate(assigned_subtype = classify_spot(cur_data(), inflection_point))
  seurat_obj@meta.data$assigned_subtype <- factor(metadata_filtered$assigned_subtype)
  
  # Return the updated Seurat object and the summary table as a list
  return(list(seurat_obj = seurat_obj, subtype_summary = table(metadata_filtered$assigned_subtype)))
}

seurat_objects <- list(
  V573_1 = V573_1, V573_2 = V573_2, V371_1 = V371_1, V371_2 = V371_2, 
  V838_1 = V838_1, V838_2 = V838_2, V763_1 = V763_1, V763_2 = V763_2, 
  V688_1 = V688_1, V688_2 = V688_2, V015_1 = V015_1, V015_2 = V015_2,
  V797_1 = V797_1, V797_2 = V797_2, GSM7058756_C1 = GSM7058756_C1,
  GSM7058757_C2 = GSM7058757_C2, GSM7058758_C3 = GSM7058758_C3, GSM7058759_C4 = GSM7058759_C4
)


processed_results<-list()
for (obj_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[obj_name]]
  # Run the process on each Seurat object
  processed_results[[obj_name]] <- process_seurat_object(seurat_obj, obj_name)
}

# Combine subtype summaries from all datasets into a single data frame for visualization
combine_subtype_summaries <- function(processed_results) {
  
  # Initialize an empty data frame
  combined_data <- data.frame()
  
  # Loop over each result in processed_results
  for (obj_name in names(processed_results)) {
    
    # Extract the subtype summary table
    subtype_summary <- processed_results[[obj_name]]$subtype_summary
    
    # Convert the summary table to a data frame with two columns: 'Subtype' and 'Count'
    subtype_df <- as.data.frame(subtype_summary)
    colnames(subtype_df) <- c("Subtype", "Count")
    
    # Add a column for the dataset name
    subtype_df$Dataset <- obj_name
    
    # Combine with the main data frame
    combined_data <- rbind(combined_data, subtype_df)
  }
  
  return(combined_data)
}

# Combine the subtype summaries
combined_subtype_data <- combine_subtype_summaries(processed_results)

# Visualize the combined data
p <- ggplot(combined_subtype_data, aes(x = Subtype, y = Count, fill = Subtype)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Dataset, scales = "free_y") +  # Create separate plots for each dataset
  theme_minimal() +
  ggtitle("CAF Subtype and EMR Summary Across Datasets") +
  xlab("Subtype") +
  ylab("Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),   # Rotate x-axis labels for readability
        legend.position = "none")  # Remove legend if not needed
ggsave('UniqueAssignemnet_AcrossDatasets.pdf',p,width = 13,height = 10)



##After a first assignment of the spots, I developed a methodology to assess the fraction of cell type per k, expressed as a Z score whose fraction sums up to 1. Note that we compute first a null distribution. This methodology becomes interesting when you want to look at the neighbour for a particular celltype (here EMR)

library(Seurat)
library(dbscan)
library(ggplot2)
library(parallel)

# Function to generate randomized neighbors by k using parallelization
generate_randomized_neighbors_by_k_parallel <- function(seurat_obj, max_k = 40, num_randomizations = 1000, n_cores = 18) {
  metadata <- seurat_obj@meta.data
  CAF_subtypes <- metadata$assigned_subtype
  coords <- seurat_obj@images$slice1$centroids@coords
  rownames(coords) <- rownames(metadata)
  
  # Use mclapply to parallelize over each k value
  results <- mclapply(1:max_k, function(k) {
    knn_spatial <- dbscan::kNN(as.matrix(coords), k = k)
    
    knn_spatial.norm <- data.frame(
      from = rep(1:nrow(knn_spatial$id), k),
      to = as.vector(knn_spatial$id),
      distance = as.vector(knn_spatial$dist)
    )
    
    knn_spatial.norm$from <- rownames(coords)[knn_spatial.norm$from]
    knn_spatial.norm$to <- rownames(coords)[knn_spatial.norm$to]
    
    # Generate randomization results for this k in parallel
    randomization_results <- replicate(num_randomizations, {
      shuffled_subtypes <- sample(CAF_subtypes)
      names(shuffled_subtypes) <- rownames(metadata)
      
      random_neighbors <- lapply(rownames(metadata), function(cell) {
        neighbors <- knn_spatial.norm[knn_spatial.norm$from == cell, "to"]
        neighbor_subtypes <- shuffled_subtypes[neighbors]
        table(factor(neighbor_subtypes, levels = unique(CAF_subtypes)))
      })
      
      random_NN <- Reduce('+', random_neighbors)
      return(random_NN)
    }, simplify = FALSE)
    
    return(randomization_results)
  }, mc.cores = n_cores)
  
  return(results)
}

# Function to calculate normalized z-scores by k using parallelization
calculate_normalized_z_scores_by_k_parallel <- function(real_NN_list, random_NN_list, n_cores = 18) {
  # Use mclapply to parallelize across different k values
  normalized_z_scores_by_k <- mclapply(1:length(real_NN_list), function(k) {
    real_NN <- real_NN_list[[k]]
    random_NN <- random_NN_list[[k]]
    
    # Calculate mean and SD for the randomized neighborhoods
    random_mean <- Reduce('+', random_NN) / length(random_NN)
    random_sd <- sqrt(Reduce('+', lapply(random_NN, function(x) (x - random_mean)^2)) / length(random_NN))
    
    # Calculate z-scores for real vs randomized
    z_scores <- (real_NN - random_mean) / random_sd
    
    # Normalize the z-scores so that they sum to 1 for each k
    z_scores_sum <- sum(z_scores, na.rm = TRUE)  # Sum z-scores, ignoring NAs
    normalized_z_scores <- z_scores / z_scores_sum
    
    # Return normalized z-scores in a data frame for this k
    return(data.frame(
      CAF_subtype = names(z_scores),
      normalized_z_score = as.numeric(normalized_z_scores),
      k = k
    ))
  }, mc.cores = n_cores)
  
  # Combine all k-values into a single data frame for plotting
  return(do.call(rbind, normalized_z_scores_by_k))
}

# Compute real EMR neighbors by k
compute_real_emr_neighbors_by_k <- function(seurat_obj, max_k = 40) {
  coords <- seurat_obj@images$slice1$centroids@coords
  if (is.null(coords) || nrow(coords) == 0) {
    stop("Coordinates are not available or are empty in the Seurat object.")
  }
  metadata <- seurat_obj@meta.data
  rownames(coords) <- rownames(metadata)
  if (!"assigned_subtype" %in% colnames(metadata)) {
    stop("assigned_subtype column not found in metadata.")
  }
  
  emr_cells <- rownames(metadata[metadata$assigned_subtype == "EMR", ])
  
  if (length(emr_cells) == 0) {
    stop("No EMR cells found in the Seurat object.")
  }
  
  results <- list()
  
  for (k in 1:max_k) {
    knn_spatial <- dbscan::kNN(as.matrix(coords), k = k)
    
    knn_spatial.norm <- data.frame(
      from = rep(1:nrow(knn_spatial$id), k),
      to = as.vector(knn_spatial$id),
      distance = as.vector(knn_spatial$dist)
    )
    
    knn_spatial.norm$from <- rownames(coords)[knn_spatial.norm$from]
    knn_spatial.norm$to <- rownames(coords)[knn_spatial.norm$to]
    
    if (nrow(knn_spatial.norm) == 0) {
      warning(paste("No neighbors found for k =", k))
      next
    }
    
    real_neighbors <- lapply(emr_cells, function(cell) {
      neighbors <- knn_spatial.norm[knn_spatial.norm$from == cell, "to"]
      neighbor_subtypes <- metadata[neighbors, "assigned_subtype"]
      return(table(factor(neighbor_subtypes)))
    })
    
    if (length(real_neighbors) == 0 || all(sapply(real_neighbors, length) == 0)) {
      warning(paste("No neighbors for EMR cells at k =", k))
      results[[k]] <- integer(0)
    } else {
      all_subtypes <- unique(unlist(lapply(real_neighbors, names)))
      
      standardized_neighbors <- lapply(real_neighbors, function(tbl) {
        standardized_tbl <- setNames(rep(0, length(all_subtypes)), all_subtypes)
        standardized_tbl[names(tbl)] <- tbl
        return(standardized_tbl)
      })
      
      real_NN <- Reduce('+', standardized_neighbors)
      results[[k]] <- real_NN
    }
  }
  
  return(results)
}



GSM7058758_C3<-processed_results$GSM7058758_C3$seurat_obj

# Example usage
emr_real_neighbors_by_k <- compute_real_emr_neighbors_by_k(GSM7058758_C3, max_k = 40)
emr_random_neighbors_by_k <- generate_randomized_neighbors_by_k_parallel(GSM7058758_C3, max_k = 40, num_randomizations = 500, n_cores = 18)
emr_normalized_z_scores <- calculate_normalized_z_scores_by_k_parallel(emr_real_neighbors_by_k, emr_random_neighbors_by_k, n_cores = 18)

pdf('trial_C3_emr.pdf')
ggplot(emr_normalized_z_scores, aes(x = k, y = normalized_z_score, color = CAF_subtype)) +
  geom_line()+
  theme_minimal() +
  labs(
    title = "Trend of Normalized Z-scores for Neighborhood Enrichment (EMR-Centric)",
    x = "k (Number of Neighbors)",
    y = "Normalized Z-score (0 to 1)"
  ) +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = rainbow(length(unique(emr_normalized_z_scores$CAF_subtype))))
dev.off()


