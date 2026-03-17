#####
#Step1. demultiplex Jannsen flowcell
#####
#load dependencies
library(vcfR)
library(data.table)
library(dplyr)
library(readr)
library(stringr)
library(readxl)
library(pheatmap)
library(tidyverse)



#date
date <- format(Sys.Date(), "%m.%d.%Y")
date<-"01.10.2026"

#read wgs vcf
gvcf_wgs<-read.vcfR("/mnt/disks/cellranger/Janssen_delivery_souporcell/CHIP_DARARVD_010925_CHIP_DARARVD_demux_sorted.vcf")
gvcf_wgs_df<-cbind.data.frame(gvcf_wgs@fix, gvcf_wgs@gt)

#souporcell folder that has all the pool directories
#select first the ones with 4 samples 
input_dir<-"~/Immune_project/Accelerator_souporcell_v8/"
# Ottieni i nomi delle directory che corrispondono al pattern
samples <- list.files(input_dir, recursive = FALSE)
remove <- samples[grep("Accelerator", samples)]

samples<-setdiff(samples,remove)
#remove other type of files(vcfs, xlsx)
samples<-samples[-c(19:length(samples))]

#load the sequencing sheet to understand how many samples per pool 
df<-read_xlsx( "/mnt/disks/cellranger/Janssen_delivery_souporcell/DaraRVD_JJ_demux_list.xlsx")
n_samples_x_pool<-df%>%group_by(Batch)%>%summarise(n=n_distinct(Sample_ID))

#loop through pools
for(k in 1:length(samples))
{
  # read souporcell vcf for pool (if 3 samples in pool, vcf will have 3 clusters based on the what was specified while running souporcell)
  pool1_vcf<-read.vcfR(paste0(input_dir, "/",samples[k], "/cluster_genotypes.vcf")) #variant count : ~125000
  pool1_df<-cbind.data.frame(pool1_vcf@fix, pool1_vcf@gt)
  
  #create new df with SNPs that are common between wgs gvcf and souporcell genotype vcf 
  new_df<-inner_join(gvcf_wgs_df, pool1_df, by=c("CHROM","POS"))
  
  # [1] "CHROM"              "POS"                "ID.x"               "REF.x"              "ALT.x"              "QUAL.x"             "FILTER.x"           "INFO.x"             "FORMAT.x"          
  # [10] "Ultra_GS6_1_2_PB"   "Ultra_GS6_1_3_PB"   "Ultra_GS6_1_5_PB"   "Ultra_GS6_2_4_PB"   "Ultra_GS6_3_3_PB"   "Ultra_GS6_3_6_PB"   "Ultra_GS7_1_2_PB"   "Ultra_GS7_1_5_PB"   "Ultra_GS7_2_2_PB"  
  # [19] "Ultra_GS7_2_4_PB"   "Ultra_GS7_3_5_PB"   "Ultra_GS7_3_6_PB"   "Ultra_GS7_3_9_PB"   "Ultra_GS8_2_3_PB"   "Ultra_GS8_2_4_PB"   "Ultra_GS8_3_2_PB"   "Ultra_GS8_3_4_PB"   "Ultra_GS8_4_3_PB"  
  # [28] "Ultra_GS8_4_4_PB"   "Ultra_RT2_31_PB"    "Ultra_RT2_32_PB"    "Ultra_RT2_33_PB"    "Ultra_RT3_1_1_PB"   "Ultra_RT3_4_3_PB"   "Ultra_RT3_4_5_PB"   "Ultra_RT4_1_2_PB"   "Ultra_RT4_1_6_PB"  
  # [37] "Ultra_RT4_2_3_PB"   "Ultra_RT4_2_4_PB"   "Ultra_RT4_3_3_PB"   "Ultra_RT4_3_4_PB"   "Ultra_RT5_1_2_PBMC" "Ultra_RT5_1_4_PBMC" "Ultra_RT5_2_2_PBMC" "Ultra_RT5_2_5_PBMC" "Ultra_RT5_3_3_PBMC"
  # [46] "Ultra_RT5_3_5_PBMC" "Ultra_RT5_4_3_PBMC" "Ultra_RT5_4_4_PBMC" "Ultra_sci_14_16_PB" "Ultra_sci_19_13_PB" "Ultra_sci_20_2_PB"  "Ultra_sci_21_10_PB" "Ultra_sci_21_4_PB"  "pM5338g"           
  # [55] "ID.y"               "REF.y"              "ALT.y"              "QUAL.y"             "FILTER.y"           "INFO.y"             "FORMAT.y"           "0"                  "1"                 
  # [64] "2"           
  
  #keep only selected columns
  new_df<-new_df[, c(1,2,4,5,10:ncol(gvcf_wgs_df), (ncol(gvcf_wgs_df)+8): ncol(new_df))] # keep the last three columns  from souporcell clusters (here: 62-64) and the columns with genotyping from wgs, remove columns with not useful information
  
  new_df2<-new_df
  nrow(new_df2) 
  
  ###only keep first three strings for genotype columns
  new_df2[,-c(1:4)]<-lapply(new_df2[, -c(1:4)], 
                            function(x) substr(x,1,3))
  
  
  #only keep rows where genotype across all column is either 0/0, 1/0, 0/1, 1/1
  new_df3<-new_df2
  new_df3[, -c(1:4)]<-lapply(new_df3[,-c(1:4)],function(x) (x=="0/1"| x=="1/1" | x=="0/0")) #apply to columns that have genotype data, gives true or false for condition
  idx <- which(apply(new_df3[,5:ncol(new_df3)], 1, all)) # get all that are true
  new_df4<-new_df2[idx,]
  
  #calculate degree of match between wgs vcf and souporcell vcf (pairwise comparisons between wgs and souporcell clusters)
  new_df5<-new_df4
  
  GenotypeConcordance<-data.frame(matrix(ncol=2, nrow=0))
  
  for (i in 5:(ncol(gvcf_wgs_df)-5)) #wgs genotypes
  {
    for(j in (ncol(gvcf_wgs_df)-4):ncol(new_df4)) #scrNA souporcell genotypes
    {
      wgs_sample_name<-colnames(new_df5)[i]
      pool_sample_cluster_name<-colnames(new_df5)[j]
      
      comb_name<-paste0(wgs_sample_name, "_soupcell_cluster",pool_sample_cluster_name)
      
      new_df5$new_column<-new_df5[,i]==new_df5[,j]
      
      out<-sum(new_df5$new_column, na.rm=TRUE)/nrow(new_df5)
      
      names(new_df5)[names(new_df5) == "new_column"]<-comb_name
      
      GenotypeConcordance<-rbind(GenotypeConcordance, list(comb_name, out))
      
      
    }
  }
  
  colnames(GenotypeConcordance)<-c("Pair", "Prcnt_concordance")
  
  saveRDS(GenotypeConcordance, paste0(input_dir, "/",samples[k],"/cluster_genotypes", date, ".rds"))
}

# Initialize an empty list to store top4 data frames
list <- list()

for (i in seq_along(samples)) {
  # Read the RDS file for the current sample
  data_i <- read_rds(paste0(input_dir, "/", samples[i], "/cluster_genotypes",date,".rds"))
  
  # Add the sample name as a new column called 'Pool'
  data_i$Pool <- gsub( "_GEX_5", "",samples[i])
  
  # Add the top4 data frame to the list
  list[[i]] <- data_i
}

# Combine all the data frames into one
final <- bind_rows(list)

final <- final %>%
  mutate(
    Genotype = str_extract(Pair, ".*(?=_soupcell)"),
    soupcell_cluster = str_extract(Pair, "\\d+$")  
  )%>%select(c(-1))

final <- final %>%
  mutate(
    Pool_cluster = paste0(Pool, "_cluster", soupcell_cluster),
    Prcnt_concordance = as.numeric(Prcnt_concordance)
  )


# ============================================================
# 1) Build wide matrix:
#    rows = Pool_cluster (e.g. B41_cluster0)
#    columns = Genotype (WGS samples)
# ============================================================

mat_df <- final %>%
  pivot_wider(
    names_from  = Genotype,
    values_from = Prcnt_concordance,
    values_fill = NA
  ) %>%
  as.data.frame()

# Set row names using Pool_cluster and remove non-numeric columns
rownames(mat_df) <- mat_df$Pool_cluster
mat_df$Pool_cluster <- NULL
mat_df$Pool <- NULL
mat_df$soupcell_cluster <- NULL

# Ensure all remaining columns are numeric
mat_df <- mat_df %>% mutate(across(everything(), as.numeric))

# Convert to matrix and express values as percentages
mat <- as.matrix(mat_df) * 100


# ============================================================
# 2) Order ROWS by Pool (B41, B42, ...) and then by cluster index
# ============================================================

# Extract pool and cluster information from row names
row_pool    <- sub("_cluster.*", "", rownames(mat))
row_cluster <- as.numeric(sub(".*_cluster", "", rownames(mat)))
pool_num    <- as.numeric(sub("^B", "", row_pool))

# Define row order: first by pool number, then by cluster number
row_order <- order(pool_num, row_cluster)
mat <- mat[row_order, , drop = FALSE]


# ============================================================
# 3) Convert matrix to long format
#    (needed to compute column ordering logic)
# ============================================================

long <- as.data.frame(mat) %>%
  mutate(Pool_cluster = rownames(mat)) %>%
  pivot_longer(
    -Pool_cluster,
    names_to  = "Genotype",
    values_to = "Concordance"
  ) %>%
  mutate(
    Pool = sub("_cluster.*", "", Pool_cluster)
  )


# ============================================================
# 4) Rule for number of genotypes per pool:
#    - B41 and B42: top 6
#    - all other pools: top 4
# ============================================================

topN_for_pool <- function(p) {
  if (p %in% c("B41", "B42")) 6 else 4
}


# ============================================================
# 5) Build COLUMN ORDER (core step for block-diagonal structure)
#    - Within each pool, rank genotypes by concordance
#    - Keep top N genotypes per pool
#    - Concatenate pools in numeric order (B41, B42, ...)
#    - Remove duplicates, keeping first occurrence
# ============================================================

geno_order <- long %>%
  group_by(Pool) %>%
  arrange(desc(Concordance), .by_group = TRUE) %>%
  mutate(rk = row_number()) %>%
  filter(rk <= topN_for_pool(Pool[1])) %>%
  ungroup() %>%
  arrange(as.numeric(sub("^B", "", Pool))) %>%
  pull(Genotype) %>%
  unique()

# Keep only genotypes that are present in the matrix
geno_order <- geno_order[geno_order %in% colnames(mat)]

# Reorder matrix columns accordingly
mat <- mat[, geno_order, drop = FALSE]


# ============================================================
# 6) (Optional but recommended) Define gaps between pools
#    for better visual separation in the heatmap
# ============================================================

# Row gaps: separate pools along rows
row_gaps <- cumsum(table(row_pool))

# Column gaps: separate pools along columns
col_pool <- sapply(colnames(mat), function(g) {
  long$Pool[match(g, long$Genotype)]
})
col_gaps <- cumsum(table(col_pool))


# ============================================================
# 7) Final heatmap: block-diagonal genotype concordance
# ============================================================

pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = row_gaps,
  gaps_col = col_gaps,
  na_col = "grey90",
  border_color = NA,
  main = "WGS vs souporcell genotype concordance (%)"
)

# ============================================================
# 8) Generate Final CSV with Best Matches per Pool
#    Using topN logic: B41/B42 = 6, others = 4
# ============================================================

message("\n=== Generating Final Best Matches CSV ===")

# Convert concordance back to 0-1 scale for consistency
final_for_export <- final %>%
  mutate(Prcnt_concordance = as.numeric(Prcnt_concordance))

# For each pool, select top N best matches based on concordance
best_matches <- final_for_export %>%
  group_by(Pool) %>%
  arrange(desc(Prcnt_concordance), .by_group = TRUE) %>%
  mutate(
    topN = topN_for_pool(Pool[1]),
    rank = row_number()
  ) %>%
  filter(rank <= topN) %>%
  ungroup()

# For each cluster, extract the single best matching genotype
best_match_per_cluster <- best_matches %>%
  group_by(Pool, Pool_cluster) %>%
  slice_max(Prcnt_concordance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Cluster = sub(".*_cluster", "", Pool_cluster)) %>%
  select(Pool, Cluster, Pool_cluster, Genotype, Concordance = Prcnt_concordance, soupcell_cluster) %>%
  arrange(
    as.numeric(sub("^B", "", Pool)),
    as.numeric(Cluster)
  )

message("Total clusters with assignments: ", nrow(best_match_per_cluster))

# Quality Control Checks
message("\n=== Quality Control ===")

# Check for low concordance matches (< 60%)
low_concordance <- best_match_per_cluster %>%
  filter(Concordance < 0.60)

if (nrow(low_concordance) > 0) {
  message("WARNING: ", nrow(low_concordance), " cluster(s) with concordance < 60%:")
  print(low_concordance)
} else {
  message("All clusters have concordance >= 60%")
}

# Check for duplicate genotype assignments within pools
duplicates <- best_match_per_cluster %>%
  group_by(Pool, Genotype) %>%
  filter(n() > 1) %>%
  ungroup()

if (nrow(duplicates) > 0) {
  message("\nWARNING: ", nrow(duplicates)/2, " genotype(s) assigned to multiple clusters:")
  print(duplicates %>% arrange(Pool, Genotype))
} else {
  message("No duplicate genotype assignments within pools")
}

# Summary statistics by pool
summary_stats <- best_match_per_cluster %>%
  group_by(Pool) %>%
  summarise(
    n_clusters = n(),
    n_genotypes = n_distinct(Genotype),
    avg_concordance = mean(Concordance, na.rm = TRUE),
    min_concordance = min(Concordance, na.rm = TRUE),
    max_concordance = max(Concordance, na.rm = TRUE),
    .groups = "drop"
  )

message("\n=== Summary by Pool ===")
print(summary_stats)

# Save final results
output_file <- paste0("~/CHIP_Immune_RNA_seq/final_best_matches_JJ_", date, ".csv")
write_csv(best_match_per_cluster, output_file)
message("\nFinal best matches saved to: ", output_file)

# Save the full CSV to match the demux_CHIP.csv format
write_csv(final_for_export, "~/CHIP_Immune_RNA_seq/demux_CHIP.csv")
message("Updated demux_CHIP.csv with all concordance data")

message("\n=== DONE ===")
message("Review the final_best_matches_JJ_", date, ".csv file for cluster assignments")