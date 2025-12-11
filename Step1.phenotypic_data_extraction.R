library(readr)

# ------------------------------------------------------------------
# Load cancer registry data extracted from UK Biobank
# ------------------------------------------------------------------
Cancer_registry_data <- read_csv("/mnt/project/extract_fields_ttyd/Cancer_registry_data.csv")

# Inspect number of rows and column names
nrow(Cancer_registry_data)
colnames(Cancer_registry_data)

# ------------------------------------------------------------------
# Identify all columns that contain ICD-10 cancer codes
# (multiple 'Instance' columns per participant)
# ------------------------------------------------------------------
icd10_cols <- grep("^Type of cancer: ICD10", names(Cancer_registry_data), value = TRUE)

# Identify all columns that contain cancer histology codes
histology_cols <- grep("^Histology", names(Cancer_registry_data), value = TRUE)

# ------------------------------------------------------------------
# Define coding for MGUS
# MGUS: ICD-10 = D47.2, Histology = 9765
# ------------------------------------------------------------------
mgus_icd10 <- "D47.2"
mgus_histology <- "9765"

# Flag participants with MGUS ICD-10 code in *any* ICD-10 column
has_mgus_icd10 <- apply(
  Cancer_registry_data[icd10_cols], 1,
  function(x) any(grepl(mgus_icd10, x))
)

# Flag participants with MGUS histology code in *any* histology column
has_mgus_histology <- apply(
  Cancer_registry_data[histology_cols], 1,
  function(x) any(grepl(mgus_histology, x))
)

# Participants with MGUS by either ICD-10 OR histology
has_mgus <- has_mgus_icd10 | has_mgus_histology

# Count MGUS cases
sum(has_mgus)

# ------------------------------------------------------------------
# Define coding for Multiple Myeloma (MM)
# MM: ICD-10 = C90.0, Histology = 9732
# ------------------------------------------------------------------
mm_icd10 <- "C90.0"
mm_histology <- "9732"

# Flag participants with MM ICD-10 code
has_mm_icd10 <- apply(
  Cancer_registry_data[icd10_cols], 1,
  function(x) any(grepl(mm_icd10, x))
)

# Flag participants with MM histology code
has_mm_histology <- apply(
  Cancer_registry_data[histology_cols], 1,
  function(x) any(grepl(mm_histology, x))
)

# Participants with MM by either ICD-10 OR histology
has_mm <- has_mm_icd10 | has_mm_histology

# Count MM cases
sum(has_mm)

# ------------------------------------------------------------------
# Extract participant IDs for MGUS, MM, and their overlap
# ------------------------------------------------------------------
mgus_ids <- Cancer_registry_data$eid[has_mgus]
mm_ids <- Cancer_registry_data$eid[has_mm]
overlap_ids <- intersect(mgus_ids, mm_ids)

mgus_ids
mm_ids
overlap_ids