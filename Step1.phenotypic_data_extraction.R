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



# -------------------------------------------------------------------------
#Load data from self-reported, primary-care, hospital admission
# -------------------------------------------------------------------------
data<-read_csv("/mnt/project/extract_fields_ttyd/Other_sources_data.csv")
nrow(data)
colnames(data)

# MGUS ICD10 code
mgus_code <- "D47.2"

# TRUE if the string contains D47.2 as a full code
has_mgus <- grepl("\\bD47\\.2\\b", data$`Diagnoses - ICD10`)

# Replace NA with FALSE
has_mgus[is.na(has_mgus)] <- FALSE

# How many MGUS
sum(has_mgus)

# If you have an ID column (e.g. eid), get their IDs
mgus_ids <- data$`Participant ID`[has_mgus]

# Subset full rows
mgus_data <- data[has_mgus, ]


# MM ICD10 code
mm_icd10 <- "C90.0"

# Flag rows where the ICD10 string contains MM (C90.0)
has_mm <- grepl("\\bC90\\.0\\b", data$`Diagnoses - ICD10`)

# Replace NA with FALSE
has_mm[is.na(has_mm)] <- FALSE

# Number of MM cases
sum(has_mm)

# Extract Participant IDs
mm_ids <- data$`Participant ID`[has_mm]

# Extract full MM rows
mm_data <- data[has_mm, ]


mgus_ids <- unique(data$`Participant ID`[has_mgus])
mm_ids <- unique(data$`Participant ID`[has_mm])
overlap_ids <- intersect(mgus_ids, mm_ids)

mgus_ids
mm_ids
overlap_ids
