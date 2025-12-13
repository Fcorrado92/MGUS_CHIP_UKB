library(dplyr)
library(tidyr)
library(stringr)

# ##################################################################
# Load cancer registry data extracted from UK Biobank
# ##################################################################
Cancer_registry_data <- read_csv("/mnt/project/extract_fields_ttyd/Cancer_registry_data.csv")

# Inspect number of rows and column names
nrow(Cancer_registry_data)
#501936

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
mgus_ids <- Cancer_registry_data$`Participant ID`[has_mgus]
mm_ids <- Cancer_registry_data$`Participant ID`[has_mm]


mgus_data_cr<-Cancer_registry_data%>%filter(`Participant ID`%in%mgus_ids)
mm_data_cr<-Cancer_registry_data%>%filter(`Participant ID`%in%mm_ids)


# ============================================================
# MGUS: find the first occurrence of ICD10 D47.2 and its date
# ============================================================

mgus_pattern <- "\\bD47\\.2\\b"  # matches "D47.2" even if followed by description

# 1) Pivot ICD10 columns to long format: one row per participant x instance
icd_long_mgus <- mgus_data_cr %>%
  select(`Participant ID`, matches("^Type of cancer: ICD10 \\| Instance \\d+$")) %>%
  pivot_longer(
    cols = -`Participant ID`,
    names_to = "icd_col",
    values_to = "icd10"
  ) %>%
  mutate(
    instance = as.integer(str_extract(icd_col, "\\d+$"))
  ) %>%
  select(`Participant ID`, instance, icd10) %>%
  filter(!is.na(icd10) & icd10 != "")

# 2) Pivot diagnosis date columns to long format: one row per participant x instance
date_long_mgus <- mgus_data_cr %>%
  select(`Participant ID`, matches("^Date of cancer diagnosis \\| Instance \\d+$")) %>%
  pivot_longer(
    cols = -`Participant ID`,
    names_to = "date_col",
    values_to = "dx_date"
  ) %>%
  mutate(
    instance = as.integer(str_extract(date_col, "\\d+$")),
    dx_date = as.Date(dx_date)  # if parsing fails: as.Date(dx_date, format="%Y-%m-%d")
  ) %>%
  select(`Participant ID`, instance, dx_date)

# 3) Join ICD10 + date by participant and instance, filter MGUS, keep earliest by date
mgus_first <- icd_long_mgus %>%
  left_join(date_long_mgus, by = c("Participant ID", "instance")) %>%
  filter(str_detect(icd10, mgus_pattern)) %>%
  arrange(`Participant ID`, dx_date, instance) %>%
  group_by(`Participant ID`) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    `Participant ID`,
    mgus_instance = instance,
    mgus_icd10    = icd10,
    mgus_date     = dx_date
  )


# ============================================================
# Multiple Myeloma: find the first occurrence of ICD10 C90.0 and its date
# ============================================================

mm_pattern <- "\\bC90\\.0\\b"  # Multiple myeloma ICD10

# 1) Pivot ICD10 columns to long format
icd_long_mm <- mm_data_cr %>%
  select(`Participant ID`, matches("^Type of cancer: ICD10 \\| Instance \\d+$")) %>%
  pivot_longer(
    cols = -`Participant ID`,
    names_to = "icd_col",
    values_to = "icd10"
  ) %>%
  mutate(
    instance = as.integer(str_extract(icd_col, "\\d+$"))
  ) %>%
  select(`Participant ID`, instance, icd10) %>%
  filter(!is.na(icd10) & icd10 != "")

# 2) Pivot diagnosis date columns to long format
date_long_mm <- mm_data_cr %>%
  select(`Participant ID`, matches("^Date of cancer diagnosis \\| Instance \\d+$")) %>%
  pivot_longer(
    cols = -`Participant ID`,
    names_to = "date_col",
    values_to = "dx_date"
  ) %>%
  mutate(
    instance = as.integer(str_extract(date_col, "\\d+$")),
    dx_date  = as.Date(dx_date)
  ) %>%
  select(`Participant ID`, instance, dx_date)

# 3) Join ICD10 + date, filter MM, keep earliest by date
mm_first <- icd_long_mm %>%
  left_join(date_long_mm, by = c("Participant ID", "instance")) %>%
  filter(str_detect(icd10, mm_pattern)) %>%
  arrange(`Participant ID`, dx_date, instance) %>%
  group_by(`Participant ID`) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    `Participant ID`,
    mm_instance = instance,
    mm_icd10    = icd10,
    mm_date     = dx_date
  )

# Combine MGUS + MM and compute time difference
final_cancer_registry_ICD10 <- mgus_first %>%
  left_join(mm_first, by = "Participant ID") %>%
  mutate(delta_days = as.numeric(mm_date - mgus_date))

colnames(final_cancer_registry_ICD10)
final_cancer_registry_ICD10<-final_cancer_registry_ICD10%>%select("Participant ID","mgus_date", "mm_date", "delta")
colnames(final_cancer_registry_ICD10)<-c("ID", "Date_MGUS",  "Date_MM"    ,"delta") 
# -------------------------------------------------------------------------
#MGUS histology extract
# -------------------------------------------------------------------------
mgus_histology <- "9765"
mm_histology   <- "9732"

# 1) Pivot histology columns to long format
histo_long <- mgus_data_cr %>%
  select(`Participant ID`, matches("^Histology of cancer tumour \\| Instance \\d+$")) %>%
  pivot_longer(
    cols = -`Participant ID`,
    names_to = "histo_col",
    values_to = "histology"
  ) %>%
  mutate(
    instance  = as.integer(str_extract(histo_col, "\\d+$")),
    histology = as.character(histology)
  ) %>%
  select(`Participant ID`, instance, histology) %>%
  filter(!is.na(histology) & histology != "")

# 2) Pivot diagnosis date columns (reuse this if you already created date_long_mm)
date_long <- mgus_data_cr %>%
  select(`Participant ID`, matches("^Date of cancer diagnosis \\| Instance \\d+$")) %>%
  pivot_longer(
    cols = -`Participant ID`,
    names_to = "date_col",
    values_to = "dx_date"
  ) %>%
  mutate(
    instance = as.integer(str_extract(date_col, "\\d+$")),
    dx_date  = as.Date(dx_date)
  ) %>%
  select(`Participant ID`, instance, dx_date)

# 3) MGUS histology: join + filter + earliest by date
mgus_first_histo <- histo_long %>%
  left_join(date_long, by = c("Participant ID", "instance")) %>%
  filter(histology == mgus_histology) %>%
  arrange(`Participant ID`, dx_date, instance) %>%
  group_by(`Participant ID`) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    `Participant ID`,
    mgus_histo_instance = instance,
    mgus_histology      = histology,
    mgus_histo_date     = dx_date
  )

mm_first_histo <- histo_long %>%
  left_join(date_long, by = c("Participant ID", "instance")) %>%
  filter(histology == mm_histology) %>%
  arrange(`Participant ID`, dx_date, instance) %>%
  group_by(`Participant ID`) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    `Participant ID`,
    mm_histo_instance = instance,
    mm_histology      = histology,
    mm_histo_date     = dx_date
  )

final_cancer_registry_HISTO <- mgus_first_histo %>%
  left_join(mm_first_histo, by = "Participant ID") %>%
  mutate(delta_days_histo = as.numeric(mm_histo_date - mgus_histo_date))

colnames(final_cancer_registry_HISTO)
final_cancer_registry_HISTO<-final_cancer_registry_HISTO%>%select("Participant ID","mgus_histo_date", "mm_histo_date", "delta_days_histo")
colnames(final_cancer_registry_HISTO)<-c("ID", "Date_MGUS",  "Date_MM"    ,"delta") 


# -------------------------------------------------------------------------
#Rbind ICD10 and HISTO codes from Cancer registry
# -------------------------------------------------------------------------
final_cancer_registry_HISTO$codes<-"histology"
final_cancer_registry_ICD10$codes<-"ICD10"
final_cancer_registry<-rbind(final_cancer_registry_HISTO, final_cancer_registry_ICD10)
#remove duplicates
final_cancer_registry<-final_cancer_registry%>%distinct(ID, .keep_all = TRUE)

write_csv(final_cancer_registry,"~/mgus_data_CR.csv")

final_mgus_data<-rbind(final_cancer_registry, mgus_data_filt_first_date)
final_mgus_data <- final_mgus_data %>%
  mutate(Date_MGUS = as.Date(Date_MGUS)) %>%   
  arrange( ID, Date_MGUS) %>%     
  group_by( ID) %>%
  slice(1) %>%
  ungroup()

write_csv(final_mgus_data,"~/mgus_data_CR_GP.csv")
