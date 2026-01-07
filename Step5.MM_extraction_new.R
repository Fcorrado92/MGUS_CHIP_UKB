library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)
library(VennDiagram)
library(grid)

# ##################################################################
# Load cancer registry data extracted from UK Biobank
# ##################################################################
Cancer_registry_data <- read_csv("/mnt/project/extract_fields_ttyd/Cancer_registry_data.csv")
mgus_patients<-read_csv("~/MGUS_screening_output.csv")
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
# Define coding for Multiple Myeloma (MM)
# MM: ICD-10 = C90.0, Histology = 9732
# ------------------------------------------------------------------
mm_icd10 <- "C90"
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
# Extract participant IDs for  MM
# ------------------------------------------------------------------
mm_ids <- Cancer_registry_data$`Participant ID`[has_mm]
mm_data_cr<-Cancer_registry_data%>%filter(`Participant ID`%in%mm_ids)


# ============================================================
# Multiple Myeloma: find the first occurrence of ICD10 C90.0 and its date
# ============================================================

mm_pattern <- "\\bC90\\b"  # Multiple myeloma ICD10

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



# -------------------------------------------------------------------------
#MM histology extract
# -------------------------------------------------------------------------
mm_histology   <- c("9731","9732", "9733", "9734")

# 1) Pivot histology columns to long format
histo_long <- mm_data_cr %>%
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
date_long <- mm_data_cr %>%
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

# -------------------------------------------------------------------------
#Rbind ICD10 and HISTO codes from Cancer registry
# -------------------------------------------------------------------------
mm_first_histo$codes<-"histology"
mm_first$codes<-"ICD10"
mm_first<-mm_first%>%select(c(1,3,4,5))
mm_first_histo<-mm_first_histo%>%select(c(1,3,4,5))

colnames(mm_first)<-c("ID", "Diagnosis","Date_MM", "codes")    
colnames(mm_first_histo)<-c("ID", "Diagnosis","Date_MM", "codes")    
final_cancer_registry<-rbind(mm_first_histo, mm_first)

#remove duplicates
#is always the same date when the patient has histology and ICD10-->NO
#select earliest occurrence
t<-final_cancer_registry%>%group_by(ID)%>%summarise(n=n_distinct(codes))%>%filter(n>1)%>%pull(ID)
sub<-final_cancer_registry%>%filter(ID%in%t)
final_cancer_registry <- final_cancer_registry %>%
  group_by(ID) %>%
  arrange(Date_MM) %>%          # ordina per data MM crescente
  slice(1) %>%                  # tieni solo la prima (earliest)
  ungroup()


final_cancer_registry<-final_cancer_registry%>%select(1,3)
# -------------------------------------------------------------------------
#Look at the overla between these and Primary Care Data
#select MM cases from GP
# -------------------------------------------------------------------------
data<-read_csv("/mnt/project/extract_fields_ttyd/GP_clinical.csv")

data%>%summarise(n=n_distinct(`Participant ID`))
nrow(data)
# 118154716

#each row is an encounter

#MM codes CTV3 and R2
mm_ctv3<-c("B63..","B630.","B6300","B63z.","Xa0SI","Xa0SL","Xa0SN","Xa36a","Xa36b","Xa36c","Xa9AA",
           "XaBB3","XaBLx","XaELI","XE20N")
mm_readv2<-c("BBn2.","BBr30" ,"BBr3z", "BBr3.","BBn0.",
             "BBnz.",
             "BBn1.",
             "BBn3.",
             "BBn..",
             "B6303",
             "B630.",
             "B63..",
             "B936.",
             "B631.",
             "B6301",
             "B936.")

#select rows with MGUS codes
mm_data_v3<-data%>%dplyr::filter(`CTV3 (Read v3) code`%in%mm_ctv3)
mm_data_v2<-data%>%dplyr::filter(`Read v2 code`%in%mm_readv2)
mm_data<-rbind(mm_data_v3, mm_data_v2)

#how many patients and how many patients have multiple encounters?
mm_data%>%summarise(n=n_distinct(`Participant ID`))
# 322

mm_data_filt<-mm_data%>%select(`Participant ID`,`Date clinical code was entered`)
colnames(mm_data_filt)<-c("ID", "Date_MM")
mm_data_filt$registry<-"GP"

mm_data_filt_first_date <- mm_data_filt %>%
  group_by(ID) %>%
  arrange(Date_MM) %>%          # sort by mm date (earliest first)
  slice(1) %>%                  # keep only the first row per participant
  ungroup()



#VEnn diagram
ids_gp <- mm_data_filt_first_date$ID

ids_cancer <- final_cancer_registry$ID


venn.plot <- draw.pairwise.venn(
  area1 = length(ids_gp),
  area2 = length(ids_cancer),
  cross.area = length(intersect(ids_gp, ids_cancer)),
  category = c("GP", "CR"),
  fill = c("gold", "lightgrey"),
  alpha = 0.7,
  cat.cex = 3,
  cex = 3,
  scaled = TRUE
)

grid.newpage()
grid.draw(venn.plot)

gp_only<-setdiff(ids_gp,ids_cancer )

gp_only<-mm_data%>%filter(`Participant ID`%in%gp_only)
ggsave(venn.plot, filename="Venn_MM_registries.pdf", width=10, height=8)


#merge the two registries

final_cancer_registry$registry<-"cancer"

merged<-rbind(final_cancer_registry, mm_data_filt_first_date)
merged_filt <- merged %>%
  group_by(ID) %>%
  arrange(Date_MM) %>%          # sort by mm date (earliest first)
  slice(1) %>%                  # keep only the first row per participant
  ungroup()

write_csv(merged_filt, "MM_screening_output.csv")
