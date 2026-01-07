library(renv)
renv::init()
install.packages("readr")
install.packages("survminer")
install.packages("cmprsk")
install.packages("tidycmprsk")

library(readr)
library(dplyr)
library(tidyverse)
library(VennDiagram)
library(grid)


# -------------------------------------------------------------------------
#extract codes from primary care data
# -------------------------------------------------------------------------
data<-read_csv("/mnt/project/extract_fields_ttyd/GP_clinical.csv")

nrow(data)
# 118154716

#each row is an encounter

#MGUS codes CTV3 and R2
mgus_ctv3<-c("C331.","C332.","C332z","Xa0l6","Xa0SJ","Xa36d","Xa36e","Xa36h","Xa36k","Xa36m","XaIt4","XE11b")
mgus_readv2<-c("BBm7.",               "C3322",               "C3310",               "C331.",
               "C332.",               "C332z")

#select rows with MGUS codes
mgus_data_v3<-data%>%dplyr::filter(`CTV3 (Read v3) code`%in%mgus_ctv3)
mgus_data_v2<-data%>%dplyr::filter(`Read v2 code`%in%mgus_readv2)
mgus_data<-rbind(mgus_data_v3, mgus_data_v2)


#how many patients and how many patients have multiple encounters?
mgus_data%>%summarise(n=n_distinct(`Participant ID`))
# 521

encounter_n<-mgus_data%>%group_by(`Participant ID`)%>%summarise(n=n_distinct(`Date clinical code was entered`))
ggplot(encounter_n, aes(x=n))+geom_histogram()
# table(encounter_n$n)
# 
# 1   2   3   4   5   6   7  10 
# 385  87  28   6   9   4   1   1 

#MGUS code is mostly inserted only once


#Select only cols ID and Date and then keep only first date
mgus_data_filt<-mgus_data%>%select(`Participant ID`,`Date clinical code was entered`, `Data provider`)
colnames(mgus_data_filt)<-c("ID", "Date_MGUS","Data_Provider")

mgus_data_filt_first_date <- mgus_data_filt %>%
  group_by(ID) %>%
  arrange(Date_MGUS) %>%          # sort by MGUS date (earliest first)
  slice(1) %>%                  # keep only the first row per participant
  ungroup()


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
# Extract participant IDs for MGUS
# ------------------------------------------------------------------
mgus_ids <- Cancer_registry_data$`Participant ID`[has_mgus]
mgus_data_cr<-Cancer_registry_data%>%filter(`Participant ID`%in%mgus_ids)


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



# -------------------------------------------------------------------------
#bind GP and CR mgus 
# -------------------------------------------------------------------------
colnames(mgus_first)
colnames(mgus_data_filt_first_date)

mgus_data_filt_first_date$registry<-"GP"
mgus_first$registry<-"cancer"


mgus_data_filt_first_date<-mgus_data_filt_first_date%>%select(c(1,2,4))

mgus_first$ID<-mgus_first$`Participant ID`
mgus_first$Date_MGUS<-mgus_first$mgus_date

mgus_first<-mgus_first%>%select(c(5:7))

final_mgus_pool<-rbind(mgus_first, mgus_data_filt_first_date)

#VEnn diagram
ids_gp <- final_mgus_pool %>%
  filter(registry == "GP") %>%
  pull(ID) %>%
  unique()

ids_cancer <- final_mgus_pool %>%
  filter(registry == "cancer") %>%
  pull(ID) %>%
  unique()


venn.plot <- draw.pairwise.venn(
  area1 = length(ids_gp),
  area2 = length(ids_cancer),
  cross.area = length(intersect(ids_gp, ids_cancer)),
  category = c("GP", "CR"),
  fill = c("gold", "lightgrey"),
  alpha = 0.7,
  cat.cex = 1.8,
  cex = 1.84,
  scaled = TRUE
)

grid.newpage()
grid.draw(venn.plot)

ggsave(venn.plot, filename="Venn_MGUS_registries.pdf", width=10, height=8)
#for duplicates take the earliest date

first_date_mgus <- final_mgus_pool %>%
  group_by(ID) %>%
  arrange(Date_MGUS) %>%          # sort by mm date (earliest first)
  slice(1) %>%                  # keep only the first row per participant
  ungroup()


df_plot <- first_date_mgus %>%
  count(registry, name = "n") %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p<-ggplot(df_plot, aes(x="1", y = prop, fill = registry)) +
  geom_col(color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title="Registry Distribution MGUS Diagnosis", y = "Percentage", fill = "Registry") +
  scale_fill_manual(values=c("lightgrey", "gold"))+
  theme_classic() +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), text=element_text(size=18), plot.title = element_text(hjust=0.5))
ggsave(p, filename="Registry_distrib_MGUS.pdf", width=8, height=6)


write_csv(first_date_mgus,"~/MGUS_screening_output.csv")

