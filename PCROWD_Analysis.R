library(readxl)
library(dplyr)
library(ggplot2)
library(kableExtra)
library(readr)
library(stringr)
library(tidyr)
library(cmprsk)
library(tidycmprsk)
library(ggsurvfit)
#load df
input_dir<-"~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/"
data<-read_xlsx("~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/Francesco CHIP samples - with annotations and diagnoses UPDATED 10.31.2025_FINAL.xlsx")
outdir<-"~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PCROWD/"

#remove those not yet annotated
data_filt<-data%>%filter(!remove_november=="1")
#number of samples
length(unique(data_filt$Prep_ID))
#455
#number of pts
length(unique(data_filt$pi_dfci_mrn))
#308
t<-data_filt%>%group_by(pi_dfci_mrn)%>%summarise(n=n_distinct(Prep_ID))
t2<-prop.table(table(t$n))
t2<-as.data.frame(t2)

#plot distribution n samples x ID
colnames(t2) <- c("n_samples_per_patient", "proportion")

# Plot stacked (una barra unica con proporzioni + numeri)
p1<-ggplot(t2, aes(x = 1, y = proportion, fill = as.factor(n_samples_per_patient))) +
  geom_col(color = "black", width=0.3) +
  geom_text(aes(label = paste0(round(proportion * 100, 1), "%")),
            position = position_stack(vjust = 0.5), size = 5, color = "black") +
  scale_fill_brewer(palette = "RdYlBu", name = "Samples per patient") +
  labs(
    x = NULL, y = "Proportion of patients",
    title = "Number of samples\n per subject(n=308)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )
ggsave(p1, filename=paste0(input_dir, "samplexpatient.pdf"), width=4, height = 4)

#Diagnosis at first sampling date
data_first <- data_filt %>%
  # Make sure the sample date is in Date format
  mutate(SAMPLEDATE = as.Date(SAMPLEDATE)) %>%
  # Sort samples by patient ID and sample date (oldest first)
  arrange(pi_dfci_mrn, SAMPLEDATE) %>%
  # Keep only the earliest (first) sample for each patient
  group_by(pi_dfci_mrn) %>%
  slice_head(n = 1) %>%
  ungroup()
#distribution intial diagnosis
table(data_first$diagnosis_at_date_of_sample)

#check if treated SMM have pi_startdate< SampleDATE
data_first %>%
  filter(diagnosis_at_date_of_sample == "Treated_SMM") %>%
  mutate(
    pi_startdate = as.Date(pi_startdate),
    SAMPLEDATE = as.Date(SAMPLEDATE),
    treated_before_sample = pi_startdate <= SAMPLEDATE
  ) %>%
  summarise(
    n_total = n(),
    n_treated_before_sample = sum(treated_before_sample, na.rm = TRUE),
    proportion = mean(treated_before_sample, na.rm = TRUE)
  )

#created DX2 where Treated and untreated are collaped, Unknown, WM, AMyloidosis are Other
# Collapse diagnosis categories into broader groups
data_first <- data_first%>%mutate(DX2=ifelse(diagnosis_at_date_of_sample %in% 
                                               c("Unknown Disease Stage", "WM", "Amyloidosis", "Other"), "Other",diagnosis_at_date_of_sample)
)
# check the new grouping
first_diagnosis<-as.data.frame(table(data_first$DX2, useNA = "ifany"))
colnames(first_diagnosis)<-c("Diagnosis", "n")
# MGUS    MM Other   SMM 
# 75    32    20   331 

first_diagnosis <- first_diagnosis %>%
  mutate(Percentage = round(100 * n / sum(n), 1))

# save as PDF
first_diagnosis$Diagnosis<-gsub("_", " ", first_diagnosis$Diagnosis)
first_diagnosis$Diagnosis<-factor(first_diagnosis$Diagnosis, levels=c("MGUS", "Untreated SMM",
                                                                      "Treated SMM", "NDMM", 
                                                                      "Treated MM", "Other"))
# reorder rows based on factor levels
first_diagnosis <- first_diagnosis %>%
  arrange(factor(Diagnosis, levels = levels(Diagnosis)))

# nicely formatted table
kbl(first_diagnosis, booktabs = TRUE, align = "lcc",
    col.names = c("Diagnosis", "Number of Patients", "% of Total"),
    caption = "Diagnosis at First Sampling Date")

#Filter Untreated SMM and MGUS then I want a swimmer plot were first data point is sample date 
#second data point shaped differently is second sample date
#third an other dtaa points are third , fourth etc sample date
#then I want to highlight pi_startdate
#and date of last follow-up
#I want a line per patient
#and I want the color of line to reflect the diagnosis associated with that sampledate
untreated_smm_first_diagnosis<-data_first%>%filter(diagnosis_at_date_of_sample%in%c("MGUS","Untreated_SMM"))%>%pull(pi_dfci_mrn)
untreat_prec <- data_filt %>%
  filter(pi_dfci_mrn %in% untreated_smm_first_diagnosis) %>%
  mutate(SAMPLEDATE = as.Date(SAMPLEDATE)) %>%
  arrange(pi_dfci_mrn, SAMPLEDATE) %>%
  group_by(pi_dfci_mrn) %>%
  mutate(
    sample_number = row_number(),
    sample_label = paste0("Sample", sample_number)
  ) %>%
  ungroup()

#how many are treated?
treatedSMM<-untreat_prec%>%filter(pi_treated=="1")%>%pull(pi_dfci_mrn)
treatedSMM<-unique(treatedSMM)

#of the ones that were not intially treated filter those who will receive treatment and calculate TFI
treated_SMM_wide<-data_first%>%filter(pi_dfci_mrn%in%treatedSMM & diagnosis_at_date_of_sample=="Untreated_SMM")%>%mutate(follow_up_until_treatment=(as.Date(pi_startdate)-as.Date(SAMPLEDATE)))

treated_SMM_wide <- treated_SMM_wide %>%
  mutate(
    followup_categories = case_when(
      follow_up_until_treatment < 30 ~ "<1 month",
      follow_up_until_treatment < 365 ~ "<1 year",
      follow_up_until_treatment < 1825 ~ "<5 years",
      follow_up_until_treatment >= 1825 ~ ">5 years",
      TRUE ~ NA_character_
    )
  )

# Count patients per follow-up category
followup_summary <- treated_SMM_wide %>%
  group_by(followup_categories) %>%
  summarise(n = n()) %>%
  mutate(followup_categories = factor(
    followup_categories,
    levels = c("<1 month", "<1 year", "<5 years", ">5 years")
  ))

# Bar plot with counts annotated
p2<-ggplot(followup_summary, aes(x = followup_categories, y = n, fill = followup_categories)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.5, size = 6) +
  scale_fill_brewer(palette = "Blues", name = "Follow-up") +
  labs(
    x = "Treatment Free Interval",
    y = "Number of patients",
    title = "Treatment Free Interval in uSMM"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
ggsave(filename=paste0(input_dir, "TFI_uSMM.pdf"), width=6, height=6, p2)


# ##now take only those uSMM that will not receive therapy and look at how many progressed
# untreated_SMM_wide<-data_first%>%filter(!pi_dfci_mrn%in%treatedSMM & diagnosis_at_date_of_sample=="Untreated_SMM")
# table(untreated_SMM_wide$MM_Progressor)
# untreated_SMM_wide<-untreated_SMM_wide%>%mutate(FU=as.Date(pi_last_appt_date)-as.Date(SAMPLEDATE))
# 
# untreated_SMM_wide <- untreated_SMM_wide %>%
#   mutate(
#     followup_categories = case_when(
#       FU < 30 ~ "<1 month",
#       FU < 365 ~ "<1 year",
#       FU < 1825 ~ "<5 years",
#       FU >= 1825 ~ ">5 years",
#       TRUE ~ NA_character_
#     )
#   )
# 
# # Count patients per follow-up category
# followup_summary <- untreated_SMM_wide %>%
#   group_by(followup_categories) %>%
#   summarise(n = n()) %>%
#   mutate(followup_categories = factor(
#     followup_categories,
#     levels = c("<1 month", "<1 year", "<5 years", ">5 years")
#   ))
# 
# # Bar plot with counts annotated
# p3<-ggplot(followup_summary, aes(x = followup_categories, y = n, fill = followup_categories)) +
#   geom_col(color = "black", width = 0.7) +
#   geom_text(aes(label = n), vjust = -0.5, size = 6) +
#   scale_fill_brewer(palette = "Blues", name = "Follow-up") +
#   labs(
#     x = "Follow-up",
#     y = "Number of patients",
#     title = "Follow-up in uSMM"
#   ) +
#   theme_classic(base_size = 16) +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5)
#   )
# ggsave(filename=paste0(input_dir, "FU_SMM_not_treated.pdf"), width=6, height=6, p3)
# 
# #MGUS
# mgus_first_diagnosis<-data_first%>%filter(diagnosis_at_date_of_sample=="MGUS")
# table(mgus_first_diagnosis$MM_Progressor)


# -------------------------------------------------------------------------
#Now load sungjae df 
# -------------------------------------------------------------------------
df<-read_csv("~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/CHIP_PCROWD_overview_SK.csv")
# Clean up the mutation string
df_clean <- df %>%
  mutate(
    # Remove curly braces, quotes, and extra spaces
    mutations = str_remove_all(INFO, "[\\{\\}'\"]"),
    mutations = str_replace_all(mutations, " ", ""),
    # Split each list of mutations by comma (creates a list column)
    mutations = str_split(mutations, ",")
  )

# Expand (unnest) to one gene per row
df_long <- df_clean %>%
  unnest(mutations) %>%
  mutate(
    # Extract the gene name (everything before the colon)
    gene = str_extract(mutations, "^[A-Za-z0-9]+"),
    # Extract the numeric value (VAF) after the colon
    VAF = as.numeric(str_extract(mutations, "[0-9.]+$"))
  ) %>%
  select(-c(mutations, INFO))

# Transform into wide format (one column per gene)
df_wide <- df_long %>%
  pivot_wider(
    names_from = gene,  # column names = gene names
    values_from = VAF   # values = VAF
  ) %>%
  mutate(
    # Calculate the minimum VAF per patient across all mutations
    min_VAF = apply(select(., where(is.numeric)), 1, min, na.rm = TRUE)
  )



# -------------------------------------------------------------------------
#merge jackie and sungjae df
# -------------------------------------------------------------------------
df_wide$Sample_ID_VB<-df_wide$sample_id
df_wide$Sample_ID_VB<-gsub("U","",df_wide$Sample_ID_VB)

#take only first sample per patient
untreat_prec_first<-untreat_prec%>%filter(sample_number=="1")
untreat_prec_first<-untreat_prec_first%>%left_join(df_wide, by="Sample_ID_VB")

#PFS event = 1 if Progression, =0 if last follow up without progression =2 if death 3= if treatment
#If progression occurs after treatment -> 3
#check within progressors how many received treatment before progression

progressors<-untreat_prec_first%>%filter(MM_Progressor=="Progressor")%>%mutate(delta=pi_mm_diagnosis_date- pi_startdate)

#if treatment precede progression use treatment date for PFS

sub <- untreat_prec_first %>%
  mutate(
    PFS_event = case_when(
      
      # 3 = progressed but MM diagnosis happened AFTER startdate (prevalent)
      MM_Progressor == "Progressor" & pi_treated=="1" &
        pi_mm_diagnosis_date > pi_startdate ~ 3,
      
      MM_Progressor == "Progressor" & pi_treated=="1" &
        pi_mm_diagnosis_date < pi_startdate ~ 1,
      
      MM_Progressor == "Progressor" & pi_treated=="2" ~ 1,
      
      # 0 = non-progressor, untreated, alive
      MM_Progressor == "Non-Progressor" &
        pi_treated == "2" & pi_surv_status == "1" ~ 0,
      
      # 2 = death (competing event)
      MM_Progressor == "Non-Progressor" &
        pi_treated == "2" & pi_surv_status == "2" ~ 2,
      
      MM_Progressor == "Non-Progressor" &
        pi_treated == "1" ~ 3,
      
      # everything else -> NA (così te ne accorgi)
      TRUE ~ NA_real_
    )
  )

#remove NA due to treatmen==1 but no date 
sub_filt<-sub%>%filter(!is.na(PFS_event))


sub_filt <- sub_filt %>%
  mutate(
    PFS_date = case_when(
      # 0 = censored → last follow-up
      PFS_event == 0 ~ pi_last_appt_date,
      
      # 1 = progression → MM diagnosis date
      PFS_event == 1 ~ pi_mm_diagnosis_date,
      
      # 2 = death → date of death
      PFS_event == 2 ~ pi_death_date,
      
      # 3 = prevalent MM / treated at baseline → treatment start
      PFS_event == 3 ~ pi_startdate,
      
      TRUE ~ as.Date(NA)
    )
  )

sub_filt<-sub_filt%>%mutate(PFS_mos=as.numeric((as.Date(PFS_date)-as.Date(pi_date_initial_dx))/30))

sub_filt2 <- sub_filt %>%
  mutate(
    # 1) forza numerico (gestisce "0"/"1"/etc)
    PFS_event_num = suppressWarnings(as.integer(as.character(PFS_event))),
    
    # 2) tieni solo 0/1/2/3
    PFS_event_num = ifelse(PFS_event_num %in% c(0,1,2,3), PFS_event_num, NA_integer_),
    
    # 3) factor con ordine: censored FIRST
    status = factor(
      PFS_event_num,
      levels = c(0, 1, 2, 3),
      labels = c("censored", "progression", "death", "treatment")
    )
  ) %>%
  filter(!is.na(PFS_mos), !is.na(status), !is.na(CHIP))

# sanity check
table(sub_filt2$status, useNA = "ifany")
#Competing risk analysis


ci <- tidycmprsk::cuminc(Surv(PFS_mos, status) ~ CHIP, data = sub_filt2)



p_prog  <- ci$cmprsk$Tests[1, "pv"]  # cause 1 = progression
p_death <- ci$cmprsk$Tests[2, "pv"]  # cause 2 = death
p_treatment<-ci$cmprsk$Tests[3, "pv"]

pcrowd_inc_plot_progression <- ggsurvfit::ggcuminc(ci, outcome = "progression") +
  add_risktable() +
  labs(
    x = "Months",
    y = "Cumulative incidence Progression",
    color = "CHIP"
  ) +
  theme_classic(base_size = 18) +
  scale_color_manual(values = c("#808080", "#1f78b4"))+
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Gray p = ", scales::pvalue(p_prog, accuracy = 0.001)),
    hjust = 1.05, vjust = -0.8,
    size = 5
  ) +
  coord_cartesian(clip = "off") 

pcrowd_inc_plot_progression

pdf("~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PCROWD/competing_risk_prog_pcrowd.pdf",
    width = 8, height = 6)

print(pcrowd_inc_plot_progression)

dev.off()


pcrowd_inc_plot_death <- ggsurvfit::ggcuminc(ci, outcome = "death") +
  add_risktable() +
  labs(
    x = "Months",
    y = "Cumulative incidence death",
    color = "CHIP"
  ) +
  theme_classic(base_size = 18) +
  scale_color_manual(values = c("#808080", "#1f78b4"))+
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Gray p = ", scales::pvalue(p_death, accuracy = 0.001)),
    hjust = 1.05, vjust = -0.8,
    size = 5
  )

pcrowd_inc_plot_death

pdf("~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PCROWD/competing_risk_death_pcrowd.pdf",
    width = 8, height = 6)

print(pcrowd_inc_plot_death)

dev.off()


pcrowd_inc_plot_treatment <- ggsurvfit::ggcuminc(ci, outcome = "treatment") +
  add_risktable() +
  labs(
    x = "Months",
    y = "Cumulative incidence treatment",
    color = "CHIP"
  ) +
  theme_classic(base_size = 18) +
  scale_color_manual(values = c("#808080", "#1f78b4"))+
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Gray p = ", scales::pvalue(p_treatment, accuracy = 0.001)),
    hjust = 1.05, vjust = -0.8,
    size = 5
  )

pcrowd_inc_plot_treatment

pdf("~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PCROWD/competing_risk_treatment_pcrowd.pdf",
    width = 8, height = 6)

print(pcrowd_inc_plot_treatment)

dev.off()



# -------------------------------------------------------------------------
#Mean VAF
# -------------------------------------------------------------------------
colnames(sub_filt2)

chip_genes <- c(
  "PPM1D", "DNMT3A", "JAK2", "TP53", "SRSF2", "TET2", "U2AF1",
  "GNB1", "KRAS", "SF3B1", "ZNF318", "IDH1", "ASXL1", "ZBTB33"
)

sub_filt2 <- sub_filt2 %>%
  rowwise() %>%
  mutate(
    mean_vaf = {
      v <- c_across(all_of(chip_genes))
      v[is.na(v)] <- 0          # 🔑 trasforma NA in 0
      v_pos <- v[v > 0]
      if (length(v_pos) == 0) 0 else mean(v_pos)
    }
  ) %>%
  ungroup()

sub_filt2 <- sub_filt2 %>%
  mutate(
    mean_vaf_group = case_when(
      CHIP %in% c("CHIP neg", "N") | is.na(mean_vaf) | mean_vaf == 0 ~ "CHIP neg",
      mean_vaf < 0.1 ~ "<10%",
      mean_vaf > 0.1 ~ ">10%",
      TRUE ~ NA_character_
    ),
    mean_vaf_group = factor(
      mean_vaf_group,
      levels = c("CHIP neg", "<10%", ">10%")
    )
  )

sub_filt2_chip_only<-sub_filt2%>%filter(CHIP=="Y")
ci <- tidycmprsk::cuminc(Surv(PFS_mos, status) ~ mean_vaf_group, data = sub_filt2_chip_only)

p_prog  <- ci$cmprsk$Tests[1, "pv"]  # cause 1 = progression
p_death <- ci$cmprsk$Tests[2, "pv"]  # cause 2 = death

chip_vaf_CR<-ggsurvfit::ggcuminc(ci, outcome = "progression") +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Gray p = ", scales::pvalue(p_prog, accuracy = 0.001)),
    hjust = 1.05, vjust = -0.8,
    size = 5
  ) +
  add_risktable()+
  labs(x="Months", y="Cumulative incidence Progression", color="Mean VAF") +
  theme_classic(base_size = 16) 

ggsave(  filename = "~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PCROWD/CR_VAF_PCROWD.pdf",
         plot = chip_vaf_CR,
         width = 8,
         height = 6)

write_csv(sub_filt2,"~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PCROWD/Competitive_risk_pcrowd.csv")

