# Load libraries and data
library(tidyverse)
library(readxl)
library(haven)
library(readr)
library(writexl)
library(scales)
library(rstatix)
library(ggpubr)
library(cmprsk)
library(tidycmprsk)
library(ggsurvfit)

chip_data_path <- "~/Partners HealthCare Dropbox/Francesco Corrado/Sabine Francesco/1_CHIP_PLCO/Analyzed Data/PLCO_CHIP_data.xlsx"
chip_data <- read_excel(chip_data_path, sheet = 1)

demographics_path<- "~/Partners HealthCare Dropbox/Francesco Corrado/Sabine Francesco/1_CHIP_PLCO/Demographics and Clinical Data/Analysis File/e_2019_1020_analysis_d112524.sas7bdat"
demographics <- read_sas(demographics_path)

#during this whole analysis correct for the year:
#all PLCO variables are computed from randomization timepoint
#re-compute the variables from the dna sample timepoint
## dna_sample_drawdays Days from Randomization to DNA Collection
## dna_sample_studyyr Study Year of DNA Collection
## dna_sample_yr Calendar Year of DNA Collection

demographics_new <- demographics %>%
  mutate(j_hema_exitdays_fromdna = j_hema_exitdays - dna_sample_drawdays)%>%
  mutate(age_atdna = round(age + dna_sample_drawdays/365.5), 0)

# Specify the mutation column names
chip_genes <- c(
  "DNMT3A",
  "TET2",
  "GNB1",
  "SF3B1",
  "U2AF1",
  "ZNF318",
  "TP53",
  "ASXL1",
  "PPM1D",
  "JAK2",
  "ZBTB33",
  "KRAS",
  "IDH1",
  "SRSF2"
)

#save demo_data
plco_chip_demo_data <- demographics_new%>%
  select(plco_id, casestat, 
    j_hema_cancer, j_hema_cancer_diagdays, j_hema_cancer_first, j_hema_candxdays,
         j_hema_behavior, j_hema_grade, j_hema_seer, j_hema_seercat, j_hema_type,
         d_dth_hema, f_dth_hema, d_cancersite, 
         d_seer_death, f_seer_death, is_dead_with_cod, is_dead, 
         mortality_exitage, mortality_exitstat,
         dth_days, mortality_exitdays, 
         j_hema_exitstat, j_hema_exitdays, j_hema_exitdays_fromdna, 
         j_ph_any_trial, 
         dna_sample_studyyr, dna_sample_drawdays, dna_masked_id, dna_sample_yr, dna_sample_time,
         age, age_atdna, race7, sex, fh_cancer, hema_fh,
         hyperten_f, hearta_f, stroke_f, arthrit_f, osteopor_f, 
         diabetes_f, polyps_f, divertic_f, gallblad_f, bronchit_f)

plco_chip_demo_data_path <- "~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PLCO_CHIP_demo_data.xlsx"
write_xlsx(plco_chip_demo_data, plco_chip_demo_data_path)

#summary of CHIP results per casestat
#summary of the data
casestat_summary <- plco_chip_demo_data %>%
  group_by(casestat) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(casestat = factor(casestat, levels = c(0, 1, 2),
                           labels = c("No MG\n No Cancer", "MG Case\n Progressed to Myeloma", "MG Case\n No Myeloma"))) %>%
  # Adding a total row
  bind_rows(summarise(., casestat = "Total", Count = sum(Count)))

print(casestat_summary)

# # A tibble: 4 × 2
# casestat                          Count
# <chr>                             <int>
# 1 No MG/n No Cancer                 315
# 2 MG Case/n Progressed to Myeloma   158
# 3 MG Case/n No Myeloma              460
# 4 Total                               933

#Plot casestat distribution

casestat_plot <- casestat_summary %>%
  filter(casestat != "Total") %>%
  mutate(
    Percent = Count / sum(Count),
    casestat = factor(
      casestat,
      levels = c(
        "No MG\n No Cancer",
        "MG Case\n No Myeloma",
        "MG Case\n Progressed to Myeloma"
      )
    )
  )

p <- ggplot(casestat_plot, aes(x = casestat, y = Count, fill = casestat)) +
  geom_col(width = 0.7, color = "black") +
  geom_text(
    aes(label = paste0(Count, "\n(", percent(Percent, accuracy = 0.1), ")")),
    vjust = -0.3,
    size = 6
  ) +
  scale_fill_manual(
    values = c(
      "#9ecae1",  # No MG / No Cancer
      "#31a354",  # MG no MM
      "#de2d26"  # Progressed to MM
      
    )
  ) +
  labs(
    x = NULL,
    y = "Number of individuals",
    title = "Distribution of case status in PLCO cohort(n=933)"
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "italic"),
    axis.text.x = element_text(size = 12)
  ) +
  ylim(0, max(casestat_plot$Count) * 1.15)

p
ggsave(p, filename="~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PLCO_casestat_distribution.pdf", width = 8, height = 6)


#Age and Follow-up distribution
table(plco_chip_demo_data$casestat, plco_chip_demo_data$j_hema_exitstat)
# 1="Confirmed Cancer" 3="Last Participant Contact Prior to Unconfirmed Report" 4="Last Participant Contact" 5="Death" 6="Date Lost/Refused" 7="Registry Search Completeness Date" 8="Cancer Free at Cutoff" 9="Post-2009 Death, Exit At 12/31/09"
# 1   3   4   5   6   7   8   9
# 0   0   0 102  96   3 109   3   2
# 1 158   0   0   0   0   0   0   0
# 2   7   6 149 142   0 152   2   2

#rename casestat
plco_chip_demo_data <- plco_chip_demo_data %>%
  mutate(
    casestat = factor(
      casestat,
      levels = c(0, 2, 1),
      labels = c(
        "No MG\nNo Cancer",
        "MG Case\nNo Myeloma",
        "MG Case\nProgressed to Myeloma"
      )
    )
  )

comparisons <- list(
  c("No MG\nNo Cancer", "MG Case\nProgressed to Myeloma"),
  c("No MG\nNo Cancer", "MG Case\nNo Myeloma"),
  c("MG Case\nProgressed to Myeloma", "MG Case\nNo Myeloma")
)


p <- ggplot(plco_chip_demo_data, aes(x = casestat, y = j_hema_exitdays, fill = casestat)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.4) +
  scale_fill_manual(values = c(
    "#9ecae1",  # No MG / No Cancer
    "#31a354",  # MG no MM
    "#de2d26"
  )) +
  labs(
    x = NULL,
    y = "Follow-up Days",
    title = ""
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )+stat_compare_means(method = "t.test", comparisons = comparisons, hide.ns = TRUE,size=6)

n_df <- plco_chip_demo_data %>%
  group_by(casestat) %>%
  summarise(n = n(), .groups = "drop")

p<-p +
  geom_text(
    data = n_df,
    aes(x = casestat, y = -Inf, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = -0.6,
    size = 6
  ) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
  )

ggsave(p, filename="~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PLCO_FU_days.pdf", width = 8, height = 6)

#AGE
p <- ggplot(plco_chip_demo_data, aes(x = casestat, y = age_atdna, fill = casestat)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.4) +
  scale_fill_manual(values = c(
    "#9ecae1",  # No MG / No Cancer
    "#de2d26",  # Progressed to MM
    "#31a354"   # MG no MM
  )) +
  labs(
    x = NULL,
    y = "Age",
    title = ""
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )+stat_compare_means(method = "t.test", comparisons = comparisons, hide.ns = TRUE,size=6)

n_df <- plco_chip_demo_data %>%
  group_by(casestat) %>%
  summarise(n = n(), .groups = "drop")

p<-p +
  geom_text(
    data = n_df,
    aes(x = casestat, y = -Inf, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = -0.6,
    size = 6
  ) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 30, l = 10)
  )
ggsave(p, filename="~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PLCO_Age_days.pdf", width = 8, height = 6)


#merge with CHIP calls from SK
chip_plco<-read_csv("~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/CHIP_PLCO_overview_SK.csv")




chip_plco_wide <- chip_plco %>%
  mutate(
    INFO2 = str_replace_all(INFO, "'", '"'),
    INFO2 = na_if(INFO2, "{}")
  ) %>%
  mutate(
    parsed = map(INFO2, ~{
      if (is.na(.x)) {
        # placeholder keeps the row through unnest/pivot
        return(tibble(gene = ".none", vaf = NA_real_))
      }
      x <- jsonlite::fromJSON(.x)
      tibble(gene = names(x), vaf = as.numeric(unname(x)))
    })
  ) %>%
  select(-c(INFO2)) %>%
  unnest(parsed) %>%
  pivot_wider(names_from = gene, values_from = vaf) %>%
  select(-any_of(".none"))   # remove placeholder column

# sanity check
nrow(chip_plco_wide)
#load linker file
linker_plco<-read_xlsx("~/Partners HealthCare Dropbox/Francesco Corrado/Sabine Francesco/1_CHIP_PLCO/Linker file PLCO Vanderbilt sequencing.xlsx")
linker_plco$Current_Label<-linker_plco$`Current Label`
chip_plco_wide<-chip_plco_wide%>%left_join(linker_plco[c("sample_id", "Current_Label")], by="sample_id")

#merge with demographics
plco_chip_demo_data$Current_Label<-plco_chip_demo_data$dna_masked_id
plco_chip_demo_data<-left_join(plco_chip_demo_data,chip_plco_wide, by="Current_Label")

#Demographics CHIP vs noCHIP
#diagnosis and progression group
plco_chip_demo_data <- plco_chip_demo_data %>%
  mutate(
    diagnosis_group = case_when(
      grepl("No MG", casestat, ignore.case = TRUE) ~ "No MG",
      TRUE ~ "MG"
    )
  )
#filter only MG
plco_mg<-plco_chip_demo_data%>%filter(diagnosis_group=="MG")

plco_mg <- plco_mg %>%
  mutate(
    MM_Progressor = case_when(
      grepl("Progressed to Myeloma", casestat, ignore.case = TRUE) ~ "Progressor",
      grepl("No Myeloma", casestat, ignore.case = TRUE) ~ "Non-Progressor",
      TRUE ~ NA_character_
    )
  )

#compute PFS_event
plco_mg<-plco_mg%>%mutate(PFS_event=ifelse(MM_Progressor=="Progressor", "1", "0"))
#compute PFS_months
plco_mg<-plco_mg%>%mutate(PFS_mos=j_hema_exitdays_fromdna/30)
#competing risk death PFS_months
plco_mg<-plco_mg%>%mutate(PFS_event=ifelse(j_hema_exitstat=="5","2", PFS_event))

plco_mg_filt<-plco_mg%>%filter(!is.na(CHIP))%>%select(c("diagnosis_group","age_atdna","CHIP", "sample_id", "MM_Progressor", "PFS_event", "PFS_mos", "race7", chip_genes))

plco_mg_filt <- plco_mg_filt %>%
  mutate(
    CHIP = factor(CHIP, levels = c("N", "Y"), labels = c("CHIP neg", "CHIP pos")),
    race7 = factor(race7)
  )

n_df <- plco_mg_filt %>%
  filter(!is.na(age_atdna), !is.na(CHIP)) %>%
  group_by(CHIP) %>%
  summarise(n = n(), .groups = "drop")


p_age <- ggplot(plco_mg_filt, aes(x = CHIP, y = age_atdna, fill = CHIP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.2) +
  scale_fill_manual(values = c("#bdbdbd", "#3182bd")) +
  labs(
    x = NULL,
    y = "Age at DNA draw",
    title = "Age distribution by CHIP status (n=609)"
  ) +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",    plot.title = element_text(hjust = 0.5, face = "italic") ) +
  geom_pwc(method = "t.test", label = "p.format", label.size=5)

p_age <- p_age +
  geom_text(
    data = n_df,
    aes(x = CHIP, y = -Inf, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = -0.6,
    size = 5
  ) +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 35, l = 10)
  )

p_age

ggsave(p_age, filename="~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/Age_CHIP_PLCO.pdf", width = 8, height = 6)



plco_mg_filt <- plco_mg_filt %>%
  mutate(
    race_bin = case_when(
      race7 == 1 ~ "White",
      race7 %in% c(2, 3, 4, 5) ~ "Non-White",
      TRUE ~ NA_character_
    ),
    race_bin = factor(race_bin, levels = c("White", "Non-White"))
  )

race_df <- plco_mg_filt %>%
  filter(!is.na(race_bin), !is.na(CHIP)) %>%
  count(CHIP, race_bin) %>%
  group_by(CHIP) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

p_race_bin <- ggplot(race_df, aes(x = CHIP, y = percent, fill = race_bin)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(
    aes(label = percent(percent, accuracy = 1)),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 5
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c(
      "White" = "#f0f0f0",
      "Non-White" = "#636363"
    )
  ) +
  labs(
    x = NULL,
    y = "Percentage",
    fill = "Race",
    title = "Race distribution by CHIP status (n=609)"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic")
  )

p_race_bin
ggsave(p_race_bin, filename="~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/Race_CHIP_PLCO.pdf", width = 8, height = 6)


# -------------------------------------------------------------------------
#competing risk analysis
# -------------------------------------------------------------------------
time<-plco_mg_filt$PFS_mos
event<-plco_mg_filt$PFS_event
chip<-plco_mg_filt$CHIP




# -------------------------------------------------------------------------
#cumulative risk of progression CHIP
# -------------------------------------------------------------------------

library(dplyr)
library(tidycmprsk)
library(survival)

plco_mg_filt_cr <- plco_mg_filt %>%
  mutate(
    PFS_mos = as.numeric(PFS_mos),
    PFS_event = as.integer(PFS_event),   
    status = factor(
      case_when(
        PFS_event == 0 ~ "censored",
        PFS_event == 1 ~ "progression",
        PFS_event == 2 ~ "death",
        TRUE ~ NA_character_
      ),
      levels = c("censored", "progression", "death")
    ),
    CHIP = factor(CHIP)  # già factor va bene
  ) %>%
  filter(!is.na(PFS_mos), !is.na(status), !is.na(CHIP), PFS_mos >= 0)

#REMOVE 1 for PFS mos <0
ci <- tidycmprsk::cuminc(Surv(PFS_mos, status) ~ CHIP, data = plco_mg_filt_cr)

p_prog  <- ci$cmprsk$Tests[1, "pv"]  # cause 1 = progression
p_death <- ci$cmprsk$Tests[2, "pv"]  # cause 2 = death

plco_inc_plot_death <- ggsurvfit::ggcuminc(ci, outcome = "death") +
  add_risktable() +
  labs(
    x = "Months",
    y = "Cumulative incidence Death",
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
  ) +
  coord_cartesian(clip = "off") 

plco_inc_plot_death

pdf("~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PLCO/competing_risk_death_plco.pdf",
    width = 8, height = 6)

print(plco_inc_plot_death)

dev.off()



plco_inc_plot_progression <- ggsurvfit::ggcuminc(ci, outcome = "progression") +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Gray p = ", scales::pvalue(p_prog, accuracy = 0.001)),
    hjust = 1.05, vjust = -0.8,
    size = 5
  ) +
  coord_cartesian(clip = "off") +  
  labs(x="Months", y="Cumulative incidence Progression", color="CHIP status") +
  theme_classic(base_size = 16) +
  add_risktable() +
  scale_color_manual(values = c("#808080", "#1f78b4"))
ggsave(
  filename = "~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PLCO/competing_risk_progression_plco.pdf",
  plot = plco_inc_plot_progression,
  width = 8,
  height = 6)



# -------------------------------------------------------------------------
#Do the same for each mutation( mut+ vs CHIP-)
# -------------------------------------------------------------------------
plco_mg_filt_cr <- plco_mg_filt_cr %>%
  mutate(
    across(all_of(chip_genes), ~ as.numeric(replace_na(.x, 0)))
  ) %>%
  mutate(
    across(all_of(chip_genes), ~ as.integer(.x > 0), .names = "{.col}_bin")
  )


gene_counts <- plco_mg_filt_cr %>% select(-c(race_bin))%>%
  summarise(across(ends_with("_bin"), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "gene", values_to = "n_pos") %>%
  mutate(gene = sub("_bin$", "", gene)) %>%
  arrange(desc(n_pos))

gene_counts

p<-ggplot(gene_counts, aes(x = reorder(gene, n_pos), y = n_pos)) +
  geom_col(color = "black", width = 0.75) +
  geom_text(aes(label = n_pos), vjust = -0.3, size = 5) +
  theme_classic(base_size = 18) +
  labs(
    x = NULL,
    y = "Number of individuals wth gene mutation",
    title = "CHIP gene prevalence (absolute counts)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(  filename = "~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/number_of_mutations.pdf",
         plot = p,
         width = 8,
         height = 6)

  
outdir<-"~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PLCO/"

plco_mg_filt_cr <- plco_mg_filt_cr %>%
  mutate(
    across(
      all_of(chip_genes),
      ~ replace_na(as.numeric(.x), 0)
    )
  )
plot_gene_vs_chipneg <- function(dat, gene, outdir) {
  
  dat2 <- dat %>%
    mutate(
      gene_pos = .data[[gene]] > 0,
      group = case_when(
        gene_pos ~ paste0(gene, "+"),
        CHIP %in% c("N", "CHIP neg", "CHIP-", "CHIP−") ~ "CHIP-",
        TRUE ~ NA_character_   # escludi gli altri CHIP+ senza questo gene
      ),
      group = factor(group, levels = c("CHIP-", paste0(gene, "+")))
    ) %>%
    filter(!is.na(group))
  
  # serve almeno 2 gruppi
  if (nrow(dat2) == 0 || length(unique(dat2$group)) < 2) return(NULL)
  
  ci <- tidycmprsk::cuminc(Surv(PFS_mos, status) ~ group, data = dat2)
  

  p <- ggsurvfit::ggcuminc(ci, outcome = "progression") +
    add_risktable() +
    labs(
      x = "Months",
      y = "Cumulative incidence Progression",
      color = NULL,
      title = paste0("Progression: ", gene, "+ vs CHIP-")
    ) +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) 
  ggsave(
    filename = paste0(outdir,"competing_risk_progression_", gene, "_vs_CHIPneg.pdf"),
    plot = p,
    width = 8, height = 6  )
  
  return(p)
}  


plots_by_gene <- lapply(chip_genes, function(g) {
  message("Processing ", g, " ...")
  plot_gene_vs_chipneg(plco_mg_filt_cr, g, outdir)
})

names(plots_by_gene) <- chip_genes





# -------------------------------------------------------------------------
#do it for mean VAF
# -------------------------------------------------------------------------

plco_mg_filt_cr <- plco_mg_filt_cr %>%
  rowwise() %>%
  mutate(
    mean_vaf = {
      v <- c_across(all_of(chip_genes))
      v_pos <- v[v > 0]
      if (length(v_pos) == 0) NA_real_ else mean(v_pos)
    }
  ) %>%
  ungroup()

plco_mg_filt_cr <- plco_mg_filt_cr %>%
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

plco_mg_filt_cr_chip_only<-plco_mg_filt_cr%>%filter(CHIP=="CHIP pos")
ci <- tidycmprsk::cuminc(Surv(PFS_mos, status) ~ mean_vaf_group, data = plco_mg_filt_cr_chip_only)

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

ggsave(  filename = "~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PLCO/CR_VAF_PLCO.pdf",
         plot = chip_vaf_CR,
         width = 8,
         height = 6)