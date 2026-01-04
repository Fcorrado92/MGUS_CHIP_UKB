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

chip_genes <- c(
  "PPM1D", "DNMT3A", "JAK2", "TP53", "SRSF2", "TET2", "U2AF1",
  "GNB1", "KRAS", "SF3B1", "ZNF318", "IDH1", "ASXL1", "ZBTB33"
)
pcrowd<-read_csv("~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PCROWD/Competitive_risk_pcrowd.csv")

pcrowd<-pcrowd%>%select(c("sample_id", chip_genes, "CHIP","PFS_event",                          "PFS_mos"                           
                                      ,        "status"                            
                          , "mean_vaf"            ,               "mean_vaf_group" ))

pcrowd<-pcrowd%>%mutate(CHIP=ifelse(CHIP=="N", "CHIP neg", "CHIP pos"))
plco<-read_csv("~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/PCROWD_PLCO_Analysis/PLCO/Competitive_risk_plco.csv")

plco<-plco%>%select(c("sample_id", chip_genes, "CHIP","PFS_event",                          "PFS_mos"                           
                                       ,        "status"                            
                          , "mean_vaf"            ,               "mean_vaf_group" ))

merged<-rbind(plco, pcrowd)
merged <- merged %>%
  mutate(
    status = case_when(
      PFS_event == 0 ~ "censored",
      PFS_event == 1 ~ "progression",
      PFS_event == 2 ~ "death",
      PFS_event == 3 ~ "treatment"
    ),
    status = factor(
      status,
      levels = c("censored", "progression", "death", "treatment")
    )
  )
ci <- tidycmprsk::cuminc(Surv(as.numeric(PFS_mos), status) ~ CHIP, data = merged)

ggsurvfit::ggcuminc(ci, outcome = "progression") +
  add_risktable() +
  labs(
    x = "Months",
    y = "Cumulative incidence treatment",
    color = "CHIP"
  ) +
  theme_classic(base_size = 18) +
  scale_color_manual(values = c("#808080", "#1f78b4"))


#HARMONYZE START DATE(SAMPLE VS DIAGNOSIS)

