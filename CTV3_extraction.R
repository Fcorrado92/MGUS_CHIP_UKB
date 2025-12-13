library(readr)
library(dplyr)

CTV3<-read_csv("~/Partners HealthCare Dropbox/Francesco Corrado/CHIP_Immune_RNAseq/UKBB/CTV3_codes.csv")

mgus_terms <- CTV3 %>%
  dplyr::filter(
    grepl("Monoclonal|Paraproteinaemia", Details, ignore.case = TRUE)
  )

mm_terms <- CTV3 %>%
  dplyr::filter(
    grepl("Multiple Myeloma|Myeloma|Plasmacytoma|Plasma cell", Details, ignore.case = TRUE)
  ) %>%
  dplyr::pull(CTV3) %>%   # estrae direttamente il vettore
  unique()                # rimuove eventuali duplicati