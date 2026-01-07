merged_filt<-read_csv("MM_screening_output.csv")
first_date_mgus<-read_csv("~/MGUS_screening_output.csv")
first_date_mgus<-first_date_mgus%>%rename(registry_mgus="registry")
merged_filt<-merged_filt%>%rename(registry_MM="registry")

#VEnn diagram
ids_MGUS <- first_date_mgus %>%
  pull(ID) %>%
  unique()

ids_MM <- merged_filt %>%
  pull(ID) %>%
  unique()


venn.plot <- draw.pairwise.venn(
  area1 = length(ids_MGUS),
  area2 = length(ids_MM),
  cross.area = length(intersect(ids_MGUS, ids_MM)),
  category = c("MG", "MM"),
  fill = c("steelblue", "firebrick"),
  alpha = 0.7,
  cat.cex = 1.8,
  cex = 1.84,
  scaled = TRUE
)

grid.newpage()
grid.draw(venn.plot)

ggsave(venn.plot, filename="Venn_MGUS_MM.pdf", width=10, height=8)

#merge mgus and mm
final_db<-left_join(first_date_mgus, merged_filt, by="ID")
#PFS_event
final_db<-final_db%>%mutate(PFS_event=ifelse(!is.na(Date_MM), "progression", "none"))
# -------------------------------------------------------------------------
#Censoring date
#For cancer it is based on cancer registry
# -------------------------------------------------------------------------
#load origin
origin<-read_csv("/mnt/project/extract_fields_ttyd/origin.csv")

origin <- origin %>%
  mutate(
    country = case_when(
      `UK Biobank assessment centre | Instance 0` %in% c("Cardiff", "Swansea", "Wrexham") ~ "Wales",
      `UK Biobank assessment centre | Instance 0` %in% c("Edinburgh", "Glasgow") ~ "Scotland",
      TRUE ~ "England"
    )
  )

origin<-origin%>%select(`Participant ID`, country)
colnames(origin)<-c("ID", "country")

final_db<-final_db%>%left_join(origin, by="ID")


# -------------------------------------------------------------------------
#Add censoring date fro Cancer registry based on country
# -------------------------------------------------------------------------

final_db <- final_db %>%
  mutate(
    censoring_date = 
      case_when(
        country == "England"  ~ as.Date("2023-05-31"),
        country == "Scotland" ~ as.Date("2023-09-30"),
        country == "Wales"    ~ as.Date("2016-12-31")
      )  )


# -------------------------------------------------------------------------
#Add Death date
# -------------------------------------------------------------------------
death_db<-read_csv("/mnt/project/extract_fields_ttyd/death_register.csv")
death_mgus<-death_db%>%filter(`Participant ID`%in%unique(final_db$ID))
death_mgus<-death_mgus[c(1,2)]
colnames(death_mgus)<-c("ID", "Date_Death")
final_db<-final_db%>%left_join(death_mgus, by="ID")

# censoring dates are registry-specific and depend on region and outcomes
# Deaths occurring after the registry censoring date are treated as censored (Event = 0),
# because cancer outcomes after the censoring date are unknown.

final_db <- final_db %>%
  mutate(
    PFS_event = case_when(
      # 1 = Multiple Myeloma diagnosis
      !is.na(Date_MM) ~ 1,
      
      # 2 = Death before administrative censoring
      is.na(Date_MM) &
        !is.na(Date_Death) &
        Date_Death <= censoring_date ~ 2,
      
      # 0 = Censored
      TRUE ~ 0
    )
  )

final_db <- final_db %>%
  mutate(
    PFS_date = case_when(
      PFS_event == 0 ~ censoring_date,
      PFS_event == 1 ~ Date_MM,
      PFS_event == 2 ~ Date_Death,
      TRUE ~ as.Date(NA)
    )
  )

final_db <- final_db %>%
  mutate(PFS_mos=(as.numeric(PFS_date-Date_MGUS))/30)


p <- ggplot(final_db, aes(x = PFS_mos)) +
  geom_histogram(
    binwidth = 6,
    fill = "#9ecae1",
    color = "black",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 12,
    linetype = "dashed",
    linewidth = 0.8,
    color = "#de2d26"
  ) +
  annotate(
    "text",
    x = 12,
    y = Inf,
    label = "12 months",
    vjust = 1.3,
    hjust = -0.05,
    color = "#de2d26",
    size = 5
  ) +
  labs(
    x = "Time from MGUS to MM, death, or end of follow-up (months)",
    y = "Number of individuals"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.margin = margin(10, 20, 10, 10)
  ) +
  coord_cartesian(clip = "off")

p
ggsave(p, filename="~/distribution_follow_up_times.pdf", width = 8, height = 6)

write_csv(final_db,"~/final_db_january26.csv")


# -------------------------------------------------------------------------
#Load CHIP calls
# -------------------------------------------------------------------------
CHIP<-read_csv("/mnt/project/extract_fields_ttyd/CHIP_calls.csv")
CHIP_filt<-CHIP%>%filter(`Participant ID`%in%unique(final_db$ID))
CHIP_filt$ID<-CHIP_filt$`Participant ID`



# -------------------------------------------------------------------------
#Format CHIP df
# -------------------------------------------------------------------------
get_k <- function(x) str_match(x, regex("Array\\s*(\\d+)", ignore_case = TRUE))[,2]

make_chip_long2 <- function(df,
                            id_col = "Participant ID",
                            nvar_col = "Clonal haematopoiesis of indeterminate potential (CHIP) number of variants") {
  
  cn <- names(df)
  
  variant_cols <- cn[str_detect(cn, regex("\\bCHIP\\b.*\\bvariant\\b.*Array", ignore_case = TRUE))] %>%
    setdiff(cn[str_detect(cn, regex("\\bVAF\\b", ignore_case = TRUE))])  # evita di prendere anche VAF
  
  vaf_cols <- cn[str_detect(cn, regex("\\bvariant allele frequency\\b|\\bVAF\\b", ignore_case = TRUE)) &
                   str_detect(cn, regex("Array", ignore_case = TRUE))]
  
  if (length(variant_cols) == 0) stop("Non trovo colonne CHIP variant | Array k.")
  if (length(vaf_cols) == 0) stop("Non trovo colonne VAF | Array k.")
  
  # ---- pivot VARIANT (character) ----
  var_long <- df %>%
    select(all_of(c(id_col, nvar_col)), all_of(variant_cols)) %>%
    pivot_longer(
      cols = all_of(variant_cols),
      names_to = "variant_col",
      values_to = "variant"
    ) %>%
    mutate(array_k = get_k(variant_col)) %>%
    select(-variant_col)
  
  # ---- pivot VAF (numeric) ----
  vaf_long <- df %>%
    select(all_of(c(id_col)), all_of(vaf_cols)) %>%
    pivot_longer(
      cols = all_of(vaf_cols),
      names_to = "vaf_col",
      values_to = "vaf"
    ) %>%
    mutate(array_k = get_k(vaf_col)) %>%
    select(-vaf_col)
  
  # ---- join + gene ----
  out <- var_long %>%
    left_join(vaf_long, by = c(id_col, "array_k")) %>%
    mutate(
      gene = str_extract(as.character(variant), "^[^:]+"),
      vaf  = suppressWarnings(as.numeric(vaf))
    ) %>%
    filter(!is.na(variant), variant != "", !is.na(gene)) %>%
    arrange(.data[[id_col]], as.integer(array_k))
  
  out
}

chip_long <- make_chip_long2(CHIP_filt)
chip_wide_summary <- chip_long %>%
  group_by(`Participant ID`) %>%
  summarise(
    n_variants = first(`Clonal haematopoiesis of indeterminate potential (CHIP) number of variants`),
    gene_vaf = paste0(gene, "=", vaf) %>% paste(collapse = "; "),
    .groups = "drop"
  )

chip_wide_summary<-chip_wide_summary%>%rename(ID="Participant ID")
final_db<-final_db%>%left_join(chip_wide_summary, by="ID")

final_db<-final_db%>%mutate(n_variants=ifelse(is.na(n_variants), 0, n_variants))
final_db<-final_db%>%mutate(CHIP_bi=ifelse(n_variants>0, "YES", "NO"))
final_db<-final_db%>%filter(PFS_mos>12)


df <- final_db %>%
  mutate(
    PFS_mos = as.numeric(PFS_mos),
    PFS_event   = as.integer(PFS_event),
    # cuminc vuole factor con 1° livello = censored
    status  = factor(case_when(
      PFS_event == 0 ~ "censored",
      PFS_event == 1 ~ "PFS_event",
      PFS_event == 2 ~ "competing",
      TRUE ~ NA_character_
    ), levels = c("censored", "PFS_event", "competing")),
    CHIP_bi = factor(CHIP_bi)   # meglio come factor
  ) %>%
  filter(!is.na(PFS_mos), !is.na(status), !is.na(CHIP_bi))

ci <- tidycmprsk::cuminc(Surv(PFS_mos, status) ~ CHIP_bi, data = df)
ci

p_prog  <- ci$cmprsk$Tests[1, "pv"]  # cause 1 = progression
p_death <- ci$cmprsk$Tests[2, "pv"]  # cause 2 = death


p<-ggsurvfit::ggcuminc(ci, outcome = "PFS_event") +
  labs(
    x = "Months",
    y = "Cumulative incidence Progression",
    color = "CHIP status"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste0("Gray p = ", scales::pvalue(p_prog, accuracy = 0.001)),
    hjust = 1.05, vjust = -0.8,
    size = 5
  ) +
  theme_classic(base_size = 16) +
  add_risktable()+
  scale_color_manual(values = c("#808080", "#1f78b4"))

ggsave(p, filename="cumulative_risk_CHIP_UKB.pdf", width=8, height=6)


# -------------------------------------------------------------------------
#look at distribution of delta
# -------------------------------------------------------------------------
distribution_delta<-ggplot(final_mgus_data, aes(x=delta/30))+geom_histogram()+
  labs(x="Months from MGUS to MM", y="n")+geom_vline(xintercept = 12, color="red", linetype="dashed")+theme_classic()
ggsave(distribution_delta, filename="~/distribution_delta_mgus.pdf", width = 8, height = 6)



