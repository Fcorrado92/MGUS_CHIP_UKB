library(renv)
renv::init()
install.packages("readr")
install.packages("survminer")
install.packages("cmprsk")
install.packages("tidycmprsk")

library(readr)
library(dplyr)

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

mgus_data_filt_first_date <- mgus_data_filt_first_date %>%
  mutate(
    censoring_date = case_when(
      Data_Provider == "England (TPP)"     ~ as.Date("2016-05-31"),
      Data_Provider == "England (Vision)"  ~ as.Date("2017-05-31"),
      Data_Provider == "Scotland"          ~ as.Date("2017-03-31"),
      Data_Provider == "Wales"             ~ as.Date("2017-08-31"),
      TRUE ~ as.Date(NA)
    )
  )



# -------------------------------------------------------------------------
#select MM cases
# -------------------------------------------------------------------------


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

encounter_n_mm<-mm_data%>%group_by(`Participant ID`)%>%summarise(n=n_distinct(`Date clinical code was entered`))
ggplot(encounter_n_mm, aes(x=n))+geom_histogram()
 table(encounter_n_mm$n)

 # 1   2   3   4   5   6   7   8   9  10  12  14  16  26  27  28  34  36  44  67 
 # 182  64  26  18   7   3   1   3   2   3   4   1   1   1   1   1   1   1   1   1  

#Select only cols ID and Date and then keep only first date
mm_data_filt<-mm_data%>%select(`Participant ID`,`Date clinical code was entered`)
colnames(mm_data_filt)<-c("ID", "Date_MM")

mm_data_filt_first_date <- mm_data_filt %>%
  group_by(ID) %>%
  arrange(Date_MM) %>%          # sort by mm date (earliest first)
  slice(1) %>%                  # keep only the first row per participant
  ungroup()


# -------------------------------------------------------------------------
#left join MGUS and MM
# -------------------------------------------------------------------------

mgus_data_filt_first_date<-mgus_data_filt_first_date%>%left_join(mm_data_filt_first_date, by="ID")

mgus_data_filt_first_date<-mgus_data_filt_first_date%>%mutate(delta=as.Date(Date_MM)-as.Date(Date_MGUS))



# -------------------------------------------------------------------------
#save dataframe and plot
# -------------------------------------------------------------------------
mgus_data_filt_first_date$codes<-"GP"

write_csv(mgus_data_filt_first_date,"mgus_data_GP.csv")



