death_db<-read_csv("/mnt/project/extract_fields_ttyd/death_register.csv")
mgus_db<-read_csv("/mnt/project/outputs/mgus_data_CR_GP.csv")
death_mgus<-death_db%>%filter(`Participant ID`%in%unique(mgus_db$ID))
death_mgus<-death_mgus[c(1,2)]
colnames(death_mgus)<-c("ID", "Date_Death")
mgus_db<-mgus_db%>%left_join(death_mgus, by="ID")

mgus_db<-mgus_db%>%mutate(Event=ifelse(is.na(Date_MM), "0", "1"))
mgus_db<-mgus_db%>%mutate(Event=ifelse(!is.na(Date_Death)&is.na(Date_MM), "2", Event))

region