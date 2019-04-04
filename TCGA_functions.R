library(GenomicDataCommons)
TCGAtranslateID = function(file_ids, legacy = TRUE) {
  info = files(legacy = T) %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  id_list = lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  barcodes_per_file = sapply(id_list,length)
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}

library(TCGAbiolinks)
#https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html#get_clinical_indexed_data
###Clinical Data
query <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Clinical", 
                  file.type = "xml")
GDCdownload(query)
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
clinical.folow_up <- GDCprepare_clinic(query, clinical.info = "follow_up")
clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
clinical.patient <- GDCprepare_clinic(query, clinical.info = "patient")
clinical.stage_event <- GDCprepare_clinic(query, clinical.info = "stage_event")
clinical.new_tumor_event <- GDCprepare_clinic(query, clinical.info = "new_tumor_event")

###Biospecimen data
query <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Biospecimen", 
                  file.type = "xml")
GDCdownload(query)
biospecimen.sample <- GDCprepare_clinic(query, clinical.info = "sample")
biospecimen.bio_patient <- GDCprepare_clinic(query, clinical.info = "bio_patient")


####Download all the clinical data 
library(data.table)
library(dplyr)
library(regexPipes)
clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA-PAAD",value=T) %>% 
  sort %>% 
  plyr::alply(1,GDCquery_clinic, .progress = "text") %>% 
  rbindlist

