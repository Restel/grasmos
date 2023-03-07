library(plyr)
library(dplyr)
library(tidyverse)

# subtiWIKI

df_subti = read.csv("data_source/Subtiwiki.csv")
df_subti <- select(df_subti, regulator.name, mode, gene.name)  %>%
          filter(!mode%in% c("[wiki|RNA switch]", 
                             "indirect effect",
                             "antisense RNA",
                             "attenuation",
                             "mRNA stability control", 
                             "RNA switch", 
                             "termination/ antitermination", 
                             "termination/antitermination",
                             "transcription termination/ antitermination",
                             "transcription termination/antitermination",
                             "translation control",
                             "unknown"))%>%
          filter(!regulator.name %in% c("stringent response")) %>%
          dplyr::rename(source = regulator.name,target = gene.name ) %>%
          mutate(mode = as.factor(mode))
df_subti = df_subti %>% mutate(mode = plyr::revalue(mode, c("activation "= "activation", 
                                                      "anti-activation" = "repression",
                                                      "anti-termination" = "activation",
                                                      "antitermination" = "activation",
                                                      "auto-repression" = "repression",
                                                      "autorepression" = "repression",
                                                      "indirect positive regulation" = "activation",
                                                      "negative autoregulation" = "repression",
                                                      "processive antitermination" = "repression",
                                                      "sigma factor" = "activation",
                                                      "Sigma factor" = "activation",
                                                      "termination" = "repression",
                                                      "transcription activation" = "activation",
                                                      "transcription antitermination" = "activation",
                                                      "transcription repression" = "repression",
                                                      "transcription termination" = "repression",
                                                      "transcriptional activation" = "activation",
                                                      "transcriptional antitermination" = "activation"
                                                      )
                                              )
                         )
df_subti$id = paste0(df_subti$source,df_subti$target,df_subti$mode)
df_subti = df_subti[!duplicated(df_subti$id),]
df_subti = select(df_subti, -c(id))
write.csv(df_subti, file = "data/real/df_subti.csv", row.names = FALSE)

#### REGULON

df_regulon = read.delim("data_source/network_tf_gene.txt", skip = 38, header = FALSE)

df_regulon <- df_regulon %>% select(V2, V5, V6)
df_regulon <- df_regulon %>% 
              dplyr::rename(source = V2, target = V5, mode = V6) %>%
              mutate(mode = plyr::revalue(mode, c("-"= "repression", "+" = "activation","?" = "unknown"))) %>%
              mutate(source = sapply(source, function(s) {
                                                  first_letter = tolower(substr(s, 1, 1)) 
                                                  return(paste0(first_letter, substr(s, 2, nchar(s))))} )) %>%
              mutate(mode = as.factor(mode))
              
df_regulon$source[df_regulon$source == 'CRP']='crp'
df_regulon$id = paste0(df_regulon$source,df_regulon$target)
df_regulon = df_regulon[!duplicated(df_regulon$id)&!duplicated(df_regulon$id, fromLast = TRUE),]
df_regulon = select(df_regulon, -c(id))
write.csv(df_regulon, file = "data/real/df_regulon.csv", row.names = FALSE)
