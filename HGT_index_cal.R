#Copyright (c) 2021 Taro Maeda
#This software is released under the MIT License, see LICENSE.

library(tidyverse)

### Usage
## Rscript HGT_index.R dataA1, dataA2, dataB1, dataB2, type_name
## dataA1-2: blast output (fmt6 format) for hA index calculation (e.g. A1=result against algal database, A2=result against animal database)
## dataB1-2: blast output (fmt6 format) for h index calculation (e.g. B1=result against Prokaryote database, B2=result against Eukaryote database)


my_tophit<- function(data)
{
    blast_ftm6<-c(
      "qseqid",
      "sseqid",
      "pident",
      "length",
      "mismatch",
      "gapopen",
      "qstart",
      "qend",
      "sstart",
      "send",
      "evalue",
      "bitscore"
    )
    
    
DATA<-
  read_tsv(data, col_names = blast_ftm6, guess_max=10000)

top1<-
  DATA %>%
  group_by(qseqid) %>%
  top_n(-1, evalue) %>%
  ungroup()

return(top1)
}

my_mod <- function(dataA1, dataA2, dataB1, dataB2, type_name)
{
  

  
  DATA_A1<-
  my_tophit(dataA1)
  
  DATA_A1_mod<-
  DATA_A1 %>%
    select(qseqid, evalue, bitscore)%>%
    rename(dataA1_evalue=evalue, dataA1_bitscore=bitscore)%>%
    distinct(qseqid, .keep_all = TRUE)
  
  
  DATA_A2<-
    my_tophit(dataA2)
  
  DATA_A2_mod<-
    DATA_A2 %>%
    select(qseqid, evalue, bitscore)%>%
    rename(dataA2_evalue=evalue, dataA2_bitscore=bitscore)%>%
    distinct(qseqid, .keep_all = TRUE)
  
  DATA_B1<-
    my_tophit(dataB1)
  
  DATA_B1_mod<-
    DATA_B1 %>%
    select(qseqid, evalue, bitscore)%>%
    rename(dataB1_evalue=evalue, dataB1_bitscore=bitscore)%>%
    distinct(qseqid, .keep_all = TRUE)
  
  
  DATA_B2<-
    my_tophit(dataB2)
  
  DATA_B2_mod<-
    DATA_B2 %>%
    select(qseqid, evalue, bitscore)%>%
    rename(dataB2_evalue=evalue, dataB2_bitscore=bitscore)%>%
    distinct(qseqid, .keep_all = TRUE)
  
  
  DATA_A1_A2_B1_B2<-
    DATA_A1_mod %>%
    full_join(DATA_A2_mod, by=c("qseqid")) %>%
    full_join(DATA_B1_mod, by=c("qseqid")) %>%
    full_join(DATA_B2_mod, by=c("qseqid")) 
    
  DATA_A1_A2_B1_B2_mod<-
  DATA_A1_A2_B1_B2%>%
    replace_na((list(
                    dataA1_bitscore = 0, 
                    dataA2_bitscore=0,
                    dataB1_bitscore=0,
                    dataB2_bitscore=0)
                ))

  out<-
  DATA_A1_A2_B1_B2_mod %>%
  mutate(hgt_dataA_index=dataA1_bitscore-dataA2_bitscore) %>%
  mutate(hgt_dataB_index=dataB1_bitscore-dataB2_bitscore) %>%
  mutate(type=type_name)

  #return(poc_animal_algae)
  return(out)
}


my_plot <- function(g_data)
{
  return(
    ggplot(g_data, aes(x=hgt_dataA_index, y=hgt_dataB_index))+
      geom_point(alpha = .5, size = .6) + theme_bw()+ 
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      geom_hline(yintercept=-100, colour="red",linetype="dashed") +
      geom_vline(xintercept=100, colour="red",linetype="dashed") +
      xlim(-1000, 1000)+
      ylim(-1000, 1000)+
      labs(x="dataA2<- hA_index ->dataA1", y="Euk (dataB2)<- h_index ->Bac (dataB1)")

  )
}


########main

args <- commandArgs(trailingOnly=TRUE)

dataA1_path <- args[1]
dataA2_path <- args[2]
dataB1_path <- args[3]
dataB2_path <- args[4]
Title <-args[5]

output_d<-
  my_mod(
    dataA1_path,
    dataA2_path,
    dataB1_path,
    dataB2_path,
    Title
  )

write_tsv(output_d, "HGT_index.out")

output_d_g<-
my_plot(output_d)+
      labs(title=Title)


ggsave(file = "HGT_index_plot.pdf", plot = output_d_g)
