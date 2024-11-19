################################################################################
### Script: Functions for figure consolidation #################################
### Author: Wouter Maessen, Ana Conesa Lab #####################################
### Date: 14-10-2024 ###########################################################
################################################################################

list_of_packages <- c("ggplot2", "tidyverse", "ggpubr", "scales", "patchwork", "gridExtra", "grid", "RColorConesa", "fmsb", "MetBrewer")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")


# Install and load packages
library(ggplot2)
library(ggpubr)
library(scales)
library(patchwork)
library(gridExtra)
library(grid)
library(RColorConesa)
library(fmsb)
library(MetBrewer)

outdir = "output/main"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

#### set theme for plots
pub_theme <- theme_pubclean(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y  = element_text(vjust=0.5, size=14) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=10), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=14)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  theme(legend.position = "none")

old.libplat.palette = c( "cDNA+ONT"="#74CDF0", "cDNA+PacBio"="#EE446F", "cDNA+Illumina"="#FFCF71", 
                         "CapTrap+ONT"="#7482F0", "R2C2+ONT"="#74F0D9", "dRNA+ONT"="#13BF5E", "CapTrap+PacBio"="#d14141")

libplat.palette = c( "cDNA-PacBio"="#c06636", "CapTrap-PacBio"="#802417", "cDNA-Illumina"="#e8b960",  "Freestyle-Freestyle"="#ce9344",
                     "cDNA-ONT"="#646e3b", "CapTrap-ONT"="#17486f", "R2C2-ONT"="#508ea2", "dRNA-ONT"="#2b5851"
)
cat.palette = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")
print(getwd())
setwd(rdata)

# Load existing data
# Load and modify ES_code
ES_code <- read.csv("Challenge3_Figures_Data/ES/ES_code.txt", sep=",", header = T )
ES_summary_table <- read.csv("Challenge3_Figures_Data/ES_challenge1/ES_challenge1_metrics.summary_table_SC.csv", header = T)
manatee_code <- read.csv("Challenge3_Figures_Data/manatee/manatee_code.txt", sep=",", header = T )
manatee_metrics <- read.csv("Challenge3_Figures_Data/manatee/manatee_challenge3_metrics.challenge3_metrics.csv", sep=",",header = T) %>% t()
ES_metrics_perc <- read.csv("Challenge3_Figures_Data/ES/ES_challenge3_metrics.challenge3_metrics_perc.csv", sep=",",header = T) %>% t()
manatee_SIRV_metrics <- read.csv("Challenge3_Figures_Data/manatee/manatee_challenge3_metrics.SIRVS_metrics.csv", sep=",",header = T) %>% t()

# Check correct structure of meta_data
if(!all(colnames(manatee_code) == colnames(meta_data))){
  print('ERROR: The column names of the meta_data is expected to match the column names of manatee_code: "pipelineCode" "Library_Preps" "Platform" "Data_Category" "Lab" "Tool" "Alias"')
  break
} 

manatee_code <- rbind(manatee_code, meta_data)

manatee_code$Lib_Plat <- apply(manatee_code, 1, function(x){
  paste(x["Library_Preps"], x["Platform"], sep = "-")
})
manatee_code$Lib_DC=apply(cbind(manatee_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
manatee_code$Label <-apply(cbind(manatee_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")

ES_code <- rbind(ES_code, meta_data)

ES_code$Lib_Plat <- apply(ES_code, 1, function(x){
  paste(x["Library_Preps"], x["Platform"], sep = "-")
})
ES_code$Lib_DC=apply(cbind(ES_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
ES_code$Label <-apply(cbind(ES_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")


### Script for figure 4a
print('Script for figure 4a does not yet exist...')

### Script for figure 4b
append_new_data_4b <- function(new_data){ # is num_of_monoexons necessary?
  
  # Check correct structure of new_data
  if(ncol(new_data) != 10){
    print('ERROR: The number of columns in the provided data is not correct: Expected number of columns = 10')
    break
  } else {
    print('Appending new data to existing data...')
  }
  
  manatee_metrics <- rbind(manatee_metrics, new_data)
  rownames(manatee_metrics)[nrow(manatee_metrics)] <- as.character(manatee_code$pipelineCode[nrow(manatee_code)])
  pipelineCode <- meta_data$pipelineCode
  tool_name <- meta_data$Tool
  
  mono_exons_df<-data.frame(Row.names=c("ONT1","PB1","ONT2", "illumina1","ONT3","ONT4","PB2","PB3","ONT5", manatee_code$pipelineCode[nrow(manatee_code)]), 
                            Num.monoexons=c(0,7,79174,635028,3703,470801,52,246565,85291, sum(temp_env$sqanti_data$subcategory == 'mono-exon')))
  rownames(mono_exons_df) <- mono_exons_df$Row.names
  manatee_metrics <- merge(manatee_metrics, mono_exons_df,by.x=0, by.y="Row.names")
  manatee_metrics$Num.multiexons <- apply(manatee_metrics,1, function(x){
    as.numeric(x["Number of transcripts"]) - as.numeric(x["Num.monoexons"])
  })
  manatee_metrics <- merge(manatee_metrics, manatee_code, by.x="Row.names", by.y="pipelineCode")
  
  manatee_metrics$Tool <- gsub("-", "\n", manatee_metrics$Tool)
  colnames(manatee_metrics) <- make.names(colnames(manatee_metrics))
  
  melted_manatee <- manatee_metrics %>%
    pivot_longer(c("Num.monoexons", "Num.multiexons"))
  melted_manatee$Tool <- gsub("-", "\n", melted_manatee$Tool)
  
  melted_manatee <- melted_manatee %>% filter(Row.names!="ONT1")
  melted_manatee$Row.names <- as.factor(melted_manatee$Row.names)

  # Extract the value from the data frame and format it as required
  number_of_tr_new_tool <- as.integer(manatee_metrics[manatee_metrics$Row.names == meta_data$pipelineCode,]$Number.of.transcripts)
  example_label <- sprintf("n=%.1fK", number_of_tr_new_tool / 1000)

  # Convert the Tool column to a factor with the desired levels
  melted_manatee$Tool <- factor(melted_manatee$Tool, levels = c(sort(setdiff(unique(melted_manatee$Tool), tool_name)), tool_name))
  
  pF4.1 <- ggplot(melted_manatee, aes(x=Row.names, y=value, alpha=name, fill=Lib_Plat))+
    geom_bar(stat="identity", position = "stack", width=0.7)+
    facet_grid( .~Tool, scales = 'free_x', space="free_x") +
    scale_fill_manual(values =  libplat.palette, name="Library-Platoform") +
    scale_alpha_manual(values=c(0.6,1), breaks=c("Num.monoexons", "Num.multiexons"), labels=c("Mono-exons", "Multi-exons"))+
    scale_x_discrete(breaks=c("PB1","ONT3","PB2","ONT2","ONT5","illumina1","ONT4", "PB3", meta_data$pipelineCode),
                     labels=c("n=1.9K","n=63K","n=25K","n=179K","n=177K","n=916K","n=543K","n=294K", example_label)) +
    pub_theme +
    ylab("Num. total isoforms") +
    xlab("")+
    theme(strip.text.x = element_text(size=18)) +
    scale_y_continuous(label = unit_format(unit = "K", scale = 0.001), expand = expansion(mult = c(0,0.1)))+
    theme(legend.position = "none") +
    theme(plot.margin=unit(c(0,0,-0.5,0), "cm"))
  
  # Count unique transcripts per gene
  transcripts_per_gene <- sqanti_data %>%
    group_by(associated_gene) %>%
    summarize(num_transcripts = n_distinct(associated_transcript))
  
  transcript_categories <- transcripts_per_gene %>%
    mutate(
      category = case_when(
        num_transcripts == 1 ~ "1",
        num_transcripts == 2 | num_transcripts == 3 ~ "2-3",
        num_transcripts == 4 | num_transcripts == 5 ~ "4-5",
        num_transcripts >= 6 ~ ">6"
      )
    ) %>%
    count(category)
  
  transcript_categories <- as.data.frame(transcript_categories)
  colnames(transcript_categories) <- c('Cat', 'Num_locus')
  
  # Process info number isoforms per locus
  trx_locus.list <- list()
  for (i in c("ONT1", "ONT2","ONT3","ONT4","ONT5","PB1","PB2","PB3","illumina1")){
    trx_file=paste0("Challenge3_Figures_Data/manatee/manatee_isoforms_per_gene/",i,"_cpm_vs_trans.tsv")
    df=read.csv(trx_file,sep="\t",header = T)
    num_trx_locus <- data.frame(Cat=c("1","2-3","4-5", ">6"), 
                                Num_locus=c(which(df$n_isoforms==1) %>% length(), 
                                            which(df$n_isoforms==2 |  df$n_isoforms==3) %>% length(),
                                            which(df$n_isoforms==4 |  df$n_isoforms==5) %>% length(),
                                            which(df$n_isoforms>5) %>% length()
                                )
    )
    trx_locus.list[[i]] <- num_trx_locus
  }
  
  trx_locus.list[[length(trx_locus.list) + 1]] <- transcript_categories
  names(trx_locus.list)[[length(trx_locus.list)]] <- pipelineCode
  
  trx_locus.df <- bind_rows(trx_locus.list, .id = "pipeline")
  trx_locus.df <- merge(trx_locus.df, manatee_code, by.x="pipeline", by.y="pipelineCode")
  trx_locus.df$Tool <- gsub("-", "\n", trx_locus.df$Tool)
  
  # Convert the Tool column to a factor with the desired levels
  trx_locus.df$Tool <- factor(trx_locus.df$Tool, levels = c(sort(setdiff(unique(trx_locus.df$Tool), tool_name)), tool_name))
  
  trx_locus.df <- trx_locus.df %>% filter(pipeline!="ONT1")
  pF4.2 <- ggplot(trx_locus.df, aes(x=pipeline, y=Num_locus, fill=Cat)) + 
    geom_bar(position = "fill", stat="identity", width = 0.7) +
    facet_grid( .~Tool , scales="free_x", space="free_x", drop=TRUE) +
    scale_fill_met_d(name = "Cassatt1", direction = -1) +
    scale_x_discrete(breaks=c("PB1","ONT3","PB2","ONT2","ONT5","illumina1","ONT4", "PB3", pipelineCode),
                     labels=c("LO", "LO", "LO", "LO", "LS","SO", "LS","LS", data_category)) +
    ylab("Number of transcripts \n per gene")+
    pub_theme +
    theme( strip.text.x = element_blank(),
           axis.title.x = element_blank(),
           legend.text = element_text(size=14),
           legend.title = element_blank()) +
    scale_y_continuous(label = unit_format(unit = "%", scale = 100), expand = expansion(mult = c(0,0.1)))+
    theme(legend.position = "bottom") +
    theme(plot.margin=unit(c(-1,0,0,0), "cm"))
  
  pF <- pF4.1 / pF4.2 +
    plot_layout(heights = c(1, 1), ncol = 1)
  ggsave(file=paste0(outdir, "/panel4b.svg"), plot=pF, width=12, height=6)
  
}

### Script for figure 4c
append_new_data_4c <- function(new_data){
  
  ######### length distribution mouse
  dist.list_ES <- list()
  for (i in c("ONT_1", "ONT_2","ONT_3","ONT_4","ONT_5",
              "ONT_6", "ONT_7","ONT_8","ONT_9","ONT_10",
              "ONT_11","PB_1","PB_2","PB_3",
              "PB_4", "PB_5", "illumina_1")){
    length_file=paste0("Challenge3_Figures_Data/ES/length_distributions/",i,"_length_dist.tsv")
    df=read.csv(length_file,sep="\t",header = T)
    dist.list_ES[[i]] <- as.numeric(df$length) %>% as.data.frame()
  }
  
  dist.list_ES[[ES_code$pipelineCode[nrow(ES_code)]]] <- as.numeric(new_data$length) %>% as.data.frame()
  dist_df_ES <- bind_rows(dist.list_ES, .id = "pipeline")
  colnames(dist_df_ES)<-c("pipeline","length")
  dist_df_ES <- merge(dist_df_ES, ES_code, by.x="pipeline", by.y="pipelineCode")
  dist_df_ES$Tool <- gsub("-", "\n", dist_df_ES$Tool)
  
  label_breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                   "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS", ES_code$Label[nrow(ES_code)])
  labels = c("SO", rep("LO", 6), rep("LS", 3), data_category)
  
  # print counts to use in caption, order to match plots
  count_df_ES <- dist_df_ES %>%
    group_by(pipeline, Label, Tool) %>%
    summarise(Count = n(), .groups = 'drop')
  count_df_ES <- count_df_ES[order(tolower(count_df_ES$Tool), tolower(count_df_ES$Label)),]
  print("== Panel 4c ES counts ===")
  print(count_df_ES)
  
  # Convert the Tool column to a factor with the desired levels
  dist_df_ES$Tool <- factor(dist_df_ES$Tool, levels = c(sort(setdiff(unique(dist_df_ES$Tool), tool_name)), tool_name))
  
  pF3 <- ggplot(dist_df_ES, aes(x=Label, y=length, fill=Lib_Plat))+
    geom_violin()+
    geom_boxplot(width=0.15, color="white", alpha=0.8, outlier.shape = NA) +
    facet_grid( .~Tool , scales="free_x", space="free_x") +
    scale_fill_manual(values =  libplat.palette, name="Library-Platoform") +
    scale_x_discrete(breaks=label_breaks, labels=labels)+
    xlab("")+ 
    ylab("Length (bp), log10")+
    pub_theme +
    theme(legend.position = "none") +
    theme(strip.text.x = element_text(size=16)) +
    scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)) )
  
  ggsave(file=paste0(outdir, "/panel4c.svg"), plot=pF3, width=12, height=6)
  
} 

### Script for figure 4d
append_new_data_4d <- function(new_data){
  
  dist.list_manatee <- list()
  for (i in c("ONT1", "ONT2","ONT3","ONT4","ONT5","PB1","PB2","PB3","illumina1")){
    length_file=paste0("Challenge3_Figures_Data/manatee/length_distributions/",i,"_length_dist.tsv")
    df=read.csv(length_file,sep="\t",header = T)
    dist.list_manatee[[i]] <- as.numeric(df$length) %>% as.data.frame()
  }
  
  dist.list_manatee[[manatee_code$pipelineCode[nrow(manatee_code)]]] <- as.numeric(new_data$length) %>% as.data.frame()
  dist_df_manatee <- bind_rows(dist.list_manatee, .id = "pipeline")
  colnames(dist_df_manatee)<-c("pipeline","length")
  dist_df_manatee <- merge(dist_df_manatee, manatee_code, by.x="pipeline", by.y="pipelineCode")
  dist_df_manatee$Tool <- gsub("-", "\n", dist_df_manatee$Tool)
  
  sample_size = dist_df_manatee %>% group_by(pipeline) %>% summarize(num=n())
  label_breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                   "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS", manatee_code$Label[nrow(manatee_code)])
  labels = c("SO", rep("LO", 6), rep("LS", 3), data_category)
  
  # print counts to use in caption, order to match plots
  count_df_manatee <- dist_df_manatee %>%
    group_by(pipeline, Label, Tool) %>%
    summarise(Count = n(), .groups = 'drop')
  count_df_manatee <- count_df_manatee[order(tolower(count_df_manatee$Tool), tolower(count_df_manatee$Label)),]
  print("== Panel 4d manatee counts ===")
  print(count_df_manatee)
  
  # Convert the Tool column to a factor with the desired levels
  dist_df_manatee$Tool <- factor(dist_df_manatee$Tool, levels = c(sort(setdiff(unique(dist_df_manatee$Tool), tool_name)), tool_name))
  
  pF2 <- ggplot(dist_df_manatee, aes(x=Label, y=length, fill=Lib_Plat))+
    geom_violin()+
    geom_boxplot(width=0.15, color="white", alpha=0.8, outlier.shape = NA) +
    facet_grid( .~Tool , scales="free_x", space="free_x") +
    scale_fill_manual(values =  libplat.palette, name="Library-Platoform") +
    scale_x_discrete(breaks=label_breaks, labels=labels)+
    xlab("")+ 
    ylab("Length (bp), log10")+
    pub_theme +
    theme(legend.position = "none") +
    theme(strip.text.x = element_text(size=16)) +
    scale_y_continuous(trans='log10', breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x)) )
  
  ggsave(file=paste0(outdir, "/panel4d.svg"), plot=pF2, width=12, height=6)
  
} 

### Script for figure 4e1

append_new_data_4e1 <- function(new_data){
  
  suppressWarnings(ES_metrics_perc <- rbind(ES_metrics_perc, new_row))
  rownames(ES_metrics_perc)[nrow(ES_metrics_perc)] <- pipelineCode
  
  ES_code <- read.csv("Challenge3_Figures_Data/ES/ES_code.txt", sep=",", header = T )
  ES_code <- rbind(ES_code, meta_data)
  ES_code$Lib_Plat <- apply(ES_code, 1, function(x){
    paste(x["Library_Preps"], x["Platform"], sep = "-")
  })
  ES_code$Lib_DC=apply(cbind(ES_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  ES_code$Label <-apply(cbind(ES_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")
  
  ES_metrics_perc <- merge(ES_metrics_perc, ES_code, by.x=0, by.y="pipelineCode")
  ES_metrics_perc$Tool <- gsub("-", "\n", ES_metrics_perc$Tool)
  colnames(ES_metrics_perc) <- make.names(colnames(ES_metrics_perc))
  ES_metrics_perc[13, "Mapping.transcripts"] <- 100
  
  # Convert the Tool column to a factor with the desired levels
  ES_metrics_perc$Tool <- factor(ES_metrics_perc$Tool, levels = c(sort(setdiff(unique(ES_metrics_perc$Tool), tool_name)), tool_name))
  
  pC_SJ <- ggplot(ES_metrics_perc, aes(x=Label, y=as.numeric(Splice.Junctions.with.short.read.coverage), color=Lib_Plat, shape=Data_Category)) + 
    geom_segment(aes(x=Label, xend=Label, y=0, yend=as.numeric(Splice.Junctions.with.short.read.coverage), color=Lib_Plat), size=2) +
    geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
    pub_theme+
    theme(axis.title.y = element_text(size=18),
          axis.text.y  = element_text(size=18) ) +
    scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
    scale_color_manual(values =  libplat.palette, name="Library-Platform") +
    scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                              "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                     labels=c("SO", rep("LO", 6), rep("LS", 3)))+
    facet_grid(.~ Tool, drop = TRUE, scales="free_x", space = "free_x") +
    xlab("") + ylab("% Splice Junctions without short read coverage")+
    scale_y_continuous(label = label_percent(scale = 1), limits=c(0,100), expand = expansion(mult = c(0,0.1)))+
    theme(legend.position = "none",
          strip.text.x = element_text(size = 18))
  
  ggsave(file=paste0(outdir, "/panel4e1.svg"), plot=pC_SJ, width=8, height=8)
  
} 

### Script for figure 4e2
print('Script for figure 4e2 does not yet exist...')

### Script for figure 4e3
append_new_data_4e3 <- function(new_data){
  
  suppressWarnings(ES_metrics_perc <- rbind(ES_metrics_perc, new_row))
  rownames(ES_metrics_perc)[nrow(ES_metrics_perc)] <- pipelineCode
  
  ES_code <- read.csv("Challenge3_Figures_Data/ES/ES_code.txt", sep=",", header = T )
  ES_code <- rbind(ES_code, meta_data)
  ES_code$Lib_Plat <- apply(ES_code, 1, function(x){
    paste(x["Library_Preps"], x["Platform"], sep = "-")
  })
  ES_code$Lib_DC=apply(cbind(ES_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  ES_code$Label <-apply(cbind(ES_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")
  
  ES_metrics_perc <- merge(ES_metrics_perc, ES_code, by.x=0, by.y="pipelineCode")
  ES_metrics_perc$Tool <- gsub("-", "\n", ES_metrics_perc$Tool)
  colnames(ES_metrics_perc) <- make.names(colnames(ES_metrics_perc))
  ES_metrics_perc[12, "Mapping.transcripts"] <- 100
  
  # Convert the Tool column to a factor with the desired levels
  ES_metrics_perc$Tool <- factor(ES_metrics_perc$Tool, levels = c(sort(setdiff(unique(ES_metrics_perc$Tool), tool_name)), tool_name))
  
  pC_nonCan <- ggplot(ES_metrics_perc, aes(x=Label, y=as.numeric(Non.canonical.Splice.Junctions), color=Lib_Plat, shape=Data_Category))  + 
    geom_segment( aes(x=Label, xend=Label, y=0, yend=as.numeric(Non.canonical.Splice.Junctions), color=Lib_Plat), size=2) +
    geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
    pub_theme+
    theme(axis.title.y = element_text(size=18),
          axis.text.y  = element_text(size=18) ) +
    scale_fill_manual(values =  libplat.palette, name="Library-Platform") +
    scale_color_manual(values =  libplat.palette, name="Library-Platform") +
    scale_x_discrete(breaks=c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO","CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                              "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS",  "cDNA-PacBio-LS"),
                     labels=c("SO", rep("LO", 6), rep("LS", 3)))+
    facet_grid(.~ Tool, drop = TRUE, scales="free_x", space = "free_x") +
    xlab("") + ylab("% Non-canonical Splice Junctions")+
    scale_y_continuous(label = label_percent(scale = 1), limits=c(0,100), expand = expansion(mult = c(0,0.1)))+
    theme(legend.position = "none",
          strip.text.x = element_text(size = 18))
  
  ggsave(file=paste0(outdir, "/panel4e3.svg"), plot=pC_nonCan, width=8, height=8)
  
} 

### Script for figure 4e4
print('Script for figure 4e4 does not yet exist...')

### Script for figure 4f
if(organism == 'mouse'){
  append_new_data_4f1_ES <- function(new_data_ES){
    
    ES_BUSCO <- read.csv("Challenge3_Figures_Data/ES/ES_challenge3_metrics.BUSCO_metrics.csv", sep=",",header = T) %>% t()
    new_ES_BUSCO <- array(NA, dim = c(nrow(ES_BUSCO) + 1, ncol(ES_BUSCO)))
    new_ES_BUSCO[1:nrow(ES_BUSCO), ] <- ES_BUSCO
    new_ES_BUSCO[nrow(new_ES_BUSCO), ] <- new_data_ES$`Absolute value`
    colnames(new_ES_BUSCO) <- colnames(ES_BUSCO)
    rownames(new_ES_BUSCO) <- c(rownames(ES_BUSCO), pipelineCode)
    ES_BUSCO <- new_ES_BUSCO
    
    ES_BUSCO <- merge(ES_BUSCO, ES_code, by.x=0, by.y="pipelineCode")
    ES_BUSCO$Tool <- gsub("-", "\n", ES_BUSCO$Tool)
    colnames(ES_BUSCO) <- make.names(colnames(ES_BUSCO))
    
    total_busco <- 11366
    get_perc_busco_found <- function(x){
      r=(total_busco - as.numeric(x["Missing.BUSCOs"]))*100/total_busco 
      round(r, digits = 2)
    }
    
    get_perc_busco_complete <- function(x){
      c <- as.numeric(x["Complete.and.single.copy.BUSCOs"]) + as.numeric(x["Complete.and.duplicated.BUSCOs"]) 
      r=c*100/total_busco 
      round(r, digits = 2)
    }
    
    get_perc_busco_fragmented <- function(x){
      c <- as.numeric(x["Fragmented.BUSCOs"]) 
      r=c*100/total_busco 
      round(r, digits = 2)
    }
    
    get_perc_busco_missing <- function(x){
      c <- as.numeric(x["Missing.BUSCOs"]) 
      r=c*100/total_busco 
      round(r, digits = 2)
    }
    
    ES_BUSCO$BUSCO_found <- apply(ES_BUSCO,1,get_perc_busco_found)
    ES_BUSCO$BUSCO_complete <- apply(ES_BUSCO,1,get_perc_busco_complete)
    
    ES_BUSCO$BUSCO_fragmented <- apply(ES_BUSCO,1,get_perc_busco_fragmented)
    ES_BUSCO$BUSCO_missing <- apply(ES_BUSCO,1,get_perc_busco_missing)
    
    ES_BUSCO$Tool <- factor(ES_BUSCO$Tool, levels = c(sort(setdiff(unique(ES_BUSCO$Tool), tool_name)), tool_name))
    
    # Modify `pD7.1` plot
    pD7.1 <- ggplot(ES_BUSCO, aes(x=Label, y=as.numeric(BUSCO_found), color=Lib_Plat, shape=Data_Category))  + 
      geom_segment(aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_found), color=Lib_Plat), size=2) +
      geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
      pub_theme +
      theme_pubclean(flip=TRUE) +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=14),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,0)  # Remove margins around the plot
      ) +
      scale_fill_manual(values = libplat.palette, name="Library-Platform") +
      scale_color_manual(values = libplat.palette, name="Library-Platform") +
      scale_x_discrete(breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO", "CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                                  "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS", "cDNA-PacBio-LS"),
                       labels = c("SO", rep("LO", 6), rep("LS", 3))) +
      facet_grid(Tool ~., drop = TRUE, scales = "free_y") +
      xlab("") + 
      theme(
        axis.title = element_text(size=16),
        strip.text.y = element_blank(), 
        axis.line.y = element_blank()
      ) +
      theme(legend.position = "none") +
      scale_y_continuous("", sec.axis = sec_axis(~ ., breaks = NULL, name = "Mouse ES"),
                         label = unit_format(unit = "%"), limits = c(-1,80), expand = expansion(mult = c(0, 0.1))) +
      coord_flip()
    
    # Calculate the appropriate height for the boxes in p.mid7 based on the number of levels
    num_tools <- length(unique(manatee_BUSCO$Tool))  # Number of unique tools
    ES_BUSCO$Tool <- factor(ES_BUSCO$Tool, levels = c(tool_name ,sort(setdiff(unique(ES_BUSCO$Tool), tool_name), decreasing = TRUE)))
    
    # Modify middle plot (p.mid7) to have consistent height for boxes
    p.mid7 <- ggplot(ES_BUSCO, aes(y = Tool)) + 
      geom_tile(aes(x = 1), fill = "gray", color = "white", width = 1.25, height = 0.95) +  # Adjust box height
      geom_text(aes(x = 1, label = Tool), size = 5) +                                   # Adjust text size
      theme_void() +                                                                    # Remove background and axes
      theme(
        plot.margin = margin(10,-60,25,0)  # Remove margins around the plot
      ) +
      scale_x_continuous(expand = c(0, 0)) +                                            # Remove space on the left
      coord_fixed(ratio = num_tools / 5)  # Adjust ratio to match the height of other plots
    
    # Combine the three plots
    gg1.7 <- ggplot_gtable(ggplot_build(pD7.1))
    g.mid <- ggplot_gtable(ggplot_build(p.mid7))
    
    # Arrange plots side by side with adjusted widths
    pD7 <- grid.arrange(g.mid, gg1.7, ncol = 2, widths = c(3/9, 6/9),
                        top = textGrob("BUSCO genes found (%)", gp = gpar(fontsize=18, font=1)))
    
    # Save the combined plot
    ggsave(file = paste0(outdir, "/panel4f1_ES.svg"), plot = pD7, width = 8, height = 5)
    
  }
  append_new_data_4f2_ES <- function(new_data_ES){
    
    ES_BUSCO$Tool <- factor(ES_BUSCO$Tool, levels = c(sort(setdiff(unique(ES_BUSCO$Tool), tool_name)), tool_name))
    
    pD7.1f <- ggplot(ES_BUSCO, aes(x=Label, y=as.numeric(BUSCO_fragmented), color=Lib_Plat, shape=Data_Category))  + 
      geom_segment(aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_fragmented), color=Lib_Plat), size=2) +
      geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
      pub_theme +
      theme_pubclean(flip=TRUE) +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=14),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,-20)  # Remove margins around the plot
      ) +
      scale_fill_manual(values = libplat.palette, name="Library-Platform") +
      scale_color_manual(values = libplat.palette, name="Library-Platform") +
      scale_x_discrete(breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO", "CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                                  "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS", "cDNA-PacBio-LS"),
                       labels = c("SO", rep("LO", 6), rep("LS", 3))) +
      facet_grid(Tool ~., drop = TRUE, scales = "free_y") +
      xlab("") + 
      theme(
        axis.title = element_text(size=16),
        strip.text.y = element_blank(), 
        axis.line.y = element_blank()
      ) +
      theme(legend.position = "none") +
      scale_y_continuous("", sec.axis = sec_axis(~ ., breaks = NULL, name = "Mouse ES"),
                         label = unit_format(unit = "%"), limits = c(-1, 20), expand = expansion(mult = c(0, 0.1))) +
      coord_flip()
    
    
    gg1.7f <- ggplot_gtable(ggplot_build(pD7.1f))
    
    pD7f <- grid.arrange(g.mid, gg1.7f,ncol=2,widths=c(3/9, 6/9),
                         top = textGrob("BUSCO genes fragmented (%)",gp=gpar(fontsize=18,font=1)))
    
    ggsave(file=paste0(outdir, "/panel4f2_ES.svg"), plot=pD7f, width=8, height=5)
    
    
  }
} else {
  # For using organism: Manatee
  append_new_data_4f1_manatee <- function(new_data_manatee){
    
    manatee_BUSCO <- read.csv("Challenge3_Figures_Data/manatee/manatee_challenge3_metrics.BUSCO_metrics.csv", sep=",",header = T) %>% t()
    new_manatee_BUSCO <- array(NA, dim = c(nrow(manatee_BUSCO) + 1, ncol(manatee_BUSCO)))
    new_manatee_BUSCO[1:nrow(manatee_BUSCO), ] <- manatee_BUSCO
    new_manatee_BUSCO[nrow(new_manatee_BUSCO), ] <- new_data_manatee$`Absolute value`
    colnames(new_manatee_BUSCO) <- colnames(manatee_BUSCO)
    rownames(new_manatee_BUSCO) <- c(rownames(manatee_BUSCO), pipelineCode)
    manatee_BUSCO <- new_manatee_BUSCO
    
    manatee_BUSCO <- merge(manatee_BUSCO, manatee_code, by.x=0, by.y="pipelineCode")
    manatee_BUSCO$Tool <- gsub("-", "\n", manatee_BUSCO$Tool)
    colnames(manatee_BUSCO) <- make.names(colnames(manatee_BUSCO))
    manatee_BUSCO <- manatee_BUSCO %>% filter(Row.names!="ONT1")
    
    total_busco <- 11366
    get_perc_busco_found <- function(x){
      r=(total_busco - as.numeric(x["Missing.BUSCOs"]))*100/total_busco 
      round(r, digits = 2)
    }
    
    get_perc_busco_complete <- function(x){
      c <- as.numeric(x["Complete.and.single.copy.BUSCOs"]) + as.numeric(x["Complete.and.duplicated.BUSCOs"]) 
      r=c*100/total_busco 
      round(r, digits = 2)
    }
    
    get_perc_busco_fragmented <- function(x){
      c <- as.numeric(x["Fragmented.BUSCOs"]) 
      r=c*100/total_busco 
      round(r, digits = 2)
    }
    
    get_perc_busco_missing <- function(x){
      c <- as.numeric(x["Missing.BUSCOs"]) 
      r=c*100/total_busco 
      round(r, digits = 2)
    }
    
    manatee_BUSCO$BUSCO_found <- apply(manatee_BUSCO,1,get_perc_busco_found)
    manatee_BUSCO$BUSCO_complete <- apply(manatee_BUSCO,1,get_perc_busco_complete)
    manatee_BUSCO$BUSCO_fragmented <- apply(manatee_BUSCO,1,get_perc_busco_fragmented)
    manatee_BUSCO$BUSCO_missing <- apply(manatee_BUSCO,1,get_perc_busco_missing)
    
    manatee_BUSCO$Tool <- factor(manatee_BUSCO$Tool, levels = c(sort(setdiff(unique(manatee_BUSCO$Tool), tool_name)), tool_name))
    
    # Modify `pD7.2` plot
    pD7.2 <- ggplot(manatee_BUSCO, aes(x=Label, y=as.numeric(BUSCO_found), color=Lib_Plat, shape=Data_Category)) + 
      geom_segment(aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_found), color=Lib_Plat), size=2) +
      geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
      pub_theme +
      theme_pubclean(flip=TRUE) +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=14),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,0)  # Remove margins around the plot
      ) +
      scale_fill_manual(values = libplat.palette, name="Library-Platform") +
      scale_color_manual(values = libplat.palette, name="Library-Platform") +
      scale_x_discrete(breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO", "CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                                  "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS", "cDNA-PacBio-LS"),
                       labels = c("SO", rep("LO", 6), rep("LS", 3))) +
      facet_grid(Tool ~., drop = TRUE, scales = "free_y") +
      xlab("") + 
      scale_y_continuous("", sec.axis = sec_axis(~ ., breaks = NULL, name = "Manatee"),
                         label = unit_format(unit = "%"), limits = c(-1, 80), expand = expansion(mult = c(0, 0.1))) +
      theme(
        strip.text.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.title = element_text(size=16)
      ) +
      theme(legend.position = "none") +
      coord_flip()
    
    # Calculate the appropriate height for the boxes in p.mid7 based on the number of levels
    num_tools <- length(unique(manatee_BUSCO$Tool))  # Number of unique tools
    manatee_BUSCO$Tool <- factor(manatee_BUSCO$Tool, levels = c(tool_name ,sort(setdiff(unique(manatee_BUSCO$Tool), tool_name), decreasing = TRUE)))
    
    # Modify middle plot (p.mid7) to have consistent height for boxes
    p.mid7 <- ggplot(manatee_BUSCO, aes(y = Tool)) + 
      geom_tile(aes(x = 1), fill = "gray", color = "white", width = 1.25, height = 0.95) +  # Adjust box height
      geom_text(aes(x = 1, label = Tool), size = 5) +                                   # Adjust text size
      theme_void() +                                                                    # Remove background and axes
      theme(
        plot.margin = margin(10,-60,25,0)  # Remove margins around the plot
      ) +
      scale_x_continuous(expand = c(0, 0)) +                                            # Remove space on the left
      coord_fixed(ratio = num_tools / 5)  # Adjust ratio to match the height of other plots
    
    # Combine the two plots
    gg2.7 <- ggplot_gtable(ggplot_build(pD7.2))
    g.mid <- ggplot_gtable(ggplot_build(p.mid7))
    
    # Arrange plots side by side with adjusted widths
    pD7 <- grid.arrange(g.mid, gg2.7, ncol = 2, widths = c(3/9, 6/9),
                        top = textGrob("BUSCO genes found (%)", gp = gpar(fontsize=18, font=1)))
    
    # Save the combined plot
    ggsave(file = paste0(outdir, "/panel4f1_manatee.svg"), plot = pD7, width = 8, height = 5)
    
  }
  append_new_data_4f2_manatee <- function(new_data_manatee){
    
    manatee_BUSCO$Tool <- factor(manatee_BUSCO$Tool, levels = c(sort(setdiff(unique(manatee_BUSCO$Tool), tool_name)), tool_name))
    
    pD7.2f <- ggplot(manatee_BUSCO, aes(x=Label, y=as.numeric(BUSCO_fragmented), color=Lib_Plat, shape=Data_Category)) + 
      geom_segment(aes(x=Label, xend=Label, y=0, yend=as.numeric(BUSCO_fragmented), color=Lib_Plat), size=2) +
      geom_point(position = position_dodge(width = 1), size=5, aes(fill=Lib_Plat)) +
      pub_theme +
      theme_pubclean(flip=TRUE) +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=14),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0,0,0,-20)  # Remove margins around the plot
      ) +
      scale_fill_manual(values = libplat.palette, name="Library-Platform") +
      scale_color_manual(values = libplat.palette, name="Library-Platform") +
      scale_x_discrete(breaks = c("cDNA-Illumina-SO", "CapTrap-ONT-LO", "R2C2-ONT-LO", "cDNA-ONT-LO", "CapTrap-PacBio-LO", "cDNA-PacBio-LO",
                                  "dRNA-ONT-LO", "dRNA-ONT-LS", "cDNA-ONT-LS", "cDNA-PacBio-LS"),
                       labels = c("SO", rep("LO", 6), rep("LS", 3))) +
      facet_grid(Tool ~., drop = TRUE, scales = "free_y") +
      xlab("") + 
      scale_y_continuous("", sec.axis = sec_axis(~ ., breaks = NULL, name = "Manatee"),
                         label = unit_format(unit = "%"), limits = c(-1, 20), expand = expansion(mult = c(0, 0.1))) +
      theme(
        strip.text.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.title = element_text(size=16)
      ) +
      theme(legend.position = "none") +
      coord_flip()
    
    gg2.7f <- ggplot_gtable(ggplot_build(pD7.2f))
    
    pD7f <- grid.arrange(g.mid, gg2.7f,ncol=2,widths=c(3/9, 6/9),
                         top = textGrob("BUSCO genes fragmented (%)",gp=gpar(fontsize=18,font=1)))
    
    ggsave(file=paste0(outdir, "/panel4f2_manatee.svg"), plot=pD7f, width=8, height=5)
    
    
  }
}


### Script for figure 4g1
print('Script for figure 4g1 does not yet exist...')


### Script for figure 4g2

append_new_data_4g2 <- function(new_data){
  
  # ///
  manatee_SIRV_metrics <- rbind(manatee_SIRV_metrics, new_data)
  rownames(manatee_SIRV_metrics)[nrow(manatee_SIRV_metrics)] <- pipelineCode
  # \\\
  
  manatee_code <- read.csv("Challenge3_Figures_Data/manatee/manatee_code.txt", sep=",", header = T )
  manatee_code <- rbind(manatee_code, meta_data)
  manatee_code$Lib_Plat <- apply(manatee_code, 1, function(x){
    paste(x["Library_Preps"], x["Platform"], sep = "-")
  })
  manatee_code$Lib_DC=apply(cbind(manatee_code[,c("Library_Preps", "Data_Category")]), 1, paste, collapse="-")
  manatee_code$Label <-apply(cbind(manatee_code[,c("Library_Preps", "Platform", "Data_Category")]), 1, paste, collapse="-")
  
  manatee_SIRV_metrics <- merge(manatee_SIRV_metrics, manatee_code, by.x=0, by.y="pipelineCode")
  manatee_SIRV_metrics[,'Redundancy'] <- as.numeric(manatee_SIRV_metrics[,'Redundancy'])
  manatee_SIRV_metrics$"1/Red" <- round(1/manatee_SIRV_metrics$Redundancy,2)
  manatee_SIRV_metrics$Label <-apply(cbind(manatee_SIRV_metrics[,c("Library_Preps", "Platform","Data_Category")]), 1, paste, collapse="_")
  
  mana_metrics2 <- data.frame (as.numeric(manatee_SIRV_metrics$Sensitivity),
                               as.numeric(manatee_SIRV_metrics$`Positive Detection Rate`),
                               as.numeric(manatee_SIRV_metrics$Precision),
                               as.numeric(manatee_SIRV_metrics$`Non Redundant Precision`),
                               as.numeric(manatee_SIRV_metrics$`False Discovery Rate`),
                               as.numeric(manatee_SIRV_metrics$`1/Red`),
                               manatee_SIRV_metrics$Tool,
                               manatee_SIRV_metrics$Label
  )
  
  colnames(mana_metrics2) <- c("Sen", "PDR", "Pre", "nrPre", "FDR", "1/Red", "Tool", "Label")
  mycolors = c( "cDNA_PacBio_LO"="#d8527c", "cDNA_PacBio_LS"="#9a133d", 
                "cDNA_ONT_LO"="#6996e3", "cDNA_ONT_LS"="#1a318b" )
  
  # remove rnaSpades as it has No results
  
  mana_metrics3 <- mana_metrics2[mana_metrics2$Tool != "rnaSPAdes",]
  
  tool = c(sort(setdiff(unique(mana_metrics3$Tool), tool_name)), tool_name)
  plat_lib <- unique(mana_metrics3$Label)
  
  pdf(paste0(outdir, "/panel4g2.pdf"))
  par(mar=c(0,0.5,1.5,0) + 0.1)
  #layout(matrix(c(1:12,13,13,13,13),  byrow = TRUE, ncol = 4,nrow = 4)) # for main figure
  layout(matrix(c(1:7,8),  byrow = TRUE, ncol = 2,nrow = 4))
  
  for ( i in 1:length(tool)) {
    sel <- mana_metrics3$Tool == tool[i]
    subms <- mana_metrics3[sel,c(1:6)]
    subms <- rbind(rep(1,ncol(subms)) , rep(0,ncol(subms)), subms)
    radarchart(subms, axistype=1 , title= tool[i], 
               #custom polygon
               pcol= mycolors[mana_metrics3[sel,"Label"]], plwd=2 , plty=1.3,
               #custom the grid
               cglcol = "grey", cglty=1, axislabcol="grey", cglwd = 0.8, caxislabels=seq(0,0.25,1), 
               #custom labels
               #vlabels = rep("",6), palcex = 0, cex.main=1.5, paxislabels = 0
               vlabels = colnames(subms) , palcex = 1, cex.main=1.2, paxislabels = 0, col = c(2,2,2,3,3,1)
    )
  }
  
  # add an empty plot
  plot(x = 1, y = 1, col = "white", xaxt='n', yaxt='n', bty ="n",axes=F,frame.plot=F)
  
  # Add a legend
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  plot_colors <- c("blue","black", "green", "orange", "pink")
  l <- legend(x = "top",inset = 0, legend = names(mycolors), horiz = TRUE, bty = "n", pch=20 , col=mycolors , text.col = "black", cex=0.8, pt.cex=2)
  dev.off()
  
} 

print('All Existing Scripts are Sucessfully Loaded!')

