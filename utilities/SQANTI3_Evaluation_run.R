#####################################
##### SQANTI3 report generation ######
#####################################



### Francisco J Pardo-Palacios
### Last Modified: 05/03/2020 by francisco.pardo.palacios@gmail.com

#********************** Taking arguments from python script

args <- commandArgs(trailingOnly = TRUE)
args <- c( '~/Downloads/lrgasp-challenge-3_benchmarking_workflow-main/lrgasp-challenge-3_full_data/output_Fabian/transcriptome_classification.txt',
           '~/Downloads/lrgasp-challenge-3_benchmarking_workflow-main/lrgasp-challenge-3_full_data/output_Fabian/transcriptome_junctions.txt',
           'output_fabian',
           '/Users/woutermaessen/PycharmProjects/lrgasp-challenge-3_benchmarking_docker/lrgasp_event2_metrics/utilities',
           'PacBio',
           '/Users/woutermaessen/Downloads/lrgasp-challenge-3_benchmarking_workflow-main/lrgasp-challenge-3_full_data/output_fabian/rna_Rdata/',
           '/Users/woutermaessen/Downloads/lrgasp-challenge-3_benchmarking_workflow-main/lrgasp-challenge-3_full_data/output_Fabian/BUSCO_results.tsv', 
           'manatee',
           'IsONform',
           'IONf',
           'LS',
           'cDNA',
           'fabian_lab',
           'IONf_cDNA'
)

class.file <- args[1]
junc.file <- args[2]
name <- args[3]
utilities.path <- args[4]
platform <- args[5]
rdata <- args[6]
busco <- args[7]
organism <- args[8]
tool_name <- args[9]
alias <- args[10]
data_category <- args[11]
library_prep <- args[12]
lab <- args[13]
pipelineCode <- args[14] # Script should be able to create pipelineCode itself

###TO DO
# Create pipelineCode self instead of with args
# check if script 4f2 for mouse is working -> it did not work first for manatee
# Create functions for missing figures

meta_data <- data.frame(pipelineCode = pipelineCode, Library_Preps = library_prep, Platform = platform, Data_Category = data_category, Lab = lab, Tool = tool_name, Alias = alias)


report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1];
report.file <- paste(report.prefix, "Evaluation_report.html", sep="_");
bam_file <- paste(report.prefix, "corrected.bam", sep="_")


#********************** Packages (install if not found)

list_of_packages <- c("ggplot2", "scales", "knitr","rmarkdown")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

library(ggplot2)
library(scales)
library(knitr)
library(rmarkdown)
library(Rsamtools)

#********************* Run Calculation scripts

setwd(utilities.path)
source("LRGASP_calculations.R")

LRGASP_calculations_challenge3(NAME = name , out.dir = rdata,
                      class.file=class.file, junc.file=junc.file,
                      platform = platform, 
                      functions.dir = utilities.path,
                      bam = bam_file,
                      organism = organism)


busco_table = read.table(busco, sep="\t", header=F)
busco_table = as.data.frame(t(busco_table))
busco_results = data.frame(row.names = busco_table[,1])
busco_results[,"Absolute value"]=apply(busco_table, 1, function(Y) as.integer(Y[2]))
total_BUSCO = sum(busco_results[,"Absolute value"])
busco_results[,"Relative value (%)"] = apply(busco_results,1, function(Z){
  round( ((Z[1]/total_BUSCO)*100), digits = 2)
})

save(busco_results , file = paste(name, "_BUSCO.RData", sep = ''))

source('/Users/woutermaessen/Desktop/ConesaInternship/ScriptsToRecreateFiguresLRGASP_Manatee/LRGASP_Benchmark_Report.R')

LRGASP_Benchmark_Report(NAME = name , out.dir = rdata,
                        class.file=class.file, junc.file=junc.file,
                        platform = platform, 
                        functions.dir = utilities.path,
                        bam = bam_file,
                        organism = organism,
                        tool_name = tool_name,
                        alias = alias,
                        data_category = data_category,
                        library_prep = library_prep,
                        lab = lab
)

RMD = paste(utilities.path, "Evaluation_metrics.Rmd", sep = "/")
RMD = '/Users/woutermaessen/Desktop/ConesaInternship/ScriptsToRecreateFiguresLRGASP_Manatee/Evaluation_test.Rmd'

rmarkdown::render(RMD, params = list(
  output.directory = rdata,
  Name = name,
  Platform = platform ), output_file = report.file
)

