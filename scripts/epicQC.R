############### packages installation ###############

# packages installation
# CRAN packages installation
cran_packages <- c("ggplot2", "data.table", "tidyr", "dplyr", "FNN", "readxl", "tibble", 
                   "devtools", "plotly", "stats", "qqman", "gaston", "tidyverse")

install.packages(cran_packages)

# load CRAN packages
lapply(cran_packages, library, character.only = TRUE)

# Bioconductor packages installation
bioconductor_packages <- c("Biobase", "sva", "minfi", "limma", "sesame", "ENmix", "ChAMP", 
                           "ChAMPdata", "FlowSorted.Blood.EPIC", "wateRmelon", "meffil","FDb.InfiniumMethylation.hg19")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(bioconductor_packages)

# Charge Bioconductor package
lapply(bioconductor_packages, library, character.only = TRUE)

# github packages
install_github("perishky/meffil")
library(meffil)

# devtools package
devtools::install_version("matrixStats", version = "1.1.0")
library(readr)
library(matrixStats)



############### parameters initialization ###############

# Present a menu to the user to choose arraytype
arraytype <- menu(c("450K", "EPICv1", "EPICv2"), title = "Choose your array type (450k are not supported):")

if (arraytype == 1) {
  arraytype <- "450K"
} else if (arraytype == 2) {
  arraytype <- "EPIC"
  arraytypeV <- "EPICv1"
} else if (arraytype == 3) {
  arraytype <- "EPIC"
  arraytypeV <- "EPICv2"
}

# Output the user's choice
cat("You selected option", arraytype, "\n")

# load manifest library according to user's choice
if (arraytypeV == 'EPICv1') {
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(IlluminaHumanMethylationEPICmanifest)
  data(EPIC.manifest.hg19)
} else if (arraytypeV == 'EPICv2') {
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  library(IlluminaHumanMethylationEPICv2manifest)
  data(EPIC.manifest.hg38)
}

# directories
# set the working and output directories
basedir <- ""
workingdir <- ""

# Set the working directory
setwd(basedir)

# and create outpout directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

############### Import and process DNA methylation data ###############

# Import IDAT files using minfi
# List all IDAT files in the specified directory
idat_files <- list.files(path = basedir, pattern = "\\.idat$", recursive = TRUE)
num_idat_files <- length(idat_files) ]
num_idat_files

# read the IDAT files into an RGChannelSet extended with the minfi package. Gives general infos, and the array type
RGsetbatch <- read.metharray.exp(basedir,  recursive = TRUE, extended = TRUE)
RGsetbatch

# load the samplesheet
Sample_sheet <- read_csv('./207705780011/SampleSheet2.csv')

# Extract raw beta values from the RGChannelSet
beta_raw <- getBeta(preprocessRaw(RGsetbatch))

# 2. Basic Quality Control
# Perform basic QC to check chip-wide medians of the Meth and Unmeth channels
# creation of a MethylSet object without normalization, to use minfi functions 
MethylSet <- preprocessRaw(RGsetbatch)

# chipwide medians of the Meth and Unmeth channels.
qc <- getQC(MethylSet)
plotQC(qc, badSampleCutoff = 10.5)

# minfi qc 
minfiQC(MethylSet)

# Check bisulfite conversion rates
# Extract bisulfite conversion percentage for each array, with watermelon package
# for 450k data 
bisConv <- as.data.frame(bscon(RGsetbatch))
bisConv$Basename <- rownames(bisConv)
colnames(bisConv) <- c("Bisulfite_conversion_%", "Basename")

# Identify samples that failed bisulfite conversion (below 80%)
failed_bisConv <- bisConv[bisConv$`Bisulfite_conversion_%` < 80,]

# Remove samples that failed bisulfite conversion and clean the dataset
Sample_sheet_clean <- Sample_sheet[!Sample_sheet$Sample_Name %in% failed_bisConv$Basename,]
RGsetbatch_Sclean <- RGsetbatch[, Sample_sheet_clean$Sample_Name]
beta_raw_Sclean <- getBeta(RGsetbatch_Sclean)




############### Noob normalization ###############

# Normalize the data using the Noob method
minfi_import_Noob <- preprocessNoob(RGsetbatch_Sclean)
beta_minfi_import_Noob <- getBeta(minfi_import_Noob)

# Normalize the data using the Noob method
minfi_import_Noob <- preprocessNoob(RGsetbatch_Sclean)
Mvalues_minfi_import_Noob <- getM(minfi_import_Noob)
beta_minfi_import_Noob <- getBeta(minfi_import_Noob)

# save the .Rdata image of the normalized data
save(beta_minfi_import_Noob,file = file.path(output_dir,"beta_minfi_import_Noob.Rdata"), compress= F)

# save the .csv file of the normalized data
write.csv(beta_minfi_import_Noob, file = file.path(output_dir,"beta_minfi_import_Noob.csv"))


############### filtering and normalization with ChAMP ###############

# Import data with ChAMP
# sample sheet must be in the same directory as the IDATs
import <- champ.import(directory = basedir,
                       offset = 100,
                       arraytype = arraytype)

import_flt <- champ.filter(beta = import$beta,
                           M = NULL,
                           pd = Sample_sheet_clean,
                           intensity = NULL,
                           Meth = NULL,
                           UnMeth = NULL,
                           detP = import$detP,
                           beadcount = import$beadcount,
                           autoimpute = TRUE,
                           filterDetP = TRUE,
                           ProbeCutoff = 0.5,
                           SampleCutoff = 0.1,
                           detPcut = 0.01,
                           filterBeads = TRUE,
                           beadCutoff = 0.05,
                           filterNoCG = TRUE,
                           filterSNPs = FALSE,
                           population = NULL,
                           filterMultiHit = TRUE,
                           filterXY = TRUE,
                           fixOutlier = TRUE,
                           arraytype = arraytype)

champ.QC(beta = import$beta,
         pheno = Sample_sheet_clean$Sample_Name,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir=output_dir)


# Generate a beta matrix with SNPs for ancestry imputation


# Extract probes marked as SNPs
SNPs_probes <- rownames(mani)[which(mani$MASK_general == TRUE)]
# Filter SNP probes to include only those present in the filtered beta matrix
valid_SNPs_probes <- SNPs_probes[SNPs_probes %in% rownames(import_flt$beta)]
# Extract beta values for the valid SNP probes
SNPs_beta <- import_flt$beta[valid_SNPs_probes, ]


# Remove SNP probes from the filtered beta matrix
import_flt$beta <- import_flt$beta[!rownames(import_flt$beta) %in% SNPs_probes,]

# Select the probes that passed CHAMP filtering in the Noob normalized beta matrix
beta_minfi_import_Noob_filtered <- beta_minfi_import_Noob[rownames(import_flt$beta),]

# 6. Perform BMIQ normalization
# Normalize the beta values using the BMIQ method
beta_minfi_import_Noob_filtered_BMIQ <- champ.norm(beta = beta_minfi_import_Noob_filtered,
                                                   resultsDir = output_dir,
                                                   method = "BMIQ",
                                                   plotBMIQ = TRUE,
                                                   arraytype = arraytype,
                                                   cores = 3)



############### Plot QC metrics and distributions ###############
# Plot and save QC plots

#Define a function to save plots (the ones that are not ggplot)
save_plot <- function(filename, plot_code) {
  png(filename, width = 2400, height = 1800, bg = "white", res=300)
  par(bg = "transparent")
  plot_code
  dev.off()
}


save_plot(file.path(output_dir, "qc_plot.png"), {
  plotQC(qc)
})

save_plot(file.path(output_dir, "raw_beta_distribution.png"), {
  densityPlot(getBeta(preprocessRaw(RGsetbatch)), main = "Raw beta distribution")
})

save_plot(file.path(output_dir, "bisulfite_conversion_percentage.png"), {
  plot(bisConv$`Bisulfite_conversion_%`, ylab = "Bisulfite conversion (%)", xlab = "Index", main = "Bisulfite Conversion Percentage", type = "p")
})

save_plot(file.path(output_dir, "raw_beta_distribution_qc_done.png"), {
  densityPlot(RGsetbatch_Sclean, main = "Raw beta distribution (QC done)")
})

save_plot(file.path(output_dir, "raw_m_values_qc_done.png"), {
  Mvalues_Sclean <- B2M(beta_raw_Sclean)
  densityPlot(Mvalues_Sclean, main = "Raw M values (QC done)")
})

save_plot(file.path(output_dir, "noob_beta_distribution_qc_done.png"), {
  densityPlot(beta_minfi_import_Noob, main = "Noob beta distribution (QC done)")
})

save_plot(file.path(output_dir, "noob_filtered_beta_distribution_qc_done.png"), {
  densityPlot(beta_minfi_import_Noob_filtered, main = "Noob filtered beta distribution (QC done)")
})

save_plot(file.path(output_dir, "noob_bmiq_filtered_beta_distribution_qc_done.png"), {
  densityPlot(beta_minfi_import_Noob_filtered_BMIQ, main = "Noob+BMIQ filtered beta distribution (QC done)")
})

save_plot(file.path(output_dir, "noob_bmiq_filtered_m_values_distribution_qc_done.png"), {
  Mvalues_Sclean_norm_fileterd <- B2M(beta_minfi_import_Noob_filtered_BMIQ)
  densityPlot(Mvalues_Sclean_norm_fileterd, main = "Noob+BMIQ filtered M values distribution (QC done)")
})

save_plot(file.path(output_dir, "snps_beta_distribution.png"), {
  densityPlot(SNPs_beta, main = "SNP probes beta distribution")
})
# Note: Consider filtering out any samples with abnormal distributions based on the plots


############### sex check ###############

# Map to genome and estimate sex with minfi
GMsetEx <- mapToGenome(RGsetbatch_Sclean)
estSex <- getSex(GMsetEx)

# Estimate sex with Watermelon
sexWater<- estimateSex(beta_minfi_import_Noob, do_plot = TRUE)

# Combine Minfi and Watermelon predictions
# Convert both data frames to standard data frames if they are not
estSex_df <- as.data.frame(estSex)
sexWater_df <- as.data.frame(sexWater)

# Add rownames as a new column to use as a merge key
# Merge the data frames based on row names
merged_df <- merge(estSex_df, sexWater_df, by = "row.names", suffixes = c("_est", "_water"))

# Rename the "Row.names" column to "SampleID" for clarity
colnames(merged_df)[1] <- "SampleID"

# Normalize the predicted sex columns
normalize_sex <- function(sex) {
  # Map "F" to "Female", "M" to "Male", and leave other values as is
  sex <- tolower(sex)  # Convert to lowercase
  if (sex == "f") return("female")
  if (sex == "m") return("male")
  return(sex)  # Keep other labels like "45,XO" unchanged
}

merged_df$normalized_est <- sapply(merged_df$predictedSex, normalize_sex)
merged_df$normalized_water <- sapply(merged_df$predicted_sex, normalize_sex)

# Check for mismatches
mismatch <- merged_df$normalized_est != merged_df$normalized_water

if (any(mismatch)) {
  print("Sex mismatch detected for the following samples:")
  print(merged_df[mismatch, c("SampleID", "predictedSex", "predicted_sex")])
} else {
  cat("No mismatches in predicted sex.\n")
}

# Save the data frame to the specified file path
write.csv(water_minfi_unmatch, file.path(output_dir, "water_minfi_unmatch_sex.csv"), row.names = FALSE)






