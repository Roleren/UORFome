
#' Bazzini 2012 Experiment info table
#' LIBRARY TYPES: RFP, RNA
#' STAGES: 64, sphere, shield (2, 4, 6)
getBazzini12 <- function() {
  return(read.experiment("/export/valenfs/data/processed_data/experiment_tables_for_R/bazzini2012.csv"))
}

#' RCP-seq standard Experiment info table
#' LIBRARY TYPES: SSU, LSU, RNA
#' STAGES: 64, shield
getRCPMZ <- function(){
  return(read.experiment("/export/valenfs/data/processed_data/experiment_tables_for_R/rcp-experiment_WTvsMZ.csv"))
}
#' RCP-seq standard Experiment info table
#' LIBRARY TYPES: fractions (9-20 (-11))
#' shield
getRCP5 <- function(){
  return(read.experiment("/export/valenfs/data/processed_data/experiment_tables_for_R/Valen5_2018.csv"))
}
#' Archer 2016 Experiment info table
#' LIBRARY TYPES: SSU, RFP
#' STAGES:
getArcher16 <- function(){
  return(read.experiment("/export/valenfs/data/processed_data/experiment_tables_for_R/archer16_TCP.csv"))
}

#' Archer 2019 Experiment info table
#' LIBRARY TYPES: SSU, RNA
#' STAGES:
getArcher19 <- function(){
  return(read.experiment("/export/valenfs/data/processed_data/experiment_tables_for_R/archer19_TCP.csv"))
}

#' Archer 2019 Experiment info table
#' LIBRARY TYPES: SSU, RNA
#' STAGES: ?
getArcher19Bed <- function(){
  return(read.experiment("/export/valenfs/data/processed_data/experiment_tables_for_R/archer19_TCP_BED.csv"))
}
#' McManus 2017 Experiment info table
#' LIBRARY TYPES: CAGE
getMcManus17 <- function(){
  return(read.experiment("/export/valenfs/data/processed_data/experiment_tables_for_R/McManus2017_CAGE.csv"))
}
#' Wery 2015 Experiment info table
#' LIBRARY TYPES: CAGE
getWery15<- function(){
  return(read.experiment("/export/valenfs/data/processed_data/experiment_tables_for_R/Wery2015.csv"))
}

getRFPChew13 <- function() {
  stop("Not ready")
  return(read.experiment("/export/valenfs/data/processed_data/experiment_tables_for_R/Wery2015.csv"))
}

#' RCP-seq SSU Experiment info table
#' LIBRARY TYPES: SSU, RNA
#' STAGES: 64, sphere, shield (2, 4, 6)
getTCPdf <- function(stage = c("64", "sphere", "shield"), 
                     type = c("WT", "WT", "WT"),
                     SSU = c("/export/valenfs/projects/adam/TCP_seq/RCP_files/64cell_SSU_reps_1_2_peaks_removed_translating_filter.bam",
                             "/export/valenfs/projects/adam/TCP_seq/RCP_files/sphere_SSU_reps_1_2_3_peaks_removed_translating_filter.bam",
                             "/export/valenfs/projects/adam/TCP_seq/RCP_files/shield_SSU_reps_1_2_3_peaks_removed_translating_filter.bam"),
                     RFP = c("/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/256Cell_trimmed.bam",
                             "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Dome_trimmed.bam",
                             "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Shield_trimmed.bam")) {
  
  return(data.frame(SSU, RFP, stage, type, stringsAsFactors = FALSE))
}
#' RCP-seq SSU Experiment info table
#' LIBRARY TYPES: SSU, LSU, RFP, RNA
#' STAGES: 64, sphere, shield (2, 4, 6)
getTCPdfAll <- function(stage = c("64", "sphere", "shield"), 
                        type = c("WT", "WT", "WT"),
                        SSU = c("/export/valenfs/projects/adam/TCP_seq/RCP_files/64cell_SSU_reps_1_2_peaks_removed_translating_filter.bam",
                                "/export/valenfs/projects/adam/TCP_seq/RCP_files/sphere_SSU_reps_1_2_3_peaks_removed_translating_filter.bam",
                                "/export/valenfs/projects/adam/TCP_seq/RCP_files/shield_SSU_reps_1_2_3_peaks_removed_translating_filter.bam"),
                        RFP = c("/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/256Cell_trimmed.bam",
                                "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Dome_trimmed.bam",
                                "/export/valenfs/data/processed_data/Ribo-seq/chew_2013_zebrafish/final_results/aligned_GRCz10/Shield_trimmed.bam"),
                        RNA = c("/export/valenfs/data/processed_data/RNA-seq/lee_2013_zebrafish/total_RNA/aligned_GRCz10/WT_2hpf_Tota_mRNA.bam",
                                "/export/valenfs/data/processed_data/RNA-seq/lee_2013_zebrafish/total_RNA/aligned_GRCz10/WT_4hpf_Total_mRNA.bam",
                                "/export/valenfs/data/processed_data/RNA-seq/lee_2013_zebrafish/total_RNA/aligned_GRCz10/WT_6hpf_Total_mRNA_merged.bam"),
                        LSU = c("/export/valenfs/projects/adam/TCP_seq/RCP_files/64cell_LSU_reps_1_2_peaks_removed_translating_filter.bam",
                                "/export/valenfs/projects/adam/TCP_seq/RCP_files/sphere_LSU_reps_1_2_3_peaks_removed_translating_filter.bam",
                                "/export/valenfs/projects/adam/TCP_seq/RCP_files/shield_LSU_reps_1_2_3_peaks_removed_translating_filter.bam")) {
  
  return(data.frame(SSU, RFP, RNA, LSU, stage, type, stringsAsFactors = FALSE))
}

#' RCP-seq SSU 2019 mapped Experiment info table
#' LIBRARY TYPES: LSU
#' STAGES: 64, sphere, shield (2, 4, 6)
getTCPNew <- function(){
  mergedF <- "/export/valenfs/data/processed_data/TCP-seq/valen_all_withrRNA/aligned/"
  p <- paste0
  df2 <- data.frame(LSU = c(p(mergedF, "64_cell_LSU_V7.bam"), p(mergedF, "64_cell_LSU_V8.bam"),
                            p(mergedF, "shield_V5_merged_LSU.bam"), p(mergedF, "shield_V6_merged_LSU.bam"),
                            p(mergedF, "shield_V15_merged_LSU.bam"), p(mergedF, "shield_all_merged_LSU.bam"),
                            p(mergedF, "sphere1_V7_merged_LSU.bam"), p(mergedF, "sphere2_V7_merged_LSU.bam"),
                            p(mergedF, "sphere3_V7_merged_LSU.bam"), p(mergedF, "64_LSU_V12_4Ei.bam")),
                    stage = c(rep("64", 2), rep("shield", 4), rep("sphere", 3), "64"),
                    type = c(rep("WT", 2), rep("WT", 3), "WT_ALL", rep("WT", 3), "4Ei"), stringsAsFactors = FALSE)
  return(df2)
}