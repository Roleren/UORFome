library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(cowplot)
library(ORFikPipeline)
library(gridExtra)
library(reshape2)
library(cowplot)
Palette1 <- c('skyblue4', 'orange')

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# START HERE (Run fivePrimePeaks_create.R if they do not exist!)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Pick shield or 64 here
# 64
res64      <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/TISU_matrix_new2020feb_64cell.rds")
res64TSS  <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/TISU_matrix_new2020feb_tss_64cell.rds")
# shield
resShi     <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/TISU_matrix_new2020feb_shield.rds")
resShiTSS <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/TISU_matrix_new2020feb_tss_shield.rds")
# With pseudo counts
res64p <- copy(res64)#; res64p$SE <- res64p$SE + 0.1 # pseudo count of 0.0001
resShip <- copy(resShi)#; resShip$SE <- resShip$SE + 0.1 # pseudo count of 0.0001
res64TSSp <- copy(res64TSS)#; res64TSSp$SE <- res64TSSp$SE + 0.1 # pseudo count of 0.0001
resShiTSSp <- copy(resShiTSS)#; resShiTSSp$SE <- resShiTSSp$SE + 0.1 # pseudo count of 0.0001
# count ratio
res <- readRDS("/export/valenfs/projects/Hakon/RCP_SEQ/matrix_new2020feb_tss_countsratio.rds")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# PLOTS
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

# Plots for suplements:

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# With pseudo counts
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#

# 4Ei vs WT (64)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res64p[stage == "64cell",]), aes(x = log2(SE), colour= as.factor(condition))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Condition (WT vs 4Ei) (64)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# 4Ei vs WT (shield)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(resShip[stage == "shield",]), aes(x = log2(SE), colour= as.factor(condition))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Condition (WT vs 4Ei) (shield)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()

# TISU vs background (64)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res64p[stage == "64cell",]), aes(x = log2(SE), colour= as.factor(paste(condition, type)))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Condition & TISU status (WT vs 4Ei) (64)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# TISU vs background (shield)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(resShip[stage == "shield",]), aes(x = log2(SE), colour= as.factor(paste(condition, type)))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Condition & TISU status (WT vs 4Ei) (shield)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()

# Box plot
# 64
summary(res64[stage == "64cell" & condition == "WT",]$SE)
summary(res64[stage == "64cell" & condition == "4Ei",]$SE)
summary(res64p[stage == "64cell" & condition == "WT",]$SE)
summary(res64p[stage == "64cell" & condition == "4Ei",]$SE)

ggplot(res64[stage == "64cell",], aes(y = log2(SE), x = condition, fill = condition)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("64 cell, with pseudo")
ggslackR()
# shield
ggplot(resShip[stage == "shield",], aes(y = log2(SE), x = condition, fill = condition)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("shield, with pseudo")
# TISU
# 64
ggplot(res64p[stage == "64cell" & type == "TISU",], aes(y = log2(SE), x = condition, fill = condition)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("TISU 64 cell, with pseudo")
#shield
ggplot(resShip[stage == "shield" & type == "TISU",], aes(y = log2(SE), x = condition, fill = condition)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("TISU shield, with pseudo")

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TISU boxplots (USING TSS filter)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# 64
res64andShiTisu <- rbind(res64TSSp[stage == "64cell" & type == "TISU",], resShiTSSp[stage == "shield" & type == "TISU",])
res64andShiTisu[stage == "64cell", stage := "64-cell"]
res64andShiTisu[stage == "shield", stage := "Shield"]
res64TSSp[stage == "64cell" & type == "TISU",][, summary(SE), by = condition]
gg64 <- ggplot(res64TSSp[stage == "64cell" & type == "TISU",], aes(y = log2(SE), x = factor(condition, levels = c("WT", "4Ei")), fill = factor(condition, levels = c("WT", "4Ei")))) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("TSS TISU 64 cell, with pseudo") +
  ylab("log2(SE)") +
  xlab("")+
  labs(fill = "Condition")
#theme(legend.position = "none")
gg64
ggslackR()
# shield
ggShi <- ggplot(resShiTSSp[stage == "shield" & type == "TISU",], aes(y = log2(SE), x = factor(condition, levels = c("WT", "4Ei")), fill = factor(condition, levels = c("WT", "4Ei")))) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("TSS TISU shield, with pseudo") +
  ylab("") +
  xlab("") +
  labs(fill = "Condition")
ggShi
ggslackR()

comb <- grid.arrange(plot_grid(gg64, ggShi, nrow = 1, ncol = 2, align = 'hv', axis  = 'tb'))
ggslackR(plot = comb, width = 200, height = 200)
# Merged plot
gg64Shi <- ggplot(res64andShiTisu, aes(y = log2(SE), x = factor(condition, levels = c("WT", "4Ei"), labels = c("WT", "10 uM 4Ei-10")), fill = factor(condition, levels = c("WT", "4Ei"), labels = c("WT", "10 uM 4Ei-10")))) +
  geom_boxplot() +
  theme_classic(base_size = 14) +
  ggtitle("") +
  ylab("log2 (TSS SSU peak / RNA)") +
  xlab("") +
  labs(fill = "Treatment") +
  facet_wrap( ~ stage) +
  scale_fill_manual(values = c("#504e52", "#53ceff")) +
  theme(
    # Change legend key size and key width
    legend.key.size = unit(1.5, "cm"),
    legend.key.width = unit(1.5,"cm"),
    strip.background = element_rect(fill="lightgray")
  )
gg64Shi
ggslackR()
summary(res64andShiTisu[variable == "64_SE"]$SE)
summary(res64andShiTisu[variable == "64_4Ei_SE"]$SE)
summary(res64andShiTisu[variable == "shield_SE"]$SE)
summary(res64andShiTisu[variable == "shield_4Ei_SE"]$SE)

wilcox.test(res64andShiTisu[variable == "64_SE"]$SE, res64andShiTisu[variable == "64_4Ei_SE"]$SE)
wilcox.test(res64andShiTisu[variable == "shield_SE"]$SE, res64andShiTisu[variable == "shield_4Ei_SE"]$SE)
# Median percentage change
median((1 - ((res64andShiTisu[condition != "WT" & stage == "64-cell",]$SE / res64andShiTisu[condition == "WT" & stage == "64-cell",]$SE)))*100)
median((1 - ((res64andShiTisu[condition != "WT" & stage == "Shield",]$SE / res64andShiTisu[condition == "WT" & stage == "Shield",]$SE)))*100)
mean((1 - ((res64andShiTisu[condition != "WT" & stage == "64-cell",]$SE / res64andShiTisu[condition == "WT" & stage == "64-cell",]$SE)))*100)
mean((1 - ((res64andShiTisu[condition != "WT" & stage == "Shield",]$SE / res64andShiTisu[condition == "WT" & stage == "Shield",]$SE)))*100)
table(res64andShiTisu[condition == "WT" & stage == "64-cell",]$SE - res64andShiTisu[condition != "WT" & stage == "64-cell",]$SE  >= 0)
table(res64andShiTisu[condition == "WT" & stage == "Shield",]$SE - res64andShiTisu[condition != "WT" & stage == "Shield",]$SE  >= 0)

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# TSS ratio count box plot for figure (S4-B)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
lvls <- c("64_WT_counts", "64_4Ei_0.1_counts", "4Ei_10_counts")
gg64 <- ggplot(res, aes(y = log2(counts), x = factor(variable, levels = lvls), fill = factor(variable, levels = lvls))) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("") +
  ylab("Log2(TSS SSU counts / downstream SSU counts") +
  xlab("")+
  labs(fill = "Treatment")
#theme(legend.position = "none")
gg64

gg64 <- ggplot(res, aes(y = counts, x = factor(variable, levels = lvls, labels = c("WT", "0.1 uM 4Ei-10", "10 uM 4Ei-10")), fill = factor(variable, levels = lvls, labels = c("WT", "0.1 uM 4Ei-10", "10 uM 4Ei-10")))) +
  geom_boxplot() +
  theme_classic(base_size = 14) +
  ggtitle("") +
  ylab("TSS SSU counts / downstream (2-100 nt) SSU counts") +
  xlab("") +
  labs(fill = "Treatment") +
  coord_cartesian(ylim = c(0, 0.4)) +
  scale_fill_manual(values = c("#504e52", "#418ec4", "#53ceff")) +
  theme(
    # Change legend key size and key width
    legend.key.size = unit(1.5, "cm"),
    legend.key.width = unit(1.5,"cm")
  )

#theme(legend.position = "none")
gg64
ggslackR()
res2 <- copy(res)
wilcox.test(res[variable == "64_WT_counts"]$counts, res[variable == "64_4Ei_0.1_counts"]$counts)
wilcox.test(res[variable == "64_WT_counts"]$counts, res[variable == "4Ei_10_counts"]$counts)
wilcox.test(res[variable == "64_4Ei_0.1_counts"]$counts, res[variable == "4Ei_10_counts"]$counts)
summary(res[variable == "64_WT_counts"]$counts)
summary(res[variable == "64_4Ei_0.1_counts"]$counts)
summary(res[variable == "4Ei_10_counts"]$counts)

median((1 - ((res[variable == "64_4Ei_0.1_counts",]$counts / res[variable == "64_WT_counts",]$counts)))*100, na.rm = T)
median((1 - ((res[variable == "4Ei_10_counts",]$counts / res[variable == "64_WT_counts",]$counts)))*100, na.rm = T)
mean((1 - ((res[variable == "64_4Ei_0.1_counts",]$counts / res[variable == "64_WT_counts",]$counts)))*100, na.rm = T)
mean((1 - ((res[variable == "4Ei_10_counts",]$counts / res[variable == "64_WT_counts",]$counts)))*100, na.rm = T)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Sanity test plots (NOT USED!!!!!!!!!!!!!!)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# Without pseudo counts (OLD FROM ARTICLE)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
# 4Ei vs WT (64)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res[stage == "64cell",]), aes(x = log2(SE), colour= as.factor(condition))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Condition (WT vs 4Ei) (64)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# 4Ei vs WT (shield)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res[stage == "shield",]), aes(x = log2(SE), colour= as.factor(condition))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Condition (WT vs 4Ei) (shield)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# TISU vs background (64)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res[stage == "64cell",]), aes(x = log2(SE), colour= as.factor(paste(condition, type)))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Condition & TISU status (WT vs 4Ei) (64)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# TISU vs background (64)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res[stage == "shield",]), aes(x = log2(SE), colour= as.factor(paste(condition, type)))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Condition & TISU status (WT vs 4Ei) (shield)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()

# Without pseudo counts
# TISU VS BACKGROUND
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res), aes(x = log2(SE), colour= as.factor(type))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("TISU vs background (stages merged)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# condition (WT vs 4Ei)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res), aes(x = log2(SE), colour= as.factor(condition))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Condition (WT vs 4Ei) (stages merged)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# Stage and condition
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res), aes(x = log2(SE), colour= as.factor(variable))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Stage and condiction")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# TISU VS BACKGROUND 64
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res[stage == "64cell",]), aes(x = log2(SE), colour= as.factor(type))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("TISU vs background (64)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# TISU VS BACKGROUND shield
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res[stage == "shield",]), aes(x = log2(SE), colour= as.factor(type))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("TISU vs background (shield)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()

coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(res), aes(x = log2(SE), colour= as.factor(variable))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf

# With pseudo
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(resp), aes(x = log2(SE), colour= as.factor(type))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("TISU vs background (stages merged)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# condition (WT vs 4Ei)
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(resp), aes(x = log2(SE), colour= as.factor(condition))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Condition (WT vs 4Ei) (stages merged)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# Stage and condition
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(resp), aes(x = log2(SE), colour= as.factor(variable))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("Stage and condiction")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# TISU VS BACKGROUND 64
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(resp[stage == "64cell",]), aes(x = log2(SE), colour= as.factor(type))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("TISU vs background (64)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()
# TISU VS BACKGROUND shield
coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(resp[stage == "shield",]), aes(x = log2(SE), colour= as.factor(type))) +
  stat_ecdf() +
  theme_classic(base_size = 13) +
  theme(legend.title = element_blank()) +
  ggtitle("TISU vs background (shield)")
coverage_plot_transcript_raw_100_ecdf
ggslackR()

coverage_plot_transcript_raw_100_ecdf <- ggplot(data=as.data.frame(resp), aes(x = log2(SE), colour= as.factor(variable))) +
  stat_ecdf()
coverage_plot_transcript_raw_100_ecdf

# With SE of WT removed if 0:

t.test(log2(res3[type == "WT",]$SE + 0.000001), log2(res3[type != "WT",]$SE + 0.000001), paired = TRUE, alternative = "two.sided")
wilcox.test(res3[type == "WT",]$SE, res3[type != "WT",]$SE, paired = TRUE)
wilcox.test(res2[type == "WT",]$SE, res2[type != "WT",]$SE, paired = TRUE)

