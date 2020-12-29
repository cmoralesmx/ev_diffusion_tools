# This script uses Facotrial-ANOVA to analyse the data produced by the EV
# simulation model product of my PhD. Several experiments were completed, all
# are related so they all are analysed here

#install.packages("AICcmodavg")
library(AICcmodavg)
library(tidyverse)
library(ggpubr)
library(rstatix)

# The data is loaded from CSV, these are the classes to use for each column
cClasses <- c(
          "factor",
          "factor",
          "factor",
          "factor",
          "factor",
          "numeric",
          "factor",
          "factor",
          "factor")

# the CSV file names encode the following information
# freq
# oviductRegion, Possible values (Ampulla|Isthmus|All)
# size, The EV frequency was recorded per EV size
# oviductVersion, v36_er3: Version 36 extended range 3
# mode of cell type determination (cd - controlled|blank for random|both modes)
# a4, a8, a4-8, Expected lifespan in hours (4,8,4-8: means the data contains both values)
# 4-10h, simulated time in hours, 4-10h contains data at 4,8,10h
# ROIs version, version of the Regions of interest used (v3.0, v4.0)

# This csv1 is located in ~/data/ev_iter/analysis/ and was produced by
# concatenating the content of the files in
# ~/data/ev_iter/analysis/freq_size_all_both_a4-8_4-10h_rois_v3
csv1 <- "freqs_size_all_v36_er3_both_a4-8_4-10h_rois_v3.0b.csv"

# Evidence of how the following files were produced is available in
# ~/data/ev_iter/ev_repeats.txt, the comment per line shows the line number
csv2 <- "freqs_all_size_v36_er3_both_a4-8_4-10h_rois_v4.0.csv"  # ln 290
csv3 <- "freqs_all_size_v36_er3_cd_a4_4-8h_rois_v3.0.csv" # ln 234
csv4 <- "freqs_all_size_v36_er3_cd_a4_4-8h_rois_v4.0.csv" # ln 277
csv5 <- "freqs_ampulla_size_v36_er3_both_a4-8_4-10h_rois_v3.0.csv" # ln 248
csv6 <- "freqs_ampulla_size_v36_er3_both_a4-8_4-10h_rois_v4.0.csv" # ln 290
csv7 <- "freqs_ampulla_size_v36_er3_cd_a4-8_4-10h_rois_v3.0.csv" # ln 242
csv8 <- "freqs_ampulla_size_v36_er3_cd_a4-8_4-10h_rois_v4.0.csv" # ln 285

csv1.df <- read.csv(csv1, colClasses = cClasses)
csv2.df <- read.csv(csv2, colClasses = cClasses)
csv3.df <- read.csv(csv3, colClasses = cClasses)
csv4.df <- read.csv(csv4, colClasses = cClasses)
csv5.df <- read.csv(csv5, colClasses = cClasses)
csv6.df <- read.csv(csv6, colClasses = cClasses)
csv7.df <- read.csv(csv7, colClasses = cClasses)
csv8.df <- read.csv(csv8, colClasses = cClasses)

# formulae for analyzing the data in csv1 & csv2
f1 <- freq ~ section * stage * roi_type * time * type
f2 <- freq ~ section * stage * roi_type * time
f3 <- freq ~ section * stage * roi_type
f4 <- freq ~ section * stage * roi_type * type
f5 <- freq ~ section + stage + roi_type + time + type
f6 <- freq ~ section * stage * roi_type * apoptosis * differentiation * type
f7 <- freq ~ section + stage + roi_type + apoptosis + differentiation + type
f8 <- freq ~ section * stage * roi_type * apoptosis * type + differentiation
# formulae for csv1 only
f9 <- freq ~ section * stage * roi_type * time * apoptosis * differentiation * type
f10 <- freq ~ section * stage * roi_type + time * apoptosis * differentiation * type
f11 <- freq ~ section * stage * roi_type * apoptosis * differentiation * type + time
f12 <- freq ~ section * stage * roi_type * differentiation * type
f13 <- freq ~ section * stage + roi_type * differentiation * apoptosis * type
f14 <- freq ~ section + stage + roi_type + time + apoptosis + differentiation + type
f15 <- freq ~ section * stage * roi_type * time * apoptosis * type + differentiation
# csv4
f16 <- freq ~ section * stage * roi_type * time * apoptosis *  type
# csv5
f17 <- freq ~ stage * roi_type * time * type
f18 <- freq ~ stage * roi_type * time
f19 <- freq ~ stage * roi_type * type * time
f20 <- freq ~ stage + roi_type + time + type
# csv6
f21 <- freq ~ stage * roi_type * time * apoptosis * type
# csv7
f22 <- freq ~ stage * roi_type * apoptosis * type
f23 <- freq ~ stage * roi_type * apoptosis
f24 <- freq ~ stage * roi_type * time * apoptosis

# linear models for csv1.df
csv1.m1  <- lm(f1, data = csv1.df)
csv1.m2  <- lm(f2, data = csv1.df)
csv1.m3  <- lm(f3, data = csv1.df)
csv1.m4  <- lm(f4, data = csv1.df)
csv1.m5  <- lm(f5, data = csv1.df)
csv1.m6  <- lm(f6, data = csv1.df) # best model?
csv1.m7  <- lm(f7, data = csv1.df)
csv1.m8  <- lm(f8, data = csv1.df)
# alternate models
csv1.m9  <- lm(f9, data = csv1.df)
csv1.m10 <- lm(f10, data = csv1.df)
csv1.m11 <- lm(f11, data = csv1.df)
csv1.m12 <- lm(f12, data = csv1.df)
csv1.m13 <- lm(f13, data = csv1.df)
csv1.m14 <- lm(f14, data = csv1.df)
# blocking
csv1.m15 <- lm(f15, data = csv1.df)

csv1.models <- list(csv1.m1, csv1.m2, csv1.m3, csv1.m4, csv1.m5, csv1.m6,
                    csv1.m7, csv1.m8, csv1.m9, csv1.m10, csv1.m11, csv1.m12,
                    csv1.m13, csv1.m14, csv1.m15)
csv1.models.names <- c("csv1.m1", "csv1.m2", "csv1.m3", "csv1.m4", "csv1.m5",
                       "csv1.m6", "csv1.m7", "csv1.m8", "csv1.m9", "csv1.m10",
                       "csv1.m11", "csv1.m12", "csv1.m13", "csv1.m14",
                       "csv1.m15")
aictab(cand.set = csv1.models, modnames = csv1.models.names)

# linear models for csv2
csv2.m1 <- lm(f1, data = csv2.df)
csv2.m2 <- lm(f2, data = csv2.df)
csv2.m3 <- lm(f3, data = csv2.df)
csv2.m4 <- lm(f4, data = csv2.df)
csv2.m5 <- lm(f5, data = csv2.df)
csv2.m6 <- lm(f6, data = csv2.df)
csv2.m7 <- lm(f7, data = csv2.df)
csv2.m8 <- lm(f8, data = csv2.df) # blocking

csv2.models <- list(csv2.m1, csv2.m2, csv2.m3, csv2.m4, csv2.m5, csv2.m6,
                    csv2.m7, csv2.m8)
csv2.models.names <- c("csv2.m1", "csv2.m2", "csv2.m3", "csv2.m4",
                       "csv2.m5", "csv2.m6", "csv2.m7", "csv2.m8")
aictab(cand.set = csv2.models, modnames = csv2.models.names)

# models - csv3.df
csv3.m1 <- lm(f1, data = csv3.df)
csv3.m2 <- lm(f2, data = csv3.df)
csv3.m3 <- lm(f3, data = csv3.df)
csv3.m4 <- lm(f4, data = csv3.df)
csv3.m5 <- lm(f5, data = csv3.df)

csv3.models <- list(csv3.m1, csv3.m2, csv3.m3, csv3.m4, csv3.m5)
csv3.models.names <- c("csv3.m1", "csv3.m2", "csv3.m3", "csv3.m4",
                       "csv3.m5")
aictab(cand.set = csv3.models, modnames = csv3.models.names)

# models - csv4.df
csv4.m1 <- lm(f1, data = csv4.df)
csv4.m2 <- lm(f2, data = csv4.df)
csv4.m3 <- lm(f3, data = csv4.df)
csv4.m4 <- lm(f4, data = csv4.df)
csv4.m5 <- lm(f5, data = csv4.df)
csv4.m6 <- lm(f16, data = csv4.df)
csv4.models <- list(csv4.m1, csv4.m2, csv4.m3, csv4.m4, csv4.m5, csv4.m6)
csv4.models.names <- c("csv4.m1", "csv4.m2", "csv4.m3", "csv4.m4", "csv4.m5",
                       "csv4.m6")
aictab(cand.set = csv4.models, modnames = csv4.models.names)

# models for csv5.df
csv5.m1 <- lm(f17, data = csv5.df)
csv5.m2 <- lm(f18, data = csv5.df)
csv5.m3 <- lm(f19, data = csv5.df)
csv5.m4 <- lm(f20, data = csv5.df)

csv5.models <- list(csv5.m1, csv5.m2, csv5.m3, csv5.m4)
csv5.models.names <- c("csv5.m1", "csv5.m2", "csv5.m3", "csv5.m4")
aictab(cand.set = csv.models, modnames = csv5.models.names)

# models for csv6.df
csv6.m1 <- lm(f17, data = csv6.df)
csv6.m2 <- lm(f18, data = csv6.df)
csv6.m3 <- lm(f19, data = csv6.df)
csv6.m4 <- lm(f20, data = csv6.df)
csv6.m5 <- lm(f21, data = csv6.df)
csv6.models <- list(csv6.m1, csv6.m2, csv6.m3, csv6.m4, csv6.m5)
csv6.models.names <- c("csv6.m1", "csv6.m2", "csv6.m3", "csv6.m4", "csv6.m5")
aictab(cand.set = csv6.models, modnames = csv6.models.names)

# models for csv7.df
csv7.m1 <- lm(f17, data = csv7.df)
csv7.m2 <- lm(f18, data = csv7.df)
csv7.m3 <- lm(f19, data = csv7.df)
csv7.m4 <- lm(f20, data = csv7.df)
csv7.m5 <- lm(f21, data = csv7.df)
#
csv7.m6 <- lm(f22, data = csv7.df)
csv7.m7 <- lm(f23, data = csv7.df)
csv7.m8 <- lm(f24, data = csv7.df)
csv7.models <- list(csv7.m1, csv7.m2, csv7.m3, csv7.m4, csv7.m5, csv7.m6,
                    csv7.m7, csv7.m8)
csv7.models.names <- c("csv7.m1", "csv7.m2", "csv7.m3", "csv7.m4", "csv7.m5",
                       "csv7.m6", "csv7.m7", "csv7.m8")
aictab(cand.set = csv7.models, modnames = csv7.models.names)

# models for csv8.df
csv8.m1 <- lm(f17, data = csv8.df)
csv8.m2 <- lm(f18, data = csv8.df)
csv8.m3 <- lm(f19, data = csv8.df)
csv8.m4 <- lm(f20, data = csv8.df)
csv8.m5 <- lm(f21, data = csv8.df)
csv8.models <- list(csv8.m1, csv8.m2, csv8.m3, csv8.m4, csv8.m5)
csv8.models.names <- c("csv8.m1", "csv8.m2", "csv8.m3", "csv8.m4", "csv8.m5")
aictab(cand.set = csv8.models, modnames = csv8.models.names)
