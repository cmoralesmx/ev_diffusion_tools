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


# formulae for analyzing the data in csv2 (freqs_data_csv2)
csv2.f1 <- freq ~ section * stage * roi_type * time * type
csv2.f2 <- freq ~ section * stage * roi_type * time
csv2.f3 <- freq ~ section * stage * roi_type
csv2.f4 <- freq ~ section * stage * roi_type * type
csv2.f5 <- freq ~ section + stage + roi_type + time + type
csv2.f6 <- freq ~ section * stage * roi_type * apoptosis * differentiation * type
csv2.f7 <- freq ~ section + stage + roi_type + apoptosis + differentiation + type
csv2.f8 <- freq ~ section * stage * roi_type * apoptosis * type + differentiation
# linear models for csv2
csv2.m1 <- lm(csv2.f1, data = csv2.df)
csv2.m2 <- lm(csv2.f2, data = csv2.df)
csv2.m3 <- lm(csv2.f3, data = csv2.df)
csv2.m4 <- lm(csv2.f4, data = csv2.df)
csv2.m5 <- lm(csv2.f5, data = csv2.df)
csv2.m6 <- lm(csv2.f6, data = csv2.df)
csv2.m7 <- lm(csv2.f7, data = csv2.df)
# blocking
csv2.m8 <- lm(csv2.f8, data = csv2.df)

csv2.models <- list(csv2.m1, csv2.m2, csv2.m3, csv2.m4,
                    csv2.m5, csv2.m6,csv2.m7, csv2.m8)
csv2.models.names <- c("csv2.m1", "csv2.m2", "csv2.m3", "csv2.m4",
                       "csv2.m5", "csv2.m6", "csv2.m7", "csv2.m8")
aictab(cand.set = csv2.models, modnames = csv2.models.names)

# models - csv1.df
csv1.m1  <- lm(csv2.f1, data = csv1.df)
csv1.m2  <- lm(csv2.f2, data = csv1.df)
csv1.m3  <- lm(csv2.f3, data = csv1.df)
csv1.m4  <- lm(csv2.f4, data = csv1.df)
csv1.m5  <- lm(csv2.f5, data = csv1.df)
csv1.m6  <- lm(csv2.f6, data = csv1.df) # best model? same formula as csv2.f6

# alternate models
csv1.m7  <- lm(freq ~ section * stage * roi_type * time * apoptosis * differentiation * type, data = csv1.df)
csv1.m8  <- lm(freq ~ section * stage * roi_type + time * apoptosis * differentiation * type, data = csv1.df)
csv1.m9  <- lm(freq ~ section * stage * roi_type * apoptosis * differentiation * type + time, data = csv1.df)
csv1.m10 <- lm(freq ~ section * stage * roi_type * differentiation * type, data = csv1.df)
csv1.m11 <- lm(freq ~ section * stage + roi_type * differentiation * apoptosis * type, data = csv1.df)
csv1.m12  <- lm(freq ~ section + stage + roi_type + time + apoptosis + differentiation + type, data = csv1.df)
# blocking

csv1.m13  <- lm(csv2.f8, data = csv1.df) # same formula as csv2.f8
csv1.m14  <- lm(freq ~ section * stage * roi_type * time * apoptosis * type + differentiation, data = csv1.df)

models_all_both_v3 <- list(csv1.m1, csv1.m2,
                           csv1.m3, csv1.m4,
                           csv1.m5, csv1.m6,
                           csv1.m7, csv1.m8, csv1.m9,
                           csv1.m10, csv1.m11, csv1.m12,
                           csv1.m13, csv1.m14)
models_all_both_v3.names <- c("csv1.m1", "csv1.m2", "csv1.m3", "csv1.m4", "csv1.m5", "csv1.m6", "csv1.m7", "csv1.m8", "csv1.m9", "csv1.m10", "csv1.m11", "csv1.m12", "csv1.m13", "csv1.m14")
aictab(cand.set = models_all_both_v3, modnames = models_all_both_v3.names)

# models - csv3.df
csv3.m1 <- lm(csv2.f1, data = csv3.df)
csv3.m2 <- lm(csv2.f2, data = csv3.df)
csv3.m3 <- lm(csv2.f3, data = csv3.df)
csv3.m4 <- lm(csv2.f4, data = csv3.df)
csv3.m5 <- lm(csv2.f5, data = csv3.df)

models_all_cd_v3 <- list(csv3.m1, csv3.m2,
                         csv3.m3, csv3.m4,
                         csv3.m5)
models_all_cd_v3.names <- c("csv3.m1", "csv3.m2", "csv3.m3", "csv3.m4", "csv3.m5")
aictab(cand.set = models_all_cd_v3, modnames = models_all_cd_v3.names)

# models - csv4.df
csv4.m1 <- lm(csv2.f1, data = csv4.df)
csv4.m2 <- lm(csv2.f2, data = csv4.df)
csv4.m3 <- lm(csv2.f3, data = csv4.df)
csv4.m4 <- lm(csv2.f4, data = csv4.df)
csv4.m5 <- lm(csv2.f5, data = csv4.df)
csv4.m1 <- lm(freq ~ section * stage * roi_type * time * apoptosis *  type, data = csv4.df)
models_all_cd_v4 <- list(csv4.m1, csv4.m2,
                         csv4.m3, csv4.m4,
                         csv4.m5)
models_all_cd_v4.names <- c("csv4.m1", "csv4.m2", "csv4.m3", "csv4.m4", "csv4.m5")
aictab(cand.set = models_all_cd_v4, modnames = models_all_cd_v4.names)

# models for csv5.df
m1 <- freq ~ stage * roi_type * time * type
m2 <- freq ~ stage * roi_type * time
m3 <- freq ~ stage * roi_type * type * type
m4 <- freq ~ stage + roi_type + time + type
csv5.m1 <- lm(m1, data = csv5.df)
csv5.m2 <- lm(m2, data = csv5.df)
csv5.m3 <- lm(m3, data = csv5.df)
csv5.m4 <- lm(m4, data = csv5.df)

models_amp_s_both_v3 <- list(csv5.m1, csv5.m2,
                             csv5.m3, csv5.m4)
models_amp_s_both_v3.names <- c("csv5.m1", "csv5.m2", "csv5.m3", "csv5.m4")
aictab(cand.set = models_amp_s_both_v3, modnames = models_amp_s_both_v3.names)

# models for csv6.df
csv6.m1 <- lm(freq ~ stage * roi_type * time * type, data = csv6.df)
csv6.m2 <- lm(freq ~ stage * roi_type * time, data = csv6.df)
csv6.m3 <- lm(freq ~ stage * roi_type * type * type, data = csv6.df)
csv6.m4 <- lm(freq ~ stage + roi_type + time + type, data = csv6.df)
csv6.m5 <- lm(freq ~ stage * roi_type * time * apoptosis * type, data = csv6.df)
models_amp_s_both_v4 <- list(csv6.m1, csv6.m2,
                             csv6.m3, csv6.m4, csv6.m5)
models_amp_s_both_v4.names <- c("csv6.m1", "csv6.m2", "csv6.m3", "csv6.m4",
                                "csv6.m5")
aictab(cand.set = models_amp_s_both_v4, modnames = models_amp_s_both_v4.names)

# models for csv7.df
csv7.m1 <- lm(freq ~ stage * roi_type * time * type, data = csv7.df)
csv7.m2 <- lm(freq ~ stage * roi_type * time, data = csv7.df)
csv7.m3 <- lm(freq ~ stage * roi_type * type * type, data = csv7.df)
csv7.m4 <- lm(freq ~ stage + roi_type + time + type, data = csv7.df)
csv7.m5 <- lm(freq ~ stage * roi_type * time * apoptosis * type, data = csv7.df)
csv7.m6 <- lm(freq ~ stage * roi_type * apoptosis * type, data = csv7.df)
csv7.m7 <- lm(freq ~ stage * roi_type * apoptosis, data = csv7.df)
csv7.m8 <- lm(freq ~ stage * roi_type * time * apoptosis, data = csv7.df)
models_amp_s_cd_v3 <- list(csv7.m1, csv7.m2,
                           csv7.m3, csv7.m4,
                           csv7.m5, csv7.m6,
                           csv7.m7, csv7.m8)
models_amp_s_cd_v3.names <- c("csv7.m1", "csv7.m2", "csv7.m3", "csv7.m4",
                              "csv7.m5", "csv7.m6", "csv7.m7", "csv7.m8")
aictab(cand.set = models_amp_s_cd_v3, modnames = models_amp_s_cd_v3.names)

# models for csv8.df
csv8.m1 <- lm(freq ~ stage * roi_type * time * type, data = csv8.df)
csv8.m2 <- lm(freq ~ stage * roi_type * time, data = csv8.df)
csv8.m3 <- lm(freq ~ stage * roi_type * type * type, data = csv8.df)
csv8.m4 <- lm(freq ~ stage + roi_type + time + type, data = csv8.df)
csv8.m5 <- lm(freq ~ stage * roi_type * time * apoptosis * type, data = csv8.df)
models_amp_s_cd_v4 <- list(csv8.m1, csv8.m2,
                             csv8.m3, csv8.m4, csv8.m5)
models_amp_s_cd_v4.names <- c("csv8.m1", "csv8.m2", "csv8.m3", "csv8.m4",
                              "csv8.m5")
aictab(cand.set = models_amp_s_cd_v4, modnames = models_amp_s_cd_v4.names)
