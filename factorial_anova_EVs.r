
#install.packages("AICcmodavg")
library(AICcmodavg)
library(tidyverse)
library(ggpubr)
library(rstatix)

# older data loading, pre Controlled Differentiation
#all_size_std_a4to8_4to10_rois_v3

# newer data
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
all_size_both_a4to8_4to10h_rois_v3     <- read.csv(csv1, colClasses = cClasses)
all_size_both_a4to8_4to10h_rois_v4     <- read.csv(csv2, colClasses = cClasses)
all_size_cd_a4_4to8h_rois_v3           <- read.csv(csv3, colClasses = cClasses)
all_size_cd_a4_4to8h_rois_v4           <- read.csv(csv4, colClasses = cClasses)
ampulla_size_both_a4to8_4to10h_rois_v3 <- read.csv(csv5, colClasses = cClasses)
ampulla_size_both_a4to8_4to10h_rois_v4 <- read.csv(csv6, colClasses = cClasses)
ampulla_size_cd_a4to8_4to10h_rois_v3   <- read.csv(csv7, colClasses = cClasses)
ampulla_size_cd_a4to8_4to10h_rois_v4   <- read.csv(csv8, colClasses = cClasses)


# models for analyzing the data in freqs_all_size_both_a4to8_4to10h_rois_v4
m1_formula <- freq ~ section * stage * roi_type * time * type
m2_formula <- freq ~ section * stage * roi_type * time
m3_formula <- freq ~ section * stage * roi_type
m4_formula <- freq ~ section * stage * roi_type * type
m5_formula <- freq ~ section + stage + roi_type + time + type
m6_formula <- freq ~ section * stage * roi_type * apoptosis * differentiation * type
m7_formula <- freq ~ section + stage + roi_type + apoptosis + differentiation + type
m8_formula <- freq ~ section * stage * roi_type * apoptosis * type + differentiation

all_both_a4to8_4to10h_rois_v4.m1 <- lm(m1_formula, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m2 <- lm(m2_formula, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m3 <- lm(m3_formula, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m4 <- lm(m4_formula, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m5 <- lm(m5_formula, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m6 <- lm(m6_formula, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m7 <- lm(m7_formula, data = all_size_both_a4to8_4to10h_rois_v4)
# blocking
all_both_a4to8_4to10h_rois_v4.m8 <- lm(m8_formula, data = all_size_both_a4to8_4to10h_rois_v4)

models_all_both_v4 <- list(all_both_a4to8_4to10h_rois_v4.m1, all_both_a4to8_4to10h_rois_v4.m2,
                           all_both_a4to8_4to10h_rois_v4.m3, all_both_a4to8_4to10h_rois_v4.m4,
                           all_both_a4to8_4to10h_rois_v4.m5, all_both_a4to8_4to10h_rois_v4.m6,
                           all_both_a4to8_4to10h_rois_v4.m7, all_both_a4to8_4to10h_rois_v4.m8)
models_all_both_v4.names <- c("a_v4_m1", "a_v4_m2", "a_v4_m3", "a_v4_m4", "a_v4_m5", "a_v4_m6", "a_v4_m7", "a_v4_m8")
aictab(cand.set = models_all_both_v4, modnames = models_all_both_v4.names)

# models - all_size_both_a4to8_4to10h_rois_v3
all_both_a4to8_4to10h_rois_v3.m1  <- lm(freq ~ section * stage * roi_type * time * type, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m2  <- lm(freq ~ section * stage * roi_type * time, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m3  <- lm(freq ~ section * stage * roi_type, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m4  <- lm(freq ~ section * stage * roi_type * type, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m5  <- lm(freq ~ section + stage + roi_type + time + type, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m6  <- lm(freq ~ section * stage * roi_type * time * apoptosis * differentiation * type, data = all_size_both_a4to8_4to10h_rois_v3)
# best model?
all_both_a4to8_4to10h_rois_v3.m7  <- lm(freq ~ section * stage * roi_type * apoptosis * differentiation * type, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m8  <- lm(freq ~ section * stage * roi_type + time * apoptosis * differentiation * type, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m9  <- lm(freq ~ section * stage * roi_type * apoptosis * differentiation * type + time, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m10 <- lm(freq ~ section * stage * roi_type * differentiation * type, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m11 <- lm(freq ~ section * stage + roi_type * differentiation * apoptosis * type, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m12  <- lm(freq ~ section + stage + roi_type + time + apoptosis + differentiation + type, data = all_size_both_a4to8_4to10h_rois_v3)
# blocking
all_both_a4to8_4to10h_rois_v3.m13  <- lm(freq ~ section * stage * roi_type * apoptosis * type + differentiation, data = all_size_both_a4to8_4to10h_rois_v3)
all_both_a4to8_4to10h_rois_v3.m14  <- lm(freq ~ section * stage * roi_type * time * apoptosis * type + differentiation, data = all_size_both_a4to8_4to10h_rois_v3)

models_all_both_v3 <- list(all_both_a4to8_4to10h_rois_v3.m1, all_both_a4to8_4to10h_rois_v3.m2,
                           all_both_a4to8_4to10h_rois_v3.m3, all_both_a4to8_4to10h_rois_v3.m4,
                           all_both_a4to8_4to10h_rois_v3.m5, all_both_a4to8_4to10h_rois_v3.m6,
                           all_both_a4to8_4to10h_rois_v3.m7, all_both_a4to8_4to10h_rois_v3.m8, all_both_a4to8_4to10h_rois_v3.m9,
                           all_both_a4to8_4to10h_rois_v3.m10, all_both_a4to8_4to10h_rois_v3.m11, all_both_a4to8_4to10h_rois_v3.m12,
                           all_both_a4to8_4to10h_rois_v3.m13, all_both_a4to8_4to10h_rois_v3.m14)
models_all_both_v3.names <- c("a_v3_m1", "a_v3_m2", "a_v3_m3", "a_v3_m4", "a_v3_m5", "a_v3_m6", "a_v3_m7", "a_v3_m8", "a_v3_m9", "a_v3_m10", "a_v3_m11", "a_v3_m12", "a_v3_m13", "a_v3_m14")
aictab(cand.set = models_all_both_v3, modnames = models_all_both_v3.names)

# models - all_size_cd_a4_4to8h_rois_v3
all_cd_a4_4to8h_rois_v3.m1 <- lm(freq ~ section * stage * roi_type * time * type, data = all_size_cd_a4_4to8h_rois_v3)
all_cd_a4_4to8h_rois_v3.m2 <- lm(freq ~ section * stage * roi_type * time, data = all_size_cd_a4_4to8h_rois_v3)
all_cd_a4_4to8h_rois_v3.m3 <- lm(freq ~ section * stage * roi_type, data = all_size_cd_a4_4to8h_rois_v3)
all_cd_a4_4to8h_rois_v3.m4 <- lm(freq ~ section * stage * roi_type * type, data = all_size_cd_a4_4to8h_rois_v3)
all_cd_a4_4to8h_rois_v3.m5 <- lm(freq ~ section + stage + roi_type + time + type, data = all_size_cd_a4_4to8h_rois_v3)

models_all_cd_v3 <- list(all_cd_a4_4to8h_rois_v3.m1, all_cd_a4_4to8h_rois_v3.m2,
                         all_cd_a4_4to8h_rois_v3.m3, all_cd_a4_4to8h_rois_v3.m4,
                         all_cd_a4_4to8h_rois_v3.m5)
models_all_cd_v3.names <- c("m1", "m2", "m3", "m4", "m5")
aictab(cand.set = models_all_cd_v3, modnames = models_all_cd_v3.names)

# models - all_size_cd_a4_4to8h_rois_v4
all_cd_a4_4to8h_rois_v4.m1 <- lm(freq ~ section * stage * roi_type * time * type, data = all_size_cd_a4_4to8h_rois_v4)
all_cd_a4_4to8h_rois_v4.m2 <- lm(freq ~ section * stage * roi_type * time, data = all_size_cd_a4_4to8h_rois_v4)
all_cd_a4_4to8h_rois_v4.m3 <- lm(freq ~ section * stage * roi_type, data = all_size_cd_a4_4to8h_rois_v4)
all_cd_a4_4to8h_rois_v4.m4 <- lm(freq ~ section * stage * roi_type * type, data = all_size_cd_a4_4to8h_rois_v4)
all_cd_a4_4to8h_rois_v4.m5 <- lm(freq ~ section + stage + roi_type + time + type, data = all_size_cd_a4_4to8h_rois_v4)
all_cd_a4_4to8h_rois_v4.m1 <- lm(freq ~ section * stage * roi_type * time * apoptosis *  type, data = all_size_cd_a4_4to8h_rois_v4)
models_all_cd_v4 <- list(all_cd_a4_4to8h_rois_v4.m1, all_cd_a4_4to8h_rois_v4.m2,
                         all_cd_a4_4to8h_rois_v4.m3, all_cd_a4_4to8h_rois_v4.m4,
                         all_cd_a4_4to8h_rois_v4.m5)
models_all_cd_v4.names <- c("m1", "m2", "m3", "m4", "m5")
aictab(cand.set = models_all_cd_v4, modnames = models_all_cd_v4.names)

# models for ampulla_size_both_a4to8_4to10h_rois_v3
m1 <- freq ~ stage * roi_type * time * type
m2 <- freq ~ stage * roi_type * time
m3 <- freq ~ stage * roi_type * type * type
m4 <- freq ~ stage + roi_type + time + type
amp_s_both_a4to8_rois_v3.m1 <- lm(m1, data = ampulla_size_both_a4to8_4to10h_rois_v3)
amp_s_both_a4to8_rois_v3.m2 <- lm(m2, data = ampulla_size_both_a4to8_4to10h_rois_v3)
amp_s_both_a4to8_rois_v3.m3 <- lm(m3, data = ampulla_size_both_a4to8_4to10h_rois_v3)
amp_s_both_a4to8_rois_v3.m4 <- lm(m4, data = ampulla_size_both_a4to8_4to10h_rois_v3)

models_amp_s_both_v3 <- list(amp_s_both_a4to8_rois_v3.m1, amp_s_both_a4to8_rois_v3.m2,
                             amp_s_both_a4to8_rois_v3.m3, amp_s_both_a4to8_rois_v3.m4)
models_amp_s_both_v3.names <- c("m1", "m2", "m3", "m4")
aictab(cand.set = models_amp_s_both_v3, modnames = models_amp_s_both_v3.names)

# models for ampulla_size_both_a4to8_4to10h_rois_v4
amp_s_both_a4to8_rois_v4.m1 <- lm(freq ~ stage * roi_type * time * type, data = ampulla_size_both_a4to8_4to10h_rois_v4)
amp_s_both_a4to8_rois_v4.m2 <- lm(freq ~ stage * roi_type * time, data = ampulla_size_both_a4to8_4to10h_rois_v4)
amp_s_both_a4to8_rois_v4.m3 <- lm(freq ~ stage * roi_type * type * type, data = ampulla_size_both_a4to8_4to10h_rois_v4)
amp_s_both_a4to8_rois_v4.m4 <- lm(freq ~ stage + roi_type + time + type, data = ampulla_size_both_a4to8_4to10h_rois_v4)
amp_s_both_a4to8_rois_v4.m5 <- lm(freq ~ stage * roi_type * time * apoptosis * type, data = ampulla_size_both_a4to8_4to10h_rois_v4)
models_amp_s_both_v4 <- list(amp_s_both_a4to8_rois_v4.m1, amp_s_both_a4to8_rois_v4.m2,
                             amp_s_both_a4to8_rois_v4.m3, amp_s_both_a4to8_rois_v4.m4, amp_s_both_a4to8_rois_v4.m5)
models_amp_s_both_v4.names <- c("m1", "m2", "m3", "m4", "m5")
aictab(cand.set = models_amp_s_both_v4, modnames = models_amp_s_both_v4.names)

# models for ampulla_size_cd_a4to8_4to10h_rois_v3
amp_s_cd_a4to8_rois_v3.m1 <- lm(freq ~ stage * roi_type * time * type, data = ampulla_size_cd_a4to8_4to10h_rois_v3)
amp_s_cd_a4to8_rois_v3.m2 <- lm(freq ~ stage * roi_type * time, data = ampulla_size_cd_a4to8_4to10h_rois_v3)
amp_s_cd_a4to8_rois_v3.m3 <- lm(freq ~ stage * roi_type * type * type, data = ampulla_size_cd_a4to8_4to10h_rois_v3)
amp_s_cd_a4to8_rois_v3.m4 <- lm(freq ~ stage + roi_type + time + type, data = ampulla_size_cd_a4to8_4to10h_rois_v3)
amp_s_cd_a4to8_rois_v3.m5 <- lm(freq ~ stage * roi_type * time * apoptosis * type, data = ampulla_size_cd_a4to8_4to10h_rois_v3)
amp_s_cd_a4to8_rois_v3.m6 <- lm(freq ~ stage * roi_type * apoptosis * type, data = ampulla_size_cd_a4to8_4to10h_rois_v3)
amp_s_cd_a4to8_rois_v3.m7 <- lm(freq ~ stage * roi_type * apoptosis, data = ampulla_size_cd_a4to8_4to10h_rois_v3)
amp_s_cd_a4to8_rois_v3.m8 <- lm(freq ~ stage * roi_type * time * apoptosis, data = ampulla_size_cd_a4to8_4to10h_rois_v3)
models_amp_s_cd_v3 <- list(amp_s_cd_a4to8_rois_v3.m1, amp_s_cd_a4to8_rois_v3.m2,
                           amp_s_cd_a4to8_rois_v3.m3, amp_s_cd_a4to8_rois_v3.m4,
                           amp_s_cd_a4to8_rois_v3.m5, amp_s_cd_a4to8_rois_v3.m6,
                           amp_s_cd_a4to8_rois_v3.m7, amp_s_cd_a4to8_rois_v3.m8)
models_amp_s_cd_v3.names <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8")
aictab(cand.set = models_amp_s_cd_v3, modnames = models_amp_s_cd_v3.names)

# models for ampulla_size_cd_a4to8_4to10h_rois_v4
amp_s_cd_a4to8_rois_v4.m1 <- lm(freq ~ stage * roi_type * time * type, data = ampulla_size_cd_a4to8_4to10h_rois_v4)
amp_s_cd_a4to8_rois_v4.m2 <- lm(freq ~ stage * roi_type * time, data = ampulla_size_cd_a4to8_4to10h_rois_v4)
amp_s_cd_a4to8_rois_v4.m3 <- lm(freq ~ stage * roi_type * type * type, data = ampulla_size_cd_a4to8_4to10h_rois_v4)
amp_s_cd_a4to8_rois_v4.m4 <- lm(freq ~ stage + roi_type + time + type, data = ampulla_size_cd_a4to8_4to10h_rois_v4)
amp_s_cd_a4to8_rois_v4.m5 <- lm(freq ~ stage * roi_type * time * apoptosis * type, data = ampulla_size_cd_a4to8_4to10h_rois_v4)
models_amp_s_cd_v4 <- list(amp_s_cd_a4to8_rois_v4.m1, amp_s_cd_a4to8_rois_v4.m2,
                             amp_s_cd_a4to8_rois_v4.m3, amp_s_cd_a4to8_rois_v4.m4, amp_s_cd_a4to8_rois_v4.m5)
models_amp_s_cd_v4.names <- c("m1", "m2", "m3", "m4", "m5")
aictab(cand.set = models_amp_s_cd_v4, modnames = models_amp_s_cd_v4.names)
