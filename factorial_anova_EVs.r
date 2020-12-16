
#install.packages("AICcmodavg")
library(AICcmodavg)
library(tidyverse)
library(ggpubr)
library(rstatix)

# older data loading, pre Controlled Differentiation
#all_size_std_a4to8_4to10_rois_v3

# newer data
clss <- c("factor", "factor", "factor", "factor", "factor", "numeric", "factor", "factor", "factor")
all_size_both_a4to8_4to10h_rois_v3 <- read.csv("freqs_size_all_v36_er3_both_a4-8_4-10h_rois_v3.0b.csv", colClasses = clss)
all_size_both_a4to8_4to10h_rois_v4 <- read.csv("freqs_all_size_v36_er3_both_a4-8_4-10h_rois_v4.0.csv", colClasses = clss)
all_size_cd_a4_4to8h_rois_v3 <- read.csv("freqs_all_size_v36_er3_cd_a4_4-8h_rois_v3.0.csv", colClasses = clss)
all_size_cd_a4_4to8h_rois_v4 <- read.csv("freqs_all_size_v36_er3_cd_a4_4-8h_rois_v4.0.csv", colClasses = clss)
ampulla_size_both_a4to8_4to10h_rois_v3 <- read.csv("freqs_ampulla_size_v36_er3_both_a4-8_4-10h_rois_v3.0.csv", colClasses = clss)
ampulla_size_both_a4to8_4to10h_rois_v4 <- read.csv("freqs_ampulla_size_v36_er3_both_a4-8_4-10h_rois_v4.0.csv", colClasses = clss)
ampulla_size_cd_a4to8_4to10h_rois_v3 <- read.csv("freqs_ampulla_size_v36_er3_cd_a4-8_4-10h_rois_v3.0.csv", colClasses = clss)
ampulla_size_cd_a4to8_4to10h_rois_v4 <- read.csv("freqs_ampulla_size_v36_er3_cd_a4-8_4-10h_rois_v4.0.csv", colClasses = clss)


# models - all_size_both_a4to8_4to10h_rois_v4
all_both_a4to8_4to10h_rois_v4.m1 <- lm(freq ~ section * stage * roi_type * time * type, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m2 <- lm(freq ~ section * stage * roi_type * time, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m3 <- lm(freq ~ section * stage * roi_type, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m4 <- lm(freq ~ section * stage * roi_type * type, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m5 <- lm(freq ~ section + stage + roi_type + time + type, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m6 <- lm(freq ~ section * stage * roi_type * apoptosis * differentiation * type, data = all_size_both_a4to8_4to10h_rois_v4)
all_both_a4to8_4to10h_rois_v4.m7 <- lm(freq ~ section + stage + roi_type + apoptosis + differentiation + type, data = all_size_both_a4to8_4to10h_rois_v4)
# blocking
all_both_a4to8_4to10h_rois_v4.m8 <- lm(freq ~ section * stage * roi_type * apoptosis * type + differentiation, data = all_size_both_a4to8_4to10h_rois_v4)

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
amp_s_both_a4to8_rois_v3.m1 <- lm(freq ~ stage * roi_type * time * type, data = ampulla_size_both_a4to8_4to10h_rois_v3)
amp_s_both_a4to8_rois_v3.m2 <- lm(freq ~ stage * roi_type * time, data = ampulla_size_both_a4to8_4to10h_rois_v3)
amp_s_both_a4to8_rois_v3.m3 <- lm(freq ~ stage * roi_type * type * type, data = ampulla_size_both_a4to8_4to10h_rois_v3)
amp_s_both_a4to8_rois_v3.m4 <- lm(freq ~ stage + roi_type + time + type, data = ampulla_size_both_a4to8_4to10h_rois_v3)
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

