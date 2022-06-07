
# Setup -------------------------------------------------------------------
setwd("~/Documents/University - Imperial/Introduction to Statistical Thinking and Data Analysis/Project 3")
library(survival); library(dplyr); library(survminer)

cancer_dat_initial <- read.csv('metabric-analytical.csv', header = T)
sum(cancer_dat_initial$age_diagnosis > 75)
cancer_dat <- cancer_dat_initial %>% filter(age_diagnosis <= 80)

head(cancer_dat)
summary(cancer_dat_initial)
cancer_dat$chemotherapy <- as.factor(cancer_dat$chemotherapy)
cancer_dat$hormone_therapy <- as.factor(cancer_dat$hormone_therapy)
cancer_dat$radio_therapy <- as.factor(cancer_dat$radio_therapy)
cancer_dat$breast_surgery <- as.factor(cancer_dat$breast_surgery)
cancer_dat$surv_status <- as.factor(cancer_dat$surv_status)
cancer_dat$cancer_type <- as.factor(cancer_dat$cancer_type)
cancer_dat$vital_status <- as.factor(cancer_dat$vital_status)
cancer_dat$cellularity <- as.factor(cancer_dat$cellularity)
# cancer_dat$cellularity <- ordered(cancer_dat$cellularity, levels = c("low", "moderate", "high"))
cancer_dat$grade <- as.factor(cancer_dat$grade)
cancer_dat$her2_status <- as.factor(cancer_dat$her2_status)
cancer_dat$pr_status <- as.factor(cancer_dat$pr_status)
cancer_dat$er_status <- as.factor(cancer_dat$er_status)

# EDA ---------------------------------------------------------------------

# How many NA values?

counter <- function(d){
  sum(is.na(d))
}
sapply(cancer_dat, counter)

# Tumor size =26
# Tumor grade = 88
# surv_statis = 0 
# Stage has 514 null values so we'll avoid analysing this
plot(cancer_dat$cellularity, cancer_dat$tumor_size)
table(cancer_dat$tumor_stage, cancer_dat$hormone_therapy)

table(cancer_dat$radio_therapy, cancer_dat$hormone_therapy)
table(cancer_dat$radio_therapy, cancer_dat$chemotherapy)
table( cancer_dat$hormone_therapy, cancer_dat$chemotherapy)
# I will look at size and grade, but one of these could be swapped for cellularity
summary(cancer_dat$vital_status)
summary(cancer_dat$age)
levels(cancer_dat$vital_status) <- c(1,0,0,NA)

summary(cancer_dat$surv_status)
summary(cancer_dat$surv_months)
summary(cancer_dat$cellularity)
summary(cancer_dat$grade)
summary(cancer_dat$tumor_size)
length(which(cancer_dat$tumor_size > 23))
hist(cancer_dat$tumor_size)
cancer_dat$size_bin <- ifelse(cancer_dat$tumor_size > 30, 1, 0)
# Making tumor size binary, bigger/ larger than 


levels(cancer_dat$surv_status) <- c(1, 0)
cancer_dat$surv_status <- as.numeric(cancer_dat$surv_status)
cancer_dat$vital_status <- as.numeric(cancer_dat$vital_status)

names(cancer_dat)
dim(cancer_dat)
cancer_dat <- cancer_dat[complete.cases(cancer_dat[ , c(2,4,8,10,14,15,16)]),]
dim(cancer_dat)

# vital_status
# Aim 1 -------------------------------------------------------------------
# ask gta
km_stage <- survfit(Surv(surv_months, vital_status) ~ 1, data = cancer_dat)
summary(km_stage)
km_stage
plot(km_stage, main = "Kaplan-Meier survivor function", xlab = "Months", ylab = "Proportion Alive")

summary(cancer_dat$surv_months)
# Maybe get the differential plot?
# Check significance with LR test
# Aim 2 -------------------------------------------------------------------

# Tumor size

size_lr <- coxph(Surv(surv_months, vital_status) ~ tumor_size, data = cancer_dat)
summary(size_lr)
ggcoxzph(cox.zph(size_lr))

# Cellularity

grade_c <- coxph(Surv(surv_months, vital_status) ~ cellularity, data = cancer_dat)
summary(grade_c)
ggcoxzph(cox.zph(grade_c))
# Neither are significant
# Not significant  


both <- coxph(Surv(surv_months, vital_status) ~ cellularity + tumor_size, data = cancer_dat)
summary(both)
ggcoxzph(cox.zph(both))

# Aim 3 -------------------------------------------------------------------

aim3 <- coxph(Surv(surv_months, vital_status) ~ chemotherapy + hormone_therapy + radio_therapy, data = cancer_dat)
summary(aim3)
cox.zph(aim3)
ggcoxzph(cox.zph(aim3))


km_rad <- survfit(Surv(surv_months, vital_status) ~ radio_therapy, data = cancer_dat)
plot(km_rad, col = 2:3, main = "Kaplan-Meier survivor function", xlab = "Months", ylab = "Proportion Alive")
legend("topright", c("No", "Yes"), col = 2:3, lty = 1)
summary(km_rad)

rad_split <- survSplit(Surv(surv_months,vital_status) ~ ., data = cancer_dat, cut = 120, episode = 'period')
rad_split$period <- factor(rad_split$period, 1:2, c('<120 months', '>120 months'))
split_fit <- coxph(Surv(surv_months, vital_status)  ~ chemotherapy + hormone_therapy + radio_therapy * strata(period), data = rad_split)

summary(split_fit)
cox.zph(split_fit)
ggcoxzph(cox.zph(split_fit))
# Lots of possibilities - could be confounded by a variable not included, it could have an interaction with other variables, try running kapalain myers of radio therapy



# Sensitivity analysis
sens <- coxph(Surv(surv_months, vital_status)  ~ chemotherapy + hormone_therapy + radio_therapy * strata(period) + tumor_size, data = rad_split)
summary(sens)
cox.zph(sens)

# Coefficients remain moderately constant, none of the significance really changes
# Interestingly grade becomes significant when before it was not 

a <- coxph(Surv(surv_months, vital_status)  ~ chemotherapy + hormone_therapy + radio_therapy * strata(period) + cellularity, data = rad_split)
summary(a)
cox.zph(a)

par(mfrow=c(2,2))
ggplot(cancer_dat, aes(x = radio_therapy, y = tumor_size)) + geom_boxplot() + coord_cartesian(ylim = c(0, 100))
ggplot(cancer_dat, aes(x = radio_therapy, y = tumor_size)) + geom_boxplot()+ coord_cartesian(ylim = c(0, 100))
ggplot(cancer_dat, aes(x = radio_therapy, y = tumor_size)) + geom_boxplot()+ coord_cartesian(ylim = c(0, 100))

plot(cancer_dat$radio_therapy, cancer_dat$tumor_size)
plot(cancer_dat$hormone_therapy, cancer_dat$tumor_size)
plot(cancer_dat$chemotherapy, cancer_dat$tumor_size)
cancer_dat$vital_status <- as.factor(cancer_dat$vital_status )
ggplot(cancer_dat, aes(x = hormone_therapy, fill = hormone_therapy)) + geom_bar()
sum(cancer_dat$hormone_therapy == 'no')

