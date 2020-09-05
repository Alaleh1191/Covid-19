# Set the Directory and read in the data
setwd("~/Desktop/Thesis/Method2/")
data = read.csv("covid_D.csv")
library(stringr) # used for data processing
library(boot) # used for Tables
library(table1) # used for Tables
library(ggplot2) # used for plotting
library(ggpubr) # used for plotting
library(foreach) # parallelizing
library(data.table) # parallelizing
library(doParallel) # parallelizing
library(splines) # for fitting a cubic spline
library(ranger) # used for random forest classification
library(PRROC) # ROC and Precision Recall 
library(pROC) # ROC and Precision Recall 
library(caret) # Used for model fine tuning

# Data Preprocessing ------------------------------------------------------
### Death marker: corresponding to bad outcome, grouping hospice care with expired patients
data$death = data$DECEASED_INDICATOR
data$death[data$DISCHARGE_LOCATON == 'EXPIRED' | data$DISCHARGE_LOCATON == 'HOSPICE CARE' | 
             data$DISCHARGE_LOCATON == 'HOSPICE DISCHARGE HOME' | 
             data$DISCHARGE_LOCATON == 'HOSPICE DISCHARGE TO INSTITUTION'] = 1

# used for converting string array entries into a numerical array
num = function(ele){
  as.numeric(str_extract_all(ele, "\\d+\\.*\\d*")[[1]])
}

# isolate duplicates; take out patients with multiple visits
MultVisit = data[duplicated(data$NEW_MASKED_MRN),]
noMult = data[-which(duplicated(data$NEW_MASKED_MRN)),]

# Covid positive or Negative? store in cc column
noMult$cc = c()
noMult$cc[noMult$COVID_RESULT == "DETECTED" | noMult$COVID_RESULT == "PRESUMPTIVE POSITIVE"] = "Detected"
noMult$cc[noMult$COVID_RESULT == "NOT DETECTED" | noMult$COVID_RESULT == "INVALID"] = "Not Detected"

# Categorize the BMI based on its numerical value
noMult$BMIc = noMult$BMI
noMult$BMIc[noMult$BMI < 18.5 & !is.na(noMult$BMI)]= 'Underweight'
noMult$BMIc[noMult$BMI >= 18.5 & noMult$BMI <= 25 & !is.na(noMult$BMI)] = 'NL'
noMult$BMIc[noMult$BMI > 25 & noMult$BMI <= 30 & !is.na(noMult$BMI)] = 'Overweight'
noMult$BMIc[noMult$BMI > 30 & !is.na(noMult$BMI)] = 'Obese'
noMult$BMIc = factor(noMult$BMIc)

# Categorize Smoking by combining Passive, quit and yes into current or previous smoker
noMult$smoke = c()
noMult$smoke[noMult$SMOKING_STATUS == "" | noMult$SMOKING_STATUS == "NOT ASKED"] = NA
noMult$smoke[noMult$SMOKING_STATUS == "PASSIVE" | noMult$SMOKING_STATUS == "QUIT" | noMult$SMOKING_STATUS == "YES"] = "CurPre"
noMult$smoke[noMult$SMOKING_STATUS == "NEVER"] = "Never"

# Categorize Respiratory Rate 
noMult$RR = c()
noMult$RR[is.na(noMult$RESPIRATORY_RATE)] = NA
noMult$RR[!is.na(noMult$RESPIRATORY_RATE) & noMult$RESPIRATORY_RATE <= 24] = "twoFour"
noMult$RR[!is.na(noMult$RESPIRATORY_RATE) & noMult$RESPIRATORY_RATE > 24 & noMult$RESPIRATORY_RATE <= 30] = "lthirty"
noMult$RR[!is.na(noMult$RESPIRATORY_RATE) & noMult$RESPIRATORY_RATE > 30] = "mthirty"

# Categorize O2 Saturation 
noMult$O2 = c()
noMult$O2[is.na(noMult$O2_SAT)] = NA
noMult$O2[!is.na(noMult$O2_SAT) & noMult$O2_SAT <= 93] = "lninthree"
noMult$O2[!is.na(noMult$O2_SAT) & noMult$O2_SAT > 93] = "mninthree"

# Categorize Heart Rate 
noMult$HR = c()
noMult$HR[is.na(noMult$HEART_RATE)] = NA
noMult$HR[!is.na(noMult$HEART_RATE) & noMult$HEART_RATE <= 100] = "l100"
noMult$HR[!is.na(noMult$HEART_RATE) & noMult$HEART_RATE > 100] = "m100"

# add first recorded value for each lab
firstLab = read.csv("firstVal.csv")
noMultLab = cbind(firstLab,noMult)

# Categorize White Blood Cell Count 
noMultLab$WBCc = c()
noMultLab$WBCc[is.na(noMultLab$WBC)] = NA
noMultLab$WBCc[!is.na(noMultLab$WBC) & noMultLab$WBC < 4] = "Low"
noMultLab$WBCc[!is.na(noMultLab$WBC) & noMultLab$WBC >=  4 & noMultLab$WBC <= 11] = "NL"
noMultLab$WBCc[!is.na(noMultLab$WBC) & noMultLab$WBC > 11] = "High"

# Categorize Creatinine
noMultLab$Creat = c()
noMultLab$Creat[is.na(noMultLab$SERUM.CREATININE)] = NA
noMultLab$Creat[!is.na(noMultLab$SERUM.CREATININE) & noMultLab$SERUM.CREATININE >1.2] = "High"
noMultLab$Creat[!is.na(noMultLab$SERUM.CREATININE) & noMultLab$SERUM.CREATININE <= 1.2] = "NL"

# Categorize Platelet
noMultLab$Pltc = c()
noMultLab$Pltc[is.na(noMultLab$PLATELET)] = NA
noMultLab$Pltc[!is.na(noMultLab$PLATELET) & noMultLab$PLATELET < 150] = "Low"
noMultLab$Pltc[!is.na(noMultLab$PLATELET) & noMultLab$PLATELET >=  150 & noMultLab$WBC <= 450] = "NL"
noMultLab$Pltc[!is.na(noMultLab$PLATELET) & noMultLab$PLATELET > 450] = "High"

# Categorize AST
noMultLab$ASTc = c()
noMultLab$ASTc[is.na(noMultLab$AST)] = NA
noMultLab$ASTc[!is.na(noMultLab$AST) & noMultLab$AST < 36] = "Low"
noMultLab$ASTc[!is.na(noMultLab$AST) & noMultLab$AST >=  36 & noMultLab$WBC <= 72] = "NL"
noMultLab$ASTc[!is.na(noMultLab$AST) & noMultLab$AST > 72] = "High"

# Categorize ALT
noMultLab$ALTc = c()
noMultLab$ALTc[is.na(noMultLab$ALT)] = NA
noMultLab$ALTc[!is.na(noMultLab$ALT) & noMultLab$ALT < 46] = "Low"
noMultLab$ALTc[!is.na(noMultLab$ALT) & noMultLab$ALT >=  46 & noMultLab$WBC <= 92] = "NL"
noMultLab$ALTc[!is.na(noMultLab$ALT) & noMultLab$ALT > 92] = "High"

# Categorize Potassium
noMultLab$Potc = c()
noMultLab$Potc[is.na(noMultLab$POTASSIUM)] = NA
noMultLab$Potc[!is.na(noMultLab$POTASSIUM) & noMultLab$POTASSIUM <= 5.2] = "NL"
noMultLab$Potc[!is.na(noMultLab$POTASSIUM) & noMultLab$POTASSIUM > 5.2] = "High"

# Categorize Troponin
noMultLab$tropc = noMultLab$TROPONIN.I
noMultLab$tropc[!is.na(noMultLab$TROPONIN.I) & 0.031 < noMultLab$TROPONIN.I  & noMultLab$TROPONIN.I <=2] = 'Mild'
noMultLab$tropc[!is.na(noMultLab$TROPONIN.I) & noMultLab$TROPONIN.I > 2] = 'High'
noMultLab$tropc[!is.na(noMultLab$TROPONIN.I) & noMultLab$TROPONIN.I <= 0.031] = 'NL'

# Categorize Fibrinogen
noMultLab$fibrc = noMultLab$FIBRINOGEN
noMultLab$fibrc[!is.na(noMultLab$FIBRINOGEN) & 516 < noMultLab$FIBRINOGEN ] = 'High'
noMultLab$fibrc[!is.na(noMultLab$FIBRINOGEN) & noMultLab$FIBRINOGEN <= 516] = 'NL'

# Categorize Ferritin
noMultLab$ferc = noMultLab$FERRITIN
noMultLab$ferc[!is.na(noMultLab$FERRITIN) & 400 < noMultLab$FERRITIN ] = 'High'
noMultLab$ferc[!is.na(noMultLab$FERRITIN) & noMultLab$FERRITIN <= 400] = 'NL'

# Categorize Albumin
noMultLab$albc = noMultLab$ALBUMIN
noMultLab$albc[!is.na(noMultLab$ALBUMIN) &  noMultLab$ALBUMIN < 3.5] = 'Low'
noMultLab$albc[!is.na(noMultLab$ALBUMIN) & noMultLab$ALBUMIN >= 3.5] = 'NL'

# Save the processed Data
write.csv(noMultLab,file="processedData.csv")
noMultLab = read.csv("processedData.csv")


# Baseline Tables ------------------------------------------------------------------

# Factor the categorical variables we're interested in
table = noMultLab
table$death <- factor(table$death, levels=c(0,1), labels=c("Survivor", "Non-Survivor"))
table$SEX <- factor(table$SEX, levels=c("FEMALE","MALE"), labels=c("Female", "Male"))
table$ICU <- factor(table$ICU, levels=c(0,1), labels=c("no ICU", "ICU"))
table$ASTHMA <- factor(table$ASTHMA, levels=c(0,1), labels=c("No","Yes"))
table$COPD <- factor(table$COPD, levels=c(0,1), labels=c("No","Yes"))
table$HTN <- factor(table$HTN, levels=c(0,1), labels=c("No","Yes"))
table$OBSTRUCTIVE_SLEEP_APNEA <- factor(table$OBSTRUCTIVE_SLEEP_APNEA, levels=c(0,1), labels=c("No","Yes"))
table$DIABETES <- factor(table$DIABETES, levels=c(0,1), labels=c("No","Yes"))
table$CHRONIC_KIDNEY_DISEASE <- factor(table$CHRONIC_KIDNEY_DISEASE, levels=c(0,1), labels=c("No","Yes"))
table$CANCER_FLAG <- factor(table$CANCER_FLAG, levels=c(0,1), labels=c("No","Yes"))
table$CORONARY_ARTERY_DISEASE <- factor(table$CORONARY_ARTERY_DISEASE, levels=c(0,1), labels=c("No","Yes"))
table$ATRIAL_FIBRILLATION <- factor(table$ATRIAL_FIBRILLATION, levels=c(0,1), labels=c("No","Yes"))
table$HEART_FAILURE <- factor(table$HEART_FAILURE, levels=c(0,1), labels=c("No","Yes"))
table$CHRONIC_VIRAL_HEPATITIS <- factor(table$CHRONIC_VIRAL_HEPATITIS, levels=c(0,1), labels=c("No","Yes"))
table$ALCOHOLIC_NONALCOHOLIC_LIVER_DISEASE <- factor(table$ALCOHOLIC_NONALCOHOLIC_LIVER_DISEASE, levels=c(0,1), labels=c("No","Yes"))
table$ARDS <- factor(table$ARDS, levels=c(0,1), labels=c("No","Yes"))
table$ACUTE_KIDNEY_INJURY <- factor(table$ACUTE_KIDNEY_INJURY, levels=c(0,1), labels=c("No","Yes"))
table$ACUTE_VENOUS_THROMBOEMBOLISM <- factor(table$ACUTE_VENOUS_THROMBOEMBOLISM, levels=c(0,1), labels=c("No","Yes"))
table$CEREBRAL_INFARCTION <- factor(table$CEREBRAL_INFARCTION, levels=c(0,1), labels=c("No","Yes"))
table$INTRACEREBRAL_HEMORRHAGE <- factor(table$INTRACEREBRAL_HEMORRHAGE, levels=c(0,1), labels=c("No","Yes"))
table$ACUTE_MI <- factor(table$ACUTE_MI, levels=c(0,1), labels=c("No","Yes"))
table$TOCILIZUMAB <- factor(table$TOCILIZUMAB, levels=c(0,1), labels=c("No","Yes"))
table$REMDESIVIR <- factor(table$REMDESIVIR, levels=c(0,1), labels=c("No","Yes"))
table$SARILUMAB <- factor(table$SARILUMAB, levels=c(0,1), labels=c("No","Yes"))
table$HYDROXYCHLOROQUINE <- factor(table$HYDROXYCHLOROQUINE, levels=c(0,1), labels=c("No","Yes"))
table$ANAKINRA <- factor(table$ANAKINRA, levels=c(0,1), labels=c("No","Yes"))
table$AZITHROMYCIN <- factor(table$AZITHROMYCIN, levels=c(0,1), labels=c("No","Yes"))

# how to display continuous vars
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Median (IQR)"=sprintf("%s (%s - %s)", MEDIAN, Q1, Q3)))
}

# how to display categorical vars
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

# Demographics + Patient History
table1(~ SEX + AGE + BMI + BMIc + FACILITY+ RACE_ETHNICITY_COMBINED + smoke + ICU + DISCHARGE_DAYS_SINCE_ENCOUNTER + ASTHMA + COPD + HTN+OBSTRUCTIVE_SLEEP_APNEA +DIABETES+CHRONIC_KIDNEY_DISEASE+CANCER_FLAG+CORONARY_ARTERY_DISEASE+ATRIAL_FIBRILLATION+HEART_FAILURE+CHRONIC_VIRAL_HEPATITIS+ALCOHOLIC_NONALCOHOLIC_LIVER_DISEASE + ARDS+ACUTE_KIDNEY_INJURY+ACUTE_VENOUS_THROMBOEMBOLISM+CEREBRAL_INFARCTION+INTRACEREBRAL_HEMORRHAGE+ACUTE_MI+BLOOD_TYPE| cc*death , data=table, topclass = "Rtable1-zebra", overall = "All Patients", render.continuous = my.render.cont, render.categorical = my.render.cat)

# Vitals
table1(~ TEMPERATURE+HEART_RATE+RESPIRATORY_RATE+O2_SAT+SYSTOLIC_BP+DIASTOLIC_BP+RR+O2+HR| cc*death , data=table, topclass = "Rtable1-zebra", overall = "All Patients", render.continuous = my.render.cont, render.categorical = my.render.cat)

#CBC
table1(~ RBC.COUNT+HEMOGLOBIN+HEMATOCRIT+WBC+WBCc+NEUTROPHIL..+EOSINOPHIL..+LYMPHOCYTE..+HCO3.VENOUS+MEAN.PLATELET.VOLUME..MPV.+PLATELET+MCHC+MCV+MCH+MONOCYTE..+ANION_GAP+PROCALCITONIN| cc*death , data=table, topclass = "Rtable1-zebra", overall = "All Patients", render.continuous = my.render.cont, render.categorical = my.render.cat)

# Arteries / Venous
table1(~ PH.ARTERIAL+PO2.VENOUS+PH.CAPILLARY+PCO2.ARTERIAL+HCO3.CAPILLARY+BASOPHIL..+PH.VENOUS+O2.SATURATION.VENOUS+O2.SATURATION.ARTERIAL+HCO3.VENOUS+O2.SATURATION.CAPILLARY+BASE.EXCESS.ARTERIAL+Pltc+PC02.VENOUS+PCO2.CAPILLARY+BASE.EXCESS.VENOUS+HCO3.ARTERIAL+PO2.ARTERIAL+LACTATE_ARTERIAL| cc*death , data=table, topclass = "Rtable1-zebra", overall = "All Patients", render.continuous = my.render.cont, render.categorical = my.render.cat)

# Kidney
table1(~ URIC.ACID+GLUCOSE+SODIUM+CHLORIDE+POTASSIUM+Potc+CALCIUM+BUN+SERUM.CREATININE+Creat+EGFR| cc*death , data=table, topclass = "Rtable1-zebra", overall = "All Patients", render.continuous = my.render.cont, render.categorical = my.render.cat)

# Liver
table1(~ GLUCOSE +ALBUMIN+ albc+ALKALINE.PHOSPHATASE+ALT+CREATINE_KINASE+AST+LDH+ASTc+ALTc+TOTAL.BILIRUBIN| cc*death , data=table, topclass = "Rtable1-zebra", overall = "All Patients", render.continuous = my.render.cont, render.categorical = my.render.cat)

# Inflammatory
table1(~ C.REACTIVE.PROTEIN +FERRITIN + ferc+ESR+INTERLEUKIN.8+INTERLEUKIN.1.BETA+TNF.ALPHA+ INTERLEUKIN.6| cc*death , data=table, topclass = "Rtable1-zebra", overall = "All Patients", render.continuous = my.render.cont, render.categorical = my.render.cat)

# Coagulation State
table1(~ D.DIMER + BETA.2.GLYCOPROTEIN.I.ANTIBODY +PTT+PROTHROMBIN.TIME+FIBRINOGEN+fibrc+INR+PHOSPHOLIPID.ANTIBODY+ANTICARDIOLIPIN.ANTIBODY| cc*death , data=table, topclass = "Rtable1-zebra", overall = "All Patients", render.continuous = my.render.cont, render.categorical = my.render.cat)

# Cardiac
table1(~ TROPONIN.I +tropc +CREATINE_KINASE_MB +BRAIN.NATRIURETIC.PROTEIN| cc*death , data=table, topclass = "Rtable1-zebra", overall = "All Patients", render.continuous = my.render.cont, render.categorical = my.render.cat)

# Medication
table1(~ TOCILIZUMAB + REMDESIVIR + SARILUMAB + HYDROXYCHLOROQUINE + ANAKINRA + AZITHROMYCIN | cc*death , data=table, topclass = "Rtable1-zebra", overall = "All Patients", render.continuous = my.render.cont, render.categorical = my.render.cat)


# Exploratory Plots -------------------------------------------------------

# Function Used For plotting trajectories 
# Input: data frame of choice (Temperature), and name of y axis, d specifies whether to plot just covid positive (if d=1)
# Output: plot
plotTraj = function(tmp,name, d = 1){
  # convert the data into a dataframe, with days, calculating the mean and standard error for 4 categories of covid + and - as well as death status 
  df = data.frame(Day = rep(1:dim(tmp)[2],each = 4), mean = rep(NA,dim(tmp)[2]*4), survived = rep(c("Covid + Survivor","Covid + non-Survivor","Covid - Survivor","Covid - non-survivor"), dim(tmp)[2]), se = rep(NA,dim(tmp)[2]*4))
  
  j = 1;
  for(i in 1:dim(tmp)[2]){ #calculated the mean and standard error for each of 4 categories
    df$mean[j] = mean(tmp[noMultLab$death == 0 & noMultLab$cc == "Detected",i],na.rm = T)
    df$mean[j + 1] = mean(tmp[noMultLab$death == 1 & noMultLab$cc == "Detected",i],na.rm = T)
    df$mean[j + 2] = mean(tmp[noMultLab$death == 0 & noMultLab$cc == "Not Detected",i],na.rm = T)
    df$mean[j + 3] = mean(tmp[noMultLab$death == 1 & noMultLab$cc == "Not Detected",i],na.rm = T)
    df$se[j] = sd(tmp[noMultLab$death == 0 & noMultLab$cc == "Detected",i],na.rm = T)/sqrt(length(tmp[noMultLab$death == 0 & noMultLab$cc == "Detected",i][!is.na(tmp[noMultLab$death == 0 & noMultLab$cc == "Detected",i])]))
    df$se[j+1] = sd(tmp[noMultLab$death == 1 & noMultLab$cc == "Detected",i],na.rm = T)/sqrt(length(tmp[noMultLab$death == 1 & noMultLab$cc == "Detected",i][!is.na(tmp[noMultLab$death == 0 & noMultLab$cc == "Detected",i])]))
    df$se[j+2] = sd(tmp[noMultLab$death == 0 & noMultLab$cc == "Not Detected",i],na.rm = T)/sqrt(length(tmp[noMultLab$death == 0 & noMultLab$cc == "Not Detected",i][!is.na(tmp[noMultLab$death == 0 & noMultLab$cc == "Not Detected",i])]))
    df$se[j+3] = sd(tmp[noMultLab$death == 1 & noMultLab$cc == "Not Detected",i],na.rm = T)/sqrt(length(tmp[noMultLab$death == 1 & noMultLab$cc == "Not Detected",i][!is.na(tmp[noMultLab$death == 0 & noMultLab$cc == "Not Detected",i])]))
    j = j+4
  }
  
  if(d == 1){
    #df[df$survived == "Covid + Survivor" | df$survived == "Covid + non-Survivor",]
    qplot(x =  Day, y  =  mean, color = survived, data =  df[df$survived == "Covid + Survivor" | df$survived == "Covid + non-Survivor",]) + 
      theme_bw() +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.15)) + geom_line(linetype="dotted") + ylab(name) +theme(text = element_text(size=25))
  } else { # only plot covid + trajectories
    qplot(x =  Day, y  =  mean, color = survived, data =  df) + 
      theme_bw() +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.15)) + geom_line(linetype="dotted") + ylab(name)
  }
}

# calculates daily averages
num2 = function(ele){
  mean(as.numeric(str_extract_all(ele, "\\d+\\.*\\d*")[[1]]))
}

# calculate daily avergaes
averages = data[-which(duplicated(data$NEW_MASKED_MRN)),]
for(i in 86:dim(data)[2]-1){
  averages[,i] = sapply(data[-which(duplicated(data$NEW_MASKED_MRN)),i], num2)
  print(i)
}
write.csv(averages,file = "averages.csv")
averages = read.csv("averages.csv")

# Plot Vitals Trajectories
tmp = averages[,grep("TEMP_MAX__0", colnames(averages)):grep("TEMP_MAX__29", colnames(averages))]
p1 = plotTraj(tmp,"Temperature (F)")
HR = averages[,grep("HEART_RATE_MAX__0", colnames(averages)):grep("HEART_RATE_MAX__29", colnames(averages))]
p2 = plotTraj(HR,"Heart Rate (bpm)")
RR = averages[,grep("RESP_RATE_MAX__0", colnames(averages)):grep("RESP_RATE_MAX__29", colnames(averages))]
p3 = plotTraj(RR,"Respiratory Rate")
Sys = averages[,grep("SYSTOL_BP_MAX__0", colnames(averages)):grep("SYSTOL_BP_MAX__29", colnames(averages))]
p4 = plotTraj(Sys,"Systolic BP (mmHg)")
Dia = averages[,grep("DIASTOL_BP_MAX__0", colnames(averages)):grep("DIASTOL_BP_MAX__29", colnames(averages))]
p5 = plotTraj(Dia,"Diastolic BP (mmHg)")
O2 = averages[,grep("O2SAT_MIN__0", colnames(averages)):grep("O2SAT_MIN__29", colnames(averages))]
p6 = plotTraj(O2,"O2 Saturation (%)")
# arrange the above plots in 1 plot vitals
ggarrange(p1, p2, p3, p4, p5, p6,
          #labels = c("A", "B", "C","D", "E","F"), # labels of each subplot
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend


# Plot CBC Trajectories
wbc = averages[,grep("WBC__0", colnames(averages)):grep("WBC__29", colnames(averages))]
p1c = plotTraj(wbc,"WBC Count (K/uL)")
rbc = averages[,grep("RBC.COUNT__0", colnames(averages)):grep("RBC.COUNT__29", colnames(averages))]
p2c = plotTraj(rbc,"RBC Count (m/uL)")
plt = averages[,grep("PLATELET__0", colnames(averages)):grep("PLATELET__29", colnames(averages))]
p3c = plotTraj(plt,"Platelet (K/uL)")
mpv = averages[,grep("MEAN.PLATELET.VOLUME..MPV.__0", colnames(averages)):grep("MEAN.PLATELET.VOLUME..MPV.__29", colnames(averages))]
p4c = plotTraj(mpv,"MPV (fL)")
neu = averages[,grep("NEUTROPHIL..__0", colnames(averages)):grep("NEUTROPHIL..__29", colnames(averages))]
p5c = plotTraj(neu,"Neutrophil (K/uL)")
eos = averages[,grep("EOSINOPHIL..__0", colnames(averages)):grep("EOSINOPHIL..__29", colnames(averages))]
p6c = plotTraj(eos,"Eosinophil (K/uL)")
mono = averages[,grep("MONOCYTE..__0", colnames(averages)):grep("MONOCYTE..__29", colnames(averages))]
p7c = plotTraj(mono,"Monocyte (K/uL)")
bas = averages[,grep("BASOPHIL..__0", colnames(averages)):grep("BASOPHIL..__29", colnames(averages))]
p8c = plotTraj(bas,"Basophil (K/uL)")
lym = averages[,grep("LYMPHOCYTE..__0", colnames(averages)):grep("LYMPHOCYTE..__29", colnames(averages))]
p9c = plotTraj(lym,"Lymphocyte (K/uL)")
hemo = averages[,grep("HEMOGLOBIN__0", colnames(averages)):grep("HEMOGLOBIN__29", colnames(averages))]
p10c = plotTraj(hemo,"Hemoglobin (g/dL)")
mch = averages[,grep("MCH__0", colnames(averages)):grep("MCH__29", colnames(averages))]
p11c = plotTraj(mch,"MCH (pg)")
mchc = averages[,grep("MCHC__0", colnames(averages)):grep("MCHC__29", colnames(averages))]
p12c = plotTraj(mchc,"MCHC (g/dL)")
mcv = averages[,grep("MCV__0", colnames(averages)):grep("MCV__29", colnames(averages))]
p13c = plotTraj(mcv,"MCV (fL)")
hemat = averages[,grep("HEMATOCRIT__0", colnames(averages)):grep("HEMATOCRIT__29", colnames(averages))]
p14c = plotTraj(hemat,"Hematocrit (%)")
an = averages[,grep("ANION_GAP__0", colnames(averages)):grep("ANION_GAP__29", colnames(averages))]
p4b = plotTraj(an,"Anion Gap (mmol/L)")
procal = averages[,grep("PROCALCITONIN__0", colnames(averages)):grep("PROCALCITONIN__29", colnames(averages))]
p16o = plotTraj(procal,"Procalcitonin (ng/mL)")
# Plot CBC results
ggarrange(p1c, p5c, p6c, p7c, p8c, p9c, 
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend

ggarrange(p2c,p10c,p14c,p12c,p13c,p4b,
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend

ggarrange(p3c,p4c,p16o,
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend


# Liver
alb = averages[,grep("ALBUMIN__0", colnames(averages)):grep("ALBUMIN__29", colnames(averages))]
p2d = plotTraj(alb,"Albumin (g/dL)")
alt = averages[,grep("ALT__0", colnames(averages)):grep("ALT__29", colnames(averages))]
p3d = plotTraj(alt,"ALT (U/L)")
ast = averages[,grep("AST__0", colnames(averages)):grep("AST__29", colnames(averages))]
p4d = plotTraj(ast,"AST (U/L)")
bili = averages[,grep("TOTAL.BILIRUBIN__0", colnames(averages)):grep("TOTAL.BILIRUBIN__29", colnames(averages))]
p5d = plotTraj(bili,"Total Bilirubin (mg/dL)")
alp = averages[,grep("ALKALINE.PHOSPHATASE__0", colnames(averages)):grep("ALKALINE.PHOSPHATASE__29", colnames(averages))]
p6d = plotTraj(alp, "ALP (U/L)")
ldh = averages[,grep("LDH__0", colnames(averages)):grep("LDH__29", colnames(averages))]
p4o = plotTraj(ldh,"LDH (U/L)")
glu = averages[,grep("GLUCOSE__0", colnames(averages)):grep("GLUCOSE__29", colnames(averages))]
p1b = plotTraj(glu,"Glucose (mg/dL)")
crtK = averages[,grep("CREATINE_KINASE__0", colnames(averages)):grep("CREATINE_KINASE__29", colnames(averages))]
p9b = plotTraj(crtK,"Creatinine Kinase (U/L)")
# Plot Liver results
ggarrange(p3d, p4d, p2d,
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend

ggarrange(p5d, p4o,p1b,p9b,
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend


# Inflammatory Markers
crp = averages[,grep("C.REACTIVE.PROTEIN__0", colnames(averages)):grep("C.REACTIVE.PROTEIN__29", colnames(averages))]
p1o = plotTraj(crp,"CRP (mg/L)")
fer = averages[,grep("FERRITIN__0", colnames(averages)):grep("FERRITIN__29", colnames(averages))]
p5o = plotTraj(fer,"Ferritin (ng/mL)")
esr = averages[,grep("ESR__0", colnames(averages)):grep("ESR__29", colnames(averages))]
p6o = plotTraj(esr,"ESR (mm/hr)")
il6 = averages[,grep("INTERLEUKIN.6__0", colnames(averages)):grep("INTERLEUKIN.6__29", colnames(averages))]
p11o = plotTraj(il6,"IL 6 (pg/mL)")
il1 = averages[,grep("INTERLEUKIN.1.BETA__0", colnames(averages)):grep("INTERLEUKIN.1.BETA__29", colnames(averages))]
p12o = plotTraj(il1,"IL 1 Beta pg/mL")
il8 = averages[,grep("INTERLEUKIN.8__0", colnames(averages)):grep("INTERLEUKIN.8__29", colnames(averages))]
p13o = plotTraj(il8,"IL 8 pg/mL")
tnf = averages[,grep("TNF.ALPHA__0", colnames(averages)):grep("TNF.ALPHA__29", colnames(averages))]
p20o = plotTraj(tnf,"TNF Alpha (pg/mL)")
# Plot Inflammatory results
ggarrange(p1o, p5o, p6o, p11o, p20o,
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend

# Coagulation
dd = averages[,grep("D.DIMER__0", colnames(averages)):grep("D.DIMER__29", colnames(averages))]
p2o = plotTraj(dd,"DDimer (ug/mL)")
b2g = averages[,grep("BETA.2.GLYCOPROTEIN.I.ANTIBODY__0", colnames(averages)):grep("BETA.2.GLYCOPROTEIN.I.ANTIBODY__8", colnames(averages))]
p7o = plotTraj(b2g,"Beta 2 glyco-protein")
ptt = averages[,grep("PTT__0", colnames(averages)):grep("PTT__29", colnames(averages))]
p8o = plotTraj(ptt,"PTT (s)")
pt = averages[,grep("PROTHROMBIN.TIME__0", colnames(averages)):grep("PROTHROMBIN.TIME__29", colnames(averages))]
p9o = plotTraj(pt,"Prothrombin Time (s)")
fib = averages[,grep("FIBRINOGEN__0", colnames(averages)):grep("FIBRINOGEN__29", colnames(averages))]
p10o = plotTraj(fib,"Fibrinogen (mg/dL)")
inr = averages[,grep("INR__0", colnames(averages)):grep("INR__29", colnames(averages))]
p15o = plotTraj(inr,"INR")
phph = averages[,grep("PHOSPHOLIPID.ANTIBODY__0", colnames(averages)):grep("PHOSPHOLIPID.ANTIBODY__15", colnames(averages))]
p17o = plotTraj(phph,"Phospholipid Antibody")
ant = averages[,grep("ANTICARDIOLIPIN.ANTIBODY__0", colnames(averages)):grep("ANTICARDIOLIPIN.ANTIBODY__8", colnames(averages))]
p18o = plotTraj(ant,"Anticardiolipin Antibody")
# Plot Coagulation Results
ggarrange(p2o,p8o, p9o,p15o, p10o,  
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend

# Kidney
eg = averages[,grep("EGFR__0", colnames(averages)):grep("EGFR__29", colnames(averages))]
p3o = plotTraj(eg,"EGFR (mL/min)")
ura = averages[,grep("URIC.ACID__0", colnames(averages)):grep("URIC.ACID__29", colnames(averages))]
p19o = plotTraj(ura,"Uric Acid (mg/dL)")
sod = averages[,grep("SODIUM__0", colnames(averages)):grep("SODIUM__29", colnames(averages))]
p2b = plotTraj(sod,"Sodium (mmol/L)")
ch = averages[,grep("CHLORIDE__0", colnames(averages)):grep("CHLORIDE__29", colnames(averages))]
p3b = plotTraj(ch,"Chloride (mmol/L)")
cal = averages[,grep("CALCIUM__0", colnames(averages)):grep("CALCIUM__29", colnames(averages))]
p5b = plotTraj(cal,"Calcium (mg/dL)")
pot = averages[,grep("POTASSIUM__0", colnames(averages)):grep("POTASSIUM__29", colnames(averages))]
p6b = plotTraj(pot,"Potassium (mmol/L)")
bun = averages[,grep("BUN__0", colnames(averages)):grep("BUN__29", colnames(averages))]
p7b = plotTraj(bun,"BUN (mg/dL)")
crt = averages[,grep("SERUM.CREATININE__0", colnames(averages)):grep("SERUM.CREATININE__29", colnames(averages))]
p8b = plotTraj(crt,"Creatinine (mg/dL)")
# Plot Kidney Results
ggarrange(p8b, p3o,p7b, p19o, p2b, p3b, 
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend

ggarrange(p5b, p6b, 
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend

# Cardiac
bnp = averages[,grep("BRAIN.NATRIURETIC.PROTEIN__0", colnames(averages)):grep("BRAIN.NATRIURETIC.PROTEIN__29", colnames(averages))]
p14o = plotTraj(bnp,"BNP (pg/mL)")
trop = averages[,grep("TROPONIN.I__0", colnames(averages)):grep("TROPONIN.I__29", colnames(averages))]
p21o = plotTraj(trop,"Troponin I (ng/mL)")
crtKmb = averages[,grep("CREATINE_KINASE_MB__0", colnames(averages)):grep("CREATINE_KINASE_MB__29", colnames(averages))]
p10b = plotTraj(crtKmb,"CK-MB (U/L)")
# Plot Cardiac Results
ggarrange(p14o, p21o, p10b,
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend


#arteries/veins
pha = averages[,grep("PH.ARTERIAL__0", colnames(averages)):grep("PH.ARTERIAL__29", colnames(averages))]
p1av = plotTraj(pha,"pH Arterial")
phv = averages[,grep("PH.VENOUS__0", colnames(averages)):grep("PH.VENOUS__29", colnames(averages))]
p2av = plotTraj(phv,"pH Venous")
phc = averages[,grep("PH.CAPILLARY__0", colnames(averages)):grep("PH.CAPILLARY__19", colnames(averages))]
p3av = plotTraj(phc,"pH Capillary")
ba = averages[,grep("BASE.EXCESS.ARTERIAL__0", colnames(averages)):grep("BASE.EXCESS.ARTERIAL__29", colnames(averages))]
p4av = plotTraj(ba,"BE Arterial (mmol/L)")
bv = averages[,grep("BASE.EXCESS.VENOUS__0", colnames(averages)):grep("BASE.EXCESS.VENOUS__29", colnames(averages))]
p5av = plotTraj(bv,"BE Venous (mmol/L)")
la = averages[,grep("LACTATE_ARTERIAL__0", colnames(averages)):grep("LACTATE_ARTERIAL__29", colnames(averages))]
p6av = plotTraj(la,"Lactate Arterial (mmol/L)")
pca = averages[,grep("PCO2.ARTERIAL__0", colnames(averages)):grep("PCO2.ARTERIAL__29", colnames(averages))]
p7av = plotTraj(pca,"PCO2 Arterial (mmHg)")
pcv = averages[,grep("PC02.VENOUS__0", colnames(averages)):grep("PC02.VENOUS__29", colnames(averages))]
p8av = plotTraj(pcv,"PCO2 Venous (mmHg)")
pcc = averages[,grep("PCO2.CAPILLARY__0", colnames(averages)):grep("PCO2.CAPILLARY__19", colnames(averages))]
p9av = plotTraj(pcc,"PCO2 Capillary (mmHg)")
poa = averages[,grep("PO2.ARTERIAL__0", colnames(averages)):grep("PO2.ARTERIAL__29", colnames(averages))]
p10av = plotTraj(poa,"PO2 Arterial (mmHg)")
pov = averages[,grep("PO2.VENOUS__0", colnames(averages)):grep("PO2.VENOUS__29", colnames(averages))]
p11av = plotTraj(pov,"PO2 Venous (mmHg)")
oa = averages[,grep("O2.SATURATION.ARTERIAL__0", colnames(averages)):grep("O2.SATURATION.ARTERIAL__29", colnames(averages))]
p12av = plotTraj(oa,"O2 Sat Arterial (%)")
ov = averages[,grep("O2.SATURATION.VENOUS__0", colnames(averages)):grep("O2.SATURATION.VENOUS__29", colnames(averages))]
p13av = plotTraj(ov,"O2 Sat Venous (%)")
oc = averages[,grep("O2.SATURATION.CAPILLARY__0", colnames(averages)):grep("O2.SATURATION.CAPILLARY__9", colnames(averages))]
p14av = plotTraj(oc,"O2 Sat Capillary (%)")
hca = averages[,grep("HCO3.ARTERIAL__0", colnames(averages)):grep("HCO3.ARTERIAL__29", colnames(averages))]
p15av = plotTraj(hca,"HCO3 Arterial (meq/L)")
hcv = averages[,grep("HCO3.VENOUS__0", colnames(averages)):grep("HCO3.VENOUS__29", colnames(averages))]
p16av = plotTraj(hcv,"HCO3 Venous (meq/L)")
hcc = averages[,grep("HCO3.CAPILLARY__0", colnames(averages)):grep("HCO3.CAPILLARY__19", colnames(averages))]
p17av = plotTraj(hcc,"HCO3 Capillary (meq/L)")
# Plot Arteries/veins Trajectories
ggarrange(p1av, p2av, p4av, 
          #labels = c("A", "B", "C","D", "E","F"), # labels of each subplot
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend

ggarrange(p5av, p6av,p7av, p8av,p10av,p11av,
          #labels = c("A", "B", "C","D", "E","F"), # labels of each subplot
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend

ggarrange(p12av,p13av,p15av,p16av,
          #labels = c("A", "B", "C","D", "E","F"), # labels of each subplot
          ncol = 3, nrow = 2,
          legend = F, # legend location
          common.legend = TRUE) # include only 1 legend


### Length of Stay, read it from the generated file
LoS = read.csv("los.csv")
LoS$cc = noMultLab$cc
LoS$death = factor(noMultLab$death)
# categorize based on both covid detected/notDetected and surviving status 
LoS$ccD = paste(noMultLab$cc, noMultLab$death)
LoS$ccD[is.na(LoS$cc) | is.na(LoS$death)] = NA
LoS$ccD = factor(LoS$ccD, levels = c("Not Detected 1","Not Detected 0","Detected 1","Detected 0"))
levels(LoS$ccD) = c("Covid - Non-Survivor", "Covid - Survivor","Covid + Non-Survivor", "Covid + Survivor")

# Plot histogram of length of stay in the hospital based on whether covid +/- and survived yes or no
ggplot(LoS[!is.na(LoS$ccD),], aes(x=LoS, fill=ccD)) +
  geom_histogram(binwidth = 2) +
  scale_fill_manual(values=c("green3","palegreen1","orangered1", "lightpink")) +
  theme_minimal() +
  labs(fill = "") + 
  theme(text = element_text(size=15))+
  xlab("Length of Stay") +
  ylab("Frequency")  +
  scale_x_continuous(breaks = c(seq(-1,19,5),seq(29,59,10)),labels = c(seq(0,20,5),seq(30,60,10)))

# Plot histogram of length of stay in the hospital based on whether covid +/- and survived yes or no
# Layout to split the screen
layout(mat = matrix(c(1:13,13),7,1, byrow=TRUE),  height = c(1,1,8))
par(mar=c(0, 3.1, 1.1, 2.1), mai=c(0,0.4,0,0.2))
boxplot(LoS$LoS[LoS$cc == "Detected" & LoS$death == 1], horizontal=TRUE, xaxt="n" , col="salmon" , frame=F)
par(mar=c(0, 3.1, 1.1, 2.1), mai=c(0,0.4,0,0.2))
boxplot(LoS$LoS[LoS$cc == "Detected" & LoS$death == 0] , horizontal=TRUE , xaxt="n" , col="steelblue1" , frame=F)
# Draw 2 histograms for education colored based on whether individual's working
par(mar=c(4.1, 3.1, 1.1, 2.1))
hist(LoS$LoS[LoS$cc == "Detected" & LoS$death == 0], breaks  = c(seq(-1,90,2)), xaxt = 'n' , col=rgb(0.12, 0.565, 1,0.4) , border=F , main="" , xlab="Length of Stay", xlim = c(0,50),cex.lab=1.5, cex.axis=1.25, ylab = "")
axis(1, -1:50+1, 1:52)
hist(LoS$LoS[LoS$cc == "Detected" & LoS$death == 1], breaks  =c(seq(-1,90,2)), xaxt = 'n',  col=rgb(1, 0.55, 0.41,0.8), border=F , main="", xlim = c(0,50), add= TRUE)

layout(mat = matrix(c(1:13,13),7,1, byrow=TRUE),  height = c(1,1,8))
par(mar=c(0, 3.1, 1.1, 2.1), mai=c(0,0.4,0,0.2))
boxplot(LoS$LoS[LoS$cc == "Not Detected" & LoS$death == 1], horizontal=TRUE, xaxt="n" , col="salmon" , frame=F)
par(mar=c(0, 3.1, 1.1, 2.1), mai=c(0,0.4,0,0.2))
boxplot(LoS$LoS[LoS$cc == "Not Detected" & LoS$death == 0] , horizontal=TRUE , xaxt="n" , col="steelblue1" , frame=F)
# Draw 2 histograms for education colored based on whether individual's working
par(mar=c(4.1, 3.1, 1.1, 2.1))
hist(LoS$LoS[LoS$cc == "Not Detected" & LoS$death == 0], breaks  = c(seq(-1,90,2)), xaxt = 'n' , col=rgb(0.12, 0.565, 1,0.4) , border=F , main="" , xlab="Length of Stay", xlim = c(0,50),cex.lab=1.5, cex.axis=1.25, ylab = "")
axis(1, -1:50+1, 1:52)
hist(LoS$LoS[LoS$cc == "Not Detected" & LoS$death == 1], breaks  =c(seq(-1,90,2)), xaxt = 'n',  col=rgb(1, 0.55, 0.41,0.8), border=F , main="", xlim = c(0,50), add= TRUE)


ggplot(LoS[LoS$cc == "Detected" & !is.na(LoS$ccD),], aes(x=LoS, fill=ccD))+ 
  geom_histogram(alpha = 0.4, position = 'identity') +
  scale_fill_manual(values=c("red","blue")) +
  theme_minimal() +
  labs(fill = "") + 
  theme(text = element_text(size=15))+
  xlab("Length of Stay") +
  ylab("Frequency")  +
  scale_x_continuous(breaks = c(seq(-1,19,5),seq(29,59,10)),labels = c(seq(0,20,5),seq(30,60,10)))


# Next try plotting trajectories for those staying in the hospital for the same amount of time
### trajectories for those staying 3 days
p1 = plotTraj(tmp[LoS$LoS == 3,1:3],"Temperature (F)",2)
p2 = plotTraj(O2[LoS$LoS == 3,1:3],"O2 Saturation (%)",2)
p3 = plotTraj(RR[LoS$LoS == 3,1:3],"Respiratory Rate",2)

p4 = plotTraj(wbc[LoS$LoS == 3,1:3],"White Blood Cell Count (K/uL)",2)
p5 = plotTraj(lym[LoS$LoS == 3,1:3],"Lymphocyte (K/uL)",2)
p6 = plotTraj(neu[LoS$LoS == 3,1:3],"Neutrophil (K/uL)",2)
p7 = plotTraj(plt[LoS$LoS == 3,1:3],"Platelet (K/uL)")

p8 = plotTraj(dd[LoS$LoS == 3,1:3],"DDimer (ug/mL)",2)
p9 = plotTraj(fer[LoS$LoS == 3,1:3],"Ferritin (ng/mL)",2)
p10 = plotTraj(crp[LoS$LoS == 3,1:3],"C-Reactive Protein (mg/L)",2) #good results
p14 = plotTraj(ldh[LoS$LoS == 3,1:3],"LDH",2)#good

p11 = plotTraj(trop[LoS$LoS == 3,1:3],"Troponin I",2)
p12 = plotTraj(crt[LoS$LoS == 3,1:3],"Serum Creatinine",2)
p13 = plotTraj(eg[LoS$LoS == 3,1:3],"EGFR (mL/min)",2)

# 3 day length of stay plots
ggarrange(p2, p3, p4, p5, p6,p8,p9,p10,p14,p11,p12,p13,
          ncol = 4, nrow = 3,
          legend = "right", # legend location
          common.legend = TRUE) # include only 1 legend


### trjactories for those staying 7 days
p1 = plotTraj(tmp[LoS$LoS == 7,1:7],"Temperature (F)",2)
p2 = plotTraj(O2[LoS$LoS == 7,1:7],"O2 Saturation (%)",2)
p3 = plotTraj(RR[LoS$LoS == 7,1:7],"Respiratory Rate",2)

p4 = plotTraj(wbc[LoS$LoS == 7,1:7],"White Blood Cell Count (K/uL)",2)
p5 = plotTraj(lym[LoS$LoS == 7,1:7],"Lymphocyte (K/uL)",2)
p6 = plotTraj(neu[LoS$LoS == 7,1:7],"Neutrophil (K/uL)",2)
p7 = plotTraj(plt[LoS$LoS == 7,1:7],"Platelet (K/uL)",2)

p8 = plotTraj(dd[LoS$LoS == 7,1:7],"DDimer (ug/mL)",2)
p9 = plotTraj(fer[LoS$LoS == 7,1:7],"Ferritin (ng/mL)",2)
p10 = plotTraj(crp[LoS$LoS == 7,1:7],"C-Reactive Protein (mg/L)",2)
p11 = plotTraj(ldh[LoS$LoS == 7,1:7],"LDH",2)#good

p12 = plotTraj(trop[LoS$LoS == 7,1:7],"Troponin I",2)
p13 = plotTraj(crt[LoS$LoS == 7,1:7],"Serum Creatinine",2)
p14 = plotTraj(eg[LoS$LoS == 7,1:7],"EGFR (mL/min)",2)

# 7 day plots
ggarrange(p2, p3, p4, p5, p6,p8,p9,p10,p11,p12,p13,p14, 
          ncol = 4, nrow = 3,
          legend = "right", # legend location
          common.legend = TRUE) # include only 1 legend


# Modelling ---------------------------------------------------------------
## Step 1 Generate Observational Units
# let's first isolate the static features
obsStatic = noMultLab[, which(names(noMultLab) %in% c("AGE","SEX","RACE","ETHNICITY","FACILITY","ICU","INPATIENT_NON_ICU","cc","SMOKING_STATUS","ASTHMA","COPD","HTN","OBSTRUCTIVE_SLEEP_APNEA","OBESITY","DIABETES","CHRONIC_KIDNEY_DISEASE","HIV_FLAG","CANCER_FLAG","CORONARY_ARTERY_DISEASE","ATRIAL_FIBRILLATION","HEART_FAILURE","CHRONIC_VIRAL_HEPATITIS","ALCOHOLIC_NONALCOHOLIC_LIVER_DISEASE","ARDS","ACUTE_KIDNEY_INJURY","ACUTE_VENOUS_THROMBOEMBOLISM","CEREBRAL_INFARCTION","INTRACEREBRAL_HEMORRHAGE","ACUTE_MI","BMI","BLOOD_TYPE","death"))]
obsStatic$LOS = LoS$LoS

# the indices of covid positive and negative patients, repeat the observational units generation for both arrays
arrCovidPos = which(obsStatic$cc == "Detected")
arrnoCovidPos = which(obsStatic$cc == "Not Detected") 

# Define the trajectory/non static features to be used
# Vitals variables
tmp2 = noMultLab[,grep("TEMP_MAX__0", colnames(noMultLab)):grep("TEMP_MAX__72", colnames(noMultLab))]
HR2 = noMultLab[,grep("HEART_RATE_MAX__0", colnames(noMultLab)):grep("HEART_RATE_MAX__72", colnames(noMultLab))]
RR2 = noMultLab[,grep("RESP_RATE_MAX__0", colnames(noMultLab)):grep("RESP_RATE_MAX__72", colnames(noMultLab))]
Sys2 = noMultLab[,grep("SYSTOL_BP_MAX__0", colnames(noMultLab)):grep("SYSTOL_BP_MAX__72", colnames(noMultLab))]
Dia2 = noMultLab[,grep("DIASTOL_BP_MAX__0", colnames(noMultLab)):grep("DIASTOL_BP_MAX__72", colnames(noMultLab))]
O22 = noMultLab[,grep("O2SAT_MIN__0", colnames(noMultLab)):grep("O2SAT_MIN__72", colnames(noMultLab))]

# CBC Trajectories
wbc2 = noMultLab[,grep("WBC__0", colnames(noMultLab)):grep("WBC__72", colnames(noMultLab))]
rbc2 = noMultLab[,grep("RBC.COUNT__0", colnames(noMultLab)):grep("RBC.COUNT__72", colnames(noMultLab))]
plt2 = noMultLab[,grep("PLATELET__0", colnames(noMultLab)):grep("PLATELET__72", colnames(noMultLab))]
mpv2 = noMultLab[,grep("MEAN.PLATELET.VOLUME..MPV.__0", colnames(noMultLab)):grep("MEAN.PLATELET.VOLUME..MPV.__72", colnames(noMultLab))]
neu2 = noMultLab[,grep("NEUTROPHIL..__0", colnames(noMultLab)):grep("NEUTROPHIL..__72", colnames(noMultLab))]
eos2 = noMultLab[,grep("EOSINOPHIL..__0", colnames(noMultLab)):grep("EOSINOPHIL..__72", colnames(noMultLab))]
mono2 = noMultLab[,grep("MONOCYTE..__0", colnames(noMultLab)):grep("MONOCYTE..__72", colnames(noMultLab))]
bas2 = noMultLab[,grep("BASOPHIL..__0", colnames(noMultLab)):grep("BASOPHIL..__72", colnames(noMultLab))]
lym2 = noMultLab[,grep("LYMPHOCYTE..__0", colnames(noMultLab)):grep("LYMPHOCYTE..__72", colnames(noMultLab))]
hemo2 = noMultLab[,grep("HEMOGLOBIN__0", colnames(noMultLab)):grep("HEMOGLOBIN__72", colnames(noMultLab))]
mch2 = noMultLab[,grep("MCH__0", colnames(noMultLab)):grep("MCH__72", colnames(noMultLab))]
mchc2 = noMultLab[,grep("MCHC__0", colnames(noMultLab)):grep("MCHC__72", colnames(noMultLab))]
mcv2 = noMultLab[,grep("MCV__0", colnames(noMultLab)):grep("MCV__72", colnames(noMultLab))]
hemat2 = noMultLab[,grep("HEMATOCRIT__0", colnames(noMultLab)):grep("HEMATOCRIT__67", colnames(noMultLab))]
an2 = noMultLab[,grep("ANION_GAP__0", colnames(noMultLab)):grep("ANION_GAP__72", colnames(noMultLab))]
procal2 = noMultLab[,grep("PROCALCITONIN__0", colnames(noMultLab)):grep("PROCALCITONIN__57", colnames(noMultLab))]


# Liver
alb2 = noMultLab[,grep("ALBUMIN__0", colnames(noMultLab)):grep("ALBUMIN__67", colnames(noMultLab))]
alt2 = noMultLab[,grep("ALT__0", colnames(noMultLab)):grep("ALT__67", colnames(noMultLab))]
ast2 = noMultLab[,grep("AST__0", colnames(noMultLab)):grep("AST__67", colnames(noMultLab))]
bili2 = noMultLab[,grep("TOTAL.BILIRUBIN__0", colnames(noMultLab)):grep("TOTAL.BILIRUBIN__67", colnames(noMultLab))]
alp2 = noMultLab[,grep("ALKALINE.PHOSPHATASE__0", colnames(noMultLab)):grep("ALKALINE.PHOSPHATASE__14", colnames(noMultLab))] # many missing
ldh2 = noMultLab[,grep("LDH__0", colnames(noMultLab)):grep("LDH__60", colnames(noMultLab))]
glu2 = noMultLab[,grep("GLUCOSE__0", colnames(noMultLab)):grep("GLUCOSE__72", colnames(noMultLab))]
crtK2 = noMultLab[,grep("CREATINE_KINASE__0", colnames(noMultLab)):grep("CREATINE_KINASE__55", colnames(noMultLab))]


# Inflammatory Markers
crp2 = noMultLab[,grep("C.REACTIVE.PROTEIN__0", colnames(noMultLab)):grep("C.REACTIVE.PROTEIN__56", colnames(noMultLab))]
fer2 = noMultLab[,grep("FERRITIN__0", colnames(noMultLab)):grep("FERRITIN__59", colnames(noMultLab))]
esr2 = noMultLab[,grep("ESR__0", colnames(noMultLab)):grep("ESR__45", colnames(noMultLab))]
il62 = noMultLab[,grep("INTERLEUKIN.6__0", colnames(noMultLab)):grep("INTERLEUKIN.6__50", colnames(noMultLab))]
il12 = noMultLab[,grep("INTERLEUKIN.1.BETA__0", colnames(noMultLab)):grep("INTERLEUKIN.1.BETA__43", colnames(noMultLab))]
il82 = noMultLab[,grep("INTERLEUKIN.8__0", colnames(noMultLab)):grep("INTERLEUKIN.8__43", colnames(noMultLab))]
tnf2 = noMultLab[,grep("TNF.ALPHA__0", colnames(noMultLab)):grep("TNF.ALPHA__43", colnames(noMultLab))]

# Coagulation
dd2 = noMultLab[,grep("D.DIMER__0", colnames(noMultLab)):grep("D.DIMER__57", colnames(noMultLab))]
b2g2 = noMultLab[,grep("BETA.2.GLYCOPROTEIN.I.ANTIBODY__0", colnames(noMultLab)):grep("BETA.2.GLYCOPROTEIN.I.ANTIBODY__8", colnames(noMultLab))] # many missing
ptt2 = noMultLab[,grep("PTT__0", colnames(noMultLab)):grep("PTT__70", colnames(noMultLab))]
pt2 = noMultLab[,grep("PROTHROMBIN.TIME__0", colnames(noMultLab)):grep("PROTHROMBIN.TIME__70", colnames(noMultLab))]
fib2 = noMultLab[,grep("FIBRINOGEN__0", colnames(noMultLab)):grep("FIBRINOGEN__60", colnames(noMultLab))]
inr2 = noMultLab[,grep("INR__0", colnames(noMultLab)):grep("INR__70", colnames(noMultLab))]
phph2 = noMultLab[,grep("PHOSPHOLIPID.ANTIBODY__0", colnames(noMultLab)):grep("PHOSPHOLIPID.ANTIBODY__15", colnames(noMultLab))] # many missing
ant2 = noMultLab[,grep("ANTICARDIOLIPIN.ANTIBODY__0", colnames(noMultLab)):grep("ANTICARDIOLIPIN.ANTIBODY__8", colnames(noMultLab))] # many missing

# Kidney
eg2 = noMultLab[,grep("EGFR__0", colnames(noMultLab)):grep("EGFR__72", colnames(noMultLab))]
ura2 = noMultLab[,grep("URIC.ACID__0", colnames(noMultLab)):grep("URIC.ACID__39", colnames(noMultLab))]
sod2 = noMultLab[,grep("SODIUM__0", colnames(noMultLab)):grep("SODIUM__67", colnames(noMultLab))]
ch2 = noMultLab[,grep("CHLORIDE__0", colnames(noMultLab)):grep("CHLORIDE__72", colnames(noMultLab))]
cal2 = noMultLab[,grep("CALCIUM__0", colnames(noMultLab)):grep("CALCIUM__72", colnames(noMultLab))]
pot2 = noMultLab[,grep("POTASSIUM__0", colnames(noMultLab)):grep("POTASSIUM__72", colnames(noMultLab))]
bun2 = noMultLab[,grep("BUN__0", colnames(noMultLab)):grep("BUN__72", colnames(noMultLab))]
crt2 = noMultLab[,grep("SERUM.CREATININE__0", colnames(noMultLab)):grep("SERUM.CREATININE__72", colnames(noMultLab))]

# Cardiac
bnp2 = noMultLab[,grep("BRAIN.NATRIURETIC.PROTEIN__0", colnames(noMultLab)):grep("BRAIN.NATRIURETIC.PROTEIN__56", colnames(noMultLab))]
trop2 = noMultLab[,grep("TROPONIN.I__0", colnames(noMultLab)):grep("TROPONIN.I__56", colnames(noMultLab))]
crtKmb2 = noMultLab[,grep("CREATINE_KINASE_MB__0", colnames(noMultLab)):grep("CREATINE_KINASE_MB__56", colnames(noMultLab))]

#arteries/veins
pha2 = noMultLab[,grep("PH.ARTERIAL__0", colnames(noMultLab)):grep("PH.ARTERIAL__56", colnames(noMultLab))]
phv2 = noMultLab[,grep("PH.VENOUS__0", colnames(noMultLab)):grep("PH.VENOUS__57", colnames(noMultLab))]
phc2 = noMultLab[,grep("PH.CAPILLARY__0", colnames(noMultLab)):grep("PH.CAPILLARY__19", colnames(noMultLab))]
ba2 = noMultLab[,grep("BASE.EXCESS.ARTERIAL__0", colnames(noMultLab)):grep("BASE.EXCESS.ARTERIAL__40", colnames(noMultLab))]
bv2 = noMultLab[,grep("BASE.EXCESS.VENOUS__0", colnames(noMultLab)):grep("BASE.EXCESS.VENOUS__57", colnames(noMultLab))]
la2 = noMultLab[,grep("LACTATE_ARTERIAL__0", colnames(noMultLab)):grep("LACTATE_ARTERIAL__56", colnames(noMultLab))]
pca2 = noMultLab[,grep("PCO2.ARTERIAL__0", colnames(noMultLab)):grep("PCO2.ARTERIAL__56", colnames(noMultLab))]
pcv2 = noMultLab[,grep("PC02.VENOUS__0", colnames(noMultLab)):grep("PC02.VENOUS__57", colnames(noMultLab))]
pcc2 = noMultLab[,grep("PCO2.CAPILLARY__0", colnames(noMultLab)):grep("PCO2.CAPILLARY__19", colnames(noMultLab))]
poa2 = noMultLab[,grep("PO2.ARTERIAL__0", colnames(noMultLab)):grep("PO2.ARTERIAL__56", colnames(noMultLab))]
pov2 = noMultLab[,grep("PO2.VENOUS__0", colnames(noMultLab)):grep("PO2.VENOUS__57", colnames(noMultLab))]
oa2 = noMultLab[,grep("O2.SATURATION.ARTERIAL__0", colnames(noMultLab)):grep("O2.SATURATION.ARTERIAL__56", colnames(noMultLab))]
ov2 = noMultLab[,grep("O2.SATURATION.VENOUS__0", colnames(noMultLab)):grep("O2.SATURATION.VENOUS__57", colnames(noMultLab))]
oc2 = noMultLab[,grep("O2.SATURATION.CAPILLARY__0", colnames(noMultLab)):grep("O2.SATURATION.CAPILLARY__9", colnames(noMultLab))]
hca2 = noMultLab[,grep("HCO3.ARTERIAL__0", colnames(noMultLab)):grep("HCO3.ARTERIAL__56", colnames(noMultLab))]
hcv2 = noMultLab[,grep("HCO3.VENOUS__0", colnames(noMultLab)):grep("HCO3.VENOUS__51", colnames(noMultLab))]
hcc2 = noMultLab[,grep("HCO3.CAPILLARY__0", colnames(noMultLab)):grep("HCO3.CAPILLARY__19", colnames(noMultLab))]

## function to create observational units using 7 day window size, options = T whether to include M3
addObs7Day = function(rowsAdd = FALSE, df,col_name,cov,test,dis, death = FALSE, deathM = NA, options = F){

  # the trajectory vector
  trj = c(NA,NA,NA,NA,NA,NA)
  # the number of numbers added to the trajectory on each day
  adds = c(1,1,1,1,1,1)
  # observational units
  obs = data.frame()
  
  # create observations from every length of stay
  for(j in 1:dis){
    arr = num2(as.factor(cov[j])) # get the measurements from that day
    if(j <= length(cov)){
      if(test[j] != "" & !is.na(test[j])){ # add the measurements to trajectories
        arr = num2(as.factor(cov[j]))
        trj = c(trj, arr)
        adds = c(adds,length(arr)) # keeps track of how many measurements were added each day
      } else{
        # append to t, append NA to traj
        trj = c(trj,NA)
        adds = c(adds,1)
      }
    } else { # append NA not measured
      trj = c(trj,NA)
      adds = c(adds,1)
    }
    
    # wind will store the last 7 day measurements
    ind = adds[(length(adds)-6):length(adds)]
    wind = trj[(length(trj)-sum(ind)+1):length(trj)]
    if(options){
      if(sum(is.na(wind)) > 2){ # if you dont have measurements from 5 days or more, model 2 and 3 get NA
        x = data.frame(wind[!is.na(wind)][length(wind[!is.na(wind)])],
                       wind[length(wind)],NA, NA, NA, NA,  NA,  NA, NA, NA,NA,NA,NA,NA)
      } else { # otherwise, fit trajectories, obtain summaries or coef
        #obtain the time vector (due to variable number of measurements per day)
        time=c()
        for(t in 1:7){
          if(ind[t] == 1){
            time = c(time,t)
          } else {
            for(q in 0:(ind[t]-1)){
              time = c(time, t+q*(1/ind[t]))
            }
          }
        }
        # fit a smoothing spline
        coefs = smooth.spline(time[!is.na(wind)],wind[!is.na(wind)],cv = TRUE, nknots = 4)$fit$coef
        x = data.frame(wind[!is.na(wind)][length(wind[!is.na(wind)])],wind[7],wind[1],wind[!is.na(wind)][1],min(wind,na.rm = T),max(wind,na.rm = T),median(wind,na.rm = T), sum(!is.na(wind)),coefs[1],coefs[2],coefs[3],coefs[4],coefs[5],coefs[6]) 
        # add column names
        colnames(x) = paste(col_name, c("last_recorded","last","first","first_recorded","min","max", "median","count","0","1","2","3","4","5"), sep = "_")
      }
    } else {
      #print("im here")
      if(sum(!is.na(wind)) == 0){ # if you dont have measurements from 5 days or more, model 2 and 3 get NA
        #print("w")
        x = data.frame(NA,NA,NA, NA, NA, NA,  NA,  NA)
        
      } else{
        x = data.frame(wind[!is.na(wind)][length(wind[!is.na(wind)])],
                       wind[7],wind[1],wind[!is.na(wind)][1],min(wind,na.rm = T),
                       max(wind,na.rm = T),median(wind,na.rm = T), sum(!is.na(wind)))
      }
      # add column names
      #print("h")
      colnames(x) = paste(col_name, c("last_recorded","last","first","first_recorded","min","max", "median","count"), sep = "_")
    }
    
    
    
    
    if(rowsAdd){ # whether to add rows for each patient(i.e. a patient with LoS of 4 days gets 4 rows), and add a death marker
      if(j!=dis){
        x = cbind(x,death = 0)
      } else if(j==dis){
        x = cbind(x,death = deathM)
      }
      obs = rbind(obs,cbind(df,x))
    } else { # otherwise, rows have already been added
      obs = rbind(obs,x)
    }
  }
  if(rowsAdd){
    return(obs)
  }
  return(cbind(df,obs))
}

# start the cluster to use 6 cores
myCluster <- makeCluster(6, # number of cores to use
                         type = "PSOCK") # type of cluster
registerDoParallel(myCluster)
clusterEvalQ(myCluster, library(splines)) # load the spline library on each core
clusterEvalQ(myCluster, library(stringr)) # load the spline library on each core

## Run this segment of code twice once for arrCovidPos, and arrnoCovidPos
spObs = c() # store the covid positive/negative observations
spObs = foreach(i = arrCovidPos) %dopar% { # change to arrnoCovidPos for covid negative
#for(i in arrCovidPos){
  los = obsStatic$LOS[i]
  if(los > 0){
    row = obsStatic[i,-which(names(obsStatic) %in% c("death"))]
    
    ### Adding Vitals
    row = addObs7Day(rowsAdd = TRUE, row,"temp",sapply(tmp2[i,],num),tmp2[i,],los,death = TRUE, deathM = noMultLab$death[i])
    row = addObs7Day(FALSE, row,"HR",sapply(HR2[i,],num),HR2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"RR",sapply(RR2[i,],num),RR2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"Sys",sapply(Sys2[i,],num),Sys2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"Dia",sapply(Dia2[i,],num),Dia2[i,],los,death = FALSE) 
    row = addObs7Day(FALSE, row,"O2",sapply(O22[i,],num),O22[i,],los,death = FALSE)
    
    ### CBC columns to be added
    row = addObs7Day(FALSE, row,"wbc",sapply(wbc2[i,],num),wbc2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"rbc",sapply(rbc2[i,],num),rbc2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"plt",sapply(plt2[i,],num),plt2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"neu",sapply(neu2[i,],num),neu2[i,],los,death = FALSE) 
    row = addObs7Day(FALSE, row,"eos",sapply(eos2[i,],num),eos2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"mono",sapply(mono2[i,],num),mono2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"bas",sapply(bas2[i,],num),bas2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"lym",sapply(lym2[i,],num),lym2[i,],los,death = FALSE) 
    row = addObs7Day(FALSE, row,"hemo",sapply(hemo2[i,],num),hemo2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"mch",sapply(mch2[i,],num),mch2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"mchc",sapply(mchc2[i,],num),mchc2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"mcv",sapply(mcv2[i,],num),mcv2[i,],los,death = FALSE) 
    row = addObs7Day(FALSE, row,"hemat",sapply(hemat2[i,],num),hemat2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"anion",sapply(an2[i,],num),an2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"ProCalcitonin",sapply(procal2[i,],num),procal2[i,],los,death = FALSE)
    
    # kidney columns to be added
    row = addObs7Day(FALSE, row,"sod",sapply(sod2[i,],num),sod2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"chl",sapply(ch2[i,],num),ch2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"cal",sapply(cal2[i,],num),cal2[i,],los,death = FALSE) 
    row = addObs7Day(FALSE, row,"pot",sapply(pot2[i,],num),pot2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"bun",sapply(bun2[i,],num),bun2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"crt",sapply(crt2[i,],num),crt2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"eGFR",sapply(eg2[i,],num),eg2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"uric",sapply(ura2[i,],num),ura2[i,],los,death = FALSE)
    
    # liver columns to be added
    row = addObs7Day(FALSE, row,"glu",sapply(glu2[i,],num),glu2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"alb",sapply(alb2[i,],num),alb2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"alt",sapply(alt2[i,],num),alt2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"ast",sapply(ast2[i,],num),ast2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"alp",sapply(alp2[i,],num),alp2[i,],los,death = FALSE) 
    row = addObs7Day(FALSE, row,"bili",sapply(bili2[i,],num),bili2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"crtK",sapply(crtK2[i,],num),crtK2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"LDH",sapply(ldh2[i,],num),ldh2[i,],los,death = FALSE) 
    
    # Arteries / Venous columns to be added
    row = addObs7Day(FALSE, row,"phArt",sapply(pha2[i,],num),pha2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"phVen",sapply(phv2[i,],num),phv2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"phCap",sapply(phc2[i,],num),phc2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"BaseArt",sapply(ba2[i,],num),ba2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"BaseVen",sapply(bv2[i,],num),bv2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"LactArt",sapply(la2[i,],num),la2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"pCO2Art",sapply(pca2[i,],num),pca2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"pCO2Ven",sapply(pcv2[i,],num),pcv2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"pCO2Cap",sapply(pcc2[i,],num),pcc2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"pO2Art",sapply(poa2[i,],num),poa2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"pO2Ven",sapply(pov2[i,],num),pov2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"O2Art",sapply(oa2[i,],num),oa2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"O2Ven",sapply(ov2[i,],num),ov2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"O2Cap",sapply(oc2[i,],num),oc2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"HCO3Art",sapply(hca2[i,],num),hca2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"HCO3Ven",sapply(hcv2[i,],num),hcv2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"HCO3Cap",sapply(hcc2[i,],num),hcc2[i,],los,death = FALSE)
    
    # Inflammatory columns to be added
    row = addObs7Day(FALSE, row,"crp",sapply(crp2[i,],num),crp2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"Ferritin",sapply(fer2[i,],num),fer2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"ESR",sapply(esr2[i,],num),esr2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"IL6",sapply(il62[i,],num),il62[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"IL8",sapply(il82[i,],num),il82[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"IL1B",sapply(il12[i,],num),il12[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"tnf",sapply(tnf2[i,],num),tnf2[i,],los,death = FALSE)
    
    # Coagulation State columns to be added
    row = addObs7Day(FALSE, row,"ddimer",sapply(dd2[i,],num),dd2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"Beta2Glycoprotein",sapply(b2g2[i,],num),b2g2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"PTT",sapply(ptt2[i,],num),ptt2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"PT",sapply(pt2[i,],num),pt2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"Fibrinogen",sapply(fib2[i,],num),fib2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"INR",sapply(inr2[i,],num),inr2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"Phospholipid",sapply(phph2[i,],num),phph2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"Anticardiolipin",sapply(ant2[i,],num),ant2[i,],los,death = FALSE)
    
    # Cardiac columns to be added
    row = addObs7Day(FALSE, row,"BNP",sapply(bnp2[i,],num),bnp2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"trop",sapply(trop2[i,],num),trop2[i,],los,death = FALSE)
    row = addObs7Day(FALSE, row,"crtKmb",sapply(crtKmb2[i,],num),crtKmb2[i,],los,death = FALSE)
    
  }
}

stopCluster(myCluster)
spObs = rbindlist(spObs) #parallel case
write.csv(spObs,file="covidPositiveM2.csv")
# store it and read it when necessary 
# write.csv(spObs,file="covidNegatice.csv") # and covidNegative.csv
spObs = read.csv("covidPositive.csv")
spObs = spObs[,2:dim(spObs)[2]]
## add Length of Stay for each observational units
losVector = c()
ids = c()
for (i in arrCovidPos){ # or arrnoCovidPos
  los =  obsStatic$LOS[i]
  if(los > 0){
    for(j in 1:los){
      ids =  c(noMultLab$NEW_MASKED_MRN[i],ids)
      losVector = c(losVector,j)    
    }
  }
}
spObs$LOS = losVector


for (i in arrCovidPos){ # or arrnoCovidPos
  if(los > 0){
    for(j in 1:los){
      ids =  c(noMultLab$NEW_MASKED_MRN[i],ids)
    }
  }
}
spObs$LOS = losVector
spObs$ids = factor(ids)
spObs$death = factor(spObs$death)

# declutter get rid of columns with 5 or less non na values (gets rid of over 100 columns)
cols = as.numeric(which(colSums(!is.na(spObs)) > 5))
spObs = data.frame(spObs)[,cols]


# numerical columns to find the boundaries 
num_cols <- unlist(lapply(spObs, is.numeric))
spObs = data.frame(spObs)
min(sapply(spObs[,num_cols],min, na.rm=T),na.rm=T) # -13342.05
max(sapply(spObs[,num_cols],max, na.rm=T),na.rm = T) # 234037.5
# Replace NA columns with -10^10
covidPos7 = do.call(data.frame,lapply(spObs, function(x) replace(x, is.na(x),-10^10)))
covidPos7 = covidPos7[,-c(grep("cc", colnames(covidPos7)))] # no need for covid pos/negative column, as all the data is 1 type

# make sure categorical variables are categories
covidPos7$ICU = factor(covidPos7$ICU)
covidPos7$INPATIENT_NON_ICU = factor(covidPos7$INPATIENT_NON_ICU)
covidPos7$ASTHMA = factor(covidPos7$ASTHMA)
covidPos7$COPD = factor(covidPos7$COPD)
covidPos7$HTN = factor(covidPos7$HTN)
covidPos7$OBSTRUCTIVE_SLEEP_APNEA = factor(covidPos7$OBSTRUCTIVE_SLEEP_APNEA)
covidPos7$OBESITY = factor(covidPos7$OBESITY)
covidPos7$DIABETES = factor(covidPos7$DIABETES)
covidPos7$CHRONIC_KIDNEY_DISEASE = factor(covidPos7$CHRONIC_KIDNEY_DISEASE)
covidPos7$HIV_FLAG = factor(covidPos7$HIV_FLAG)
covidPos7$CANCER_FLAG = factor(covidPos7$CANCER_FLAG)
covidPos7$CORONARY_ARTERY_DISEASE = factor(covidPos7$CORONARY_ARTERY_DISEASE)
covidPos7$ATRIAL_FIBRILLATION = factor(covidPos7$ATRIAL_FIBRILLATION)
covidPos7$HEART_FAILURE = factor(covidPos7$HEART_FAILURE)
covidPos7$CHRONIC_VIRAL_HEPATITIS = factor(covidPos7$CHRONIC_VIRAL_HEPATITIS)
covidPos7$ALCOHOLIC_NONALCOHOLIC_LIVER_DISEASE = factor(covidPos7$ALCOHOLIC_NONALCOHOLIC_LIVER_DISEASE)
covidPos7$ARDS = factor(covidPos7$ARDS)
covidPos7$ACUTE_KIDNEY_INJURY = factor(covidPos7$ACUTE_KIDNEY_INJURY)
covidPos7$ACUTE_VENOUS_THROMBOEMBOLISM = factor(covidPos7$ACUTE_VENOUS_THROMBOEMBOLISM)
covidPos7$CEREBRAL_INFARCTION = factor(covidPos7$CEREBRAL_INFARCTION)
covidPos7$INTRACEREBRAL_HEMORRHAGE = factor(covidPos7$INTRACEREBRAL_HEMORRHAGE)
covidPos7$ACUTE_MI = factor(covidPos7$ACUTE_MI)

# separate Model 1, Model 2 and Model 3 observational units, use only patients who are staying 5 or more days
observationM1 = covidPos7[covidPos7$LOS >= 5, c(1:31,grep("ids", colnames(covidPos7)), grep("\\last_recorded$", colnames(covidPos7)),grep("death", colnames(covidPos7)))]
observationM2 = covidPos7[covidPos7$LOS >= 5, -c(grep("\\_1$", colnames(covidPos7)),grep("\\_2$", colnames(covidPos7)),grep("\\_3$", colnames(covidPos7)),grep("\\_0$", colnames(covidPos7)),grep("\\_4$", colnames(covidPos7)),grep("\\_5$", colnames(covidPos7)),grep("_first$", colnames(covidPos7)),grep("\\last$", colnames(covidPos7)))]
observationM3 = covidPos7[covidPos7$LOS >= 5, c(1:31,grep("\\last_recorded$", colnames(covidPos7)),grep("ids", colnames(covidPos7)),grep("\\_1$", colnames(covidPos7)),grep("\\_2$", colnames(covidPos7)),grep("\\_3$", colnames(covidPos7)),grep("\\_0$", colnames(covidPos7)),grep("\\_4$", colnames(covidPos7)),grep("\\_5$", colnames(covidPos7)),grep("death", colnames(covidPos7)))]

# separate Model 1 vs Model 2 only (all Length of stays)
observationM1 = covidPos7[, c(1:31, grep("\\last_recorded$", colnames(covidPos7)),grep("death", colnames(covidPos7)))]
observationM2 = covidPos7[, -c(grep("\\_1$", colnames(covidPos7)),grep("\\_2$", colnames(covidPos7)),grep("\\_3$", colnames(covidPos7)),grep("\\_0$", colnames(covidPos7)),grep("\\_4$", colnames(covidPos7)),grep("\\_5$", colnames(covidPos7)),grep("_first$", colnames(covidPos7)),grep("\\last$", colnames(covidPos7)))]

### A classification using random forest, we use Ranger Library for its faster performance compared to Random Forest. 
# only has last recorded (most recent lab)
fineTune(c(10:20),observationM1) # pick the best mtry
observationM1RF <- ranger(formula = death ~ ., data = observationM1, importance = "impurity",verbose = T, probability = T)
print(observationM1RF)
# has trajectory summaries (min/max/last/first/median)
fineTune(seq(20,75,by=2),observationM2) # pick the best mtry
observationM2RF <- ranger(formula = death ~ ., data = observationM2, importance = "impurity",verbose = T, probability = T, mtry = 26)
print(observationM2RF)
# has trajectory coefficients (Splines)
fineTune(seq(20,75,by=2),observationM3) # pick the best mtry
observationM3RF <- ranger(formula = death ~ ., data = observationM3, importance = "impurity",verbose = T, probability = T, mtry = 24)
print(observationM3RF)


# Model Evaluation --------------------------------------------------------
# Can be Repeated for Covid Negative
# Variable Importance using Gini Score
sort(importance(observationM1RF),decreasing = T)[1:30]
sort(importance(observationM2RF),decreasing = T)[1:30]
sort(importance(observationM3RF),decreasing = T)[1:30]

# Importance Plot for Model 3
par(mai=c(1,2.4,1,1))
vals = (sort(importance(observationM3RF),decreasing = T)[1:30])[order(sort(importance(observationM3RF),decreasing = T)[1:30], decreasing = F)]
names = c("O2 Saturation last","Systolic BP last","Diastolic BP last","Potassium last","O2 5th", "Heart Rate last","Potassium 5th","Anion Gap Last","Lactate Arterial last","Temperature last","BUN last","Respiratory Rate last","Heart Rate 5th","Systolic BP 5th","Chloride last","AST last","Temperature 5th","Creatinine last","O2 4th","Glucose last","O2 3rd","White Blood Cell last","Age", "Systolic BP 3rd", "Heart Rate 3rd", "Heart Rate 4th","Temperature 4th","Diastolic BP 5th", "Anion Gap 5th","Troponin last")
barplot(vals,col = "lavender", border = F, horiz = T,cex.names = 1.2,las=1 , xlim =c(0,90),names.arg = rev(names))

# Importance Plot for Model 2
par(mai=c(1,2.6,1,1))
vals = (sort(importance(observationM2RF),decreasing = T)[1:30])[order(sort(importance(observationM2RF),decreasing = T)[1:30], decreasing = F)]
# covid Pos names
names = c("O2 Saturation last","Systolic BP last","Diastolic BP last","Potassium last","Heart Rate last","Temperature last", "Anion Gap Last","Lactate Arterial last","BUN last","Respiratory Rate last", "AST last","Chloride last","Creatinine last","White Blood Cell last","Glucose last","O2 min","Age","Heart Rate max","Heart Rate Median", "Systolic BP min","PCO2 Aterial last","Calcium last","Systolic median","Neutrophil last","Troponin last","MCV last","Length of Stay","Hemoglobin last","Temperature max","Systolic BP max") 
# covid Neg names
names = c("O2 Saturation last","Systolic BP last","Diastolic BP last","Temperature last", "Anion Gap Last","Lactate Arterial last", "AST last","BUN last","Heart Rate last","Glucose last","White Blood Cell last","Length of Stay","Respiratory Rate last","Age","Chloride last","Neutrophil last","Temperature first","Systolic BP first","Temperature max","Potassium last","Calcium last","Platelet last","ALT last","Systolic max","Glucose Min","LDH last","MCHC last","Systolic BP Median","MCV last","Creatinine last") 
barplot(vals,col = "lavender", border = F, horiz = T,cex.names = 1.2,las=1 , xlim =c(0,75))#,names.arg = rev(names))


##### ROC Plot
# plots the ROC curves for our 3 models
obj3 = roc(observationM3$death~observationM3RF$predictions[,2], data = observationM3, plot = TRUE, main = "ROC CURVE", col= "blue", ci=TRUE)
obj2 = roc(observationM2$death~observationM2RF$predictions[,2], data = observationM2, plot = TRUE, main = "ROC CURVE", col= "blue", ci=TRUE)
obj1 = roc(observationM1$death~observationM1RF$predictions[,2], data = observationM1, plot = TRUE, main = "ROC CURVE", col= "blue", ci=TRUE)
# use bootstrap to sample 2000 times and create 2000 ROC curves
ciobj3 = ci.se(obj3, specificities=seq(0, 1, l=25))
ciobj2 = ci.se(obj2, specificities=seq(0, 1, l=25))
ciobj1 = ci.se(obj1, specificities=seq(0, 1, l=25))

# reformat the result
dat.ci3 = data.frame(x = as.numeric(rownames(ciobj3)),
                     lower = ciobj3[, 1],
                     med = ciobj3[, 2],
                     upper = ciobj3[, 3])
dat.ci2 = data.frame(x = as.numeric(rownames(ciobj2)),
                     lower = ciobj2[, 1],
                     med = ciobj2[, 2],
                     upper = ciobj2[, 3])
dat.ci1 = data.frame(x = as.numeric(rownames(ciobj1)),
                     lower = ciobj1[, 1],
                     med = ciobj1[, 2],
                     upper = ciobj1[, 3])

# plot the CI ROCs 
ggplot(data=dat.ci2, aes(x=(1-x), y=med)) +
  geom_line(color = "limegreen")+
  theme_minimal() + 
  geom_abline(slope=1, intercept = 0, linetype = "dashed", 
              alpha=0.7, color = "red") + # the line for a random predictor
  coord_equal() + # square plot
  geom_ribbon(data = dat.ci2, aes(x = (1-x), ymin = lower, ymax = upper), 
              fill = "limegreen", alpha= 0.2) + # the ribbon representing CI
  #ggtitle("95% Confidence Interval ROC curve") + 
  geom_ribbon(data = dat.ci3, aes(x = (1-x), ymin = lower, ymax = upper), 
              fill = "dodgerblue", alpha= 0.2)+
  geom_ribbon(data = dat.ci1, aes(x = (1-x), ymin = lower, ymax = upper), 
              fill = "firebrick1", alpha= 0.2)+
  geom_line(data = dat.ci3, aes(x = (1-x), y = med), color = "dodgerblue")+
  geom_line(data = dat.ci1, aes(x = (1-x), y = med), color = "firebrick1")+
  theme(text = element_text(size=15), legend.position = "bottomright")+
  ylab("Sensitivity")+
  xlab("Specificity")# font size


#### Precision Recall Plot
pr3 = pr.curve(scores.class0 = observationM3RF$predictions[observationM3$death == 1,2], scores.class1 = observationM3RF$predictions[observationM3$death == 0,2],curve = T)
pr2 = pr.curve(scores.class0 = observationM2RF$predictions[observationM2$death == 1,2], scores.class1 = observationM2RF$predictions[observationM2$death == 0,2],curve = T)
pr1 = pr.curve(scores.class0 = observationM1RF$predictions[observationM1$death == 1,2], scores.class1 = observationM1RF$predictions[observationM1$death == 0,2],curve = T)

plot(smooth.spline(pr2$curve[,1],pr2$curve[,2]),type = "l", xlab = "Recall", ylab="Average Precision", lwd = 2,col = "limegreen", cex.lab = 1.25,las = 1, axes = T)
lines(smooth.spline(pr1$curve[,1],pr1$curve[,2]),type = "l",col="firebrick1", lwd = 2)
lines(smooth.spline(pr3$curve[,1],pr3$curve[,2]),type = "l",col="dodgerblue1", lwd = 2)
segments(0.5,0,0.5,0.3,lty=2)
segments(-0.1,0.296,0.5,0.296,lty=2)
segments(-0.1,0.296,0.5,0.296,lty=3)
segments(-0.1,0.205,0.5,0.205,lty=4)
legend("topright", legend = c("M1","M2"), col = c("firebrick1","limegreen"),lwd = 2, border = F)
ticks = c(0,0.20,0.29,0.3,0.35)
axis(side = 2, at = ticks, las = 1, cex.lab = 0.8)
axis(side = 1)


### Monte Carlo Cross Validation 20 roc,pr values for Wilcoxon Paired test
roc = data.frame(NA, nrow = 20, ncol = 3)
pr = data.frame(NA, nrow = 20, ncol = 3)
for(i in (1:20)){
  # obtain training/test subsets
  # train.idx <- sample(x = dim(covidPos7[covidPos7$LOS >= 5, ])[1], size = ceiling(dim(covidPos7[covidPos7$LOS >= 5, ])[1]*0.7))
  # train.idx <- sample(x = dim(covidPos7)[1], size = ceiling(dim(covidPos7)[1]*0.7)) (observationM1)[train.idx,] 
  
  # make test set different set of patients
  trainLevels = sample(levels(observationM1$ids), size = ceiling(length(levels(observationM1$ids))*0.7))
  train1 = (observationM1)[observationM1$ids %in% trainLevels,] 
  test1 = (observationM1)[!(observationM1$ids %in% trainLevels),] 
  train2 = (observationM2)[observationM2$ids %in% trainLevels,] 
  test2 = (observationM2)[!(observationM1$ids %in% trainLevels),] 
  train3 = observationM3[observationM3$ids %in% trainLevels,] 
  test3 = observationM3[!(observationM3$ids %in% trainLevels),] 
 
  
  # fit ranger random forests 
  rf1 <- ranger(death ~ ., data = train1, write.forest = TRUE,probability = TRUE)
  rf2 <- ranger(death ~ ., data = train2, write.forest = TRUE,probability = TRUE, mtry = 50)
  rf3 <- ranger(death ~ ., data = train3, write.forest = TRUE,probability = TRUE, mtry = 50)
  pred1 <- predict(rf1, data = test1)$predictions
  pred2 <- predict(rf2, data = test2)$predictions
  pred3 <- predict(rf3, data = test3)$predictions
  
  # Obtain ROC, PR AUC values
  roc[i,1] <- roc.curve(scores.class0 = pred1[test1$death == 1,2], scores.class1 = pred1[test1$death == 0,2])$auc
  pr[i,1] <- pr.curve(scores.class0 = pred1[test1$death == 1,2], scores.class1 = pred1[test1$death == 0,2])$auc.davis.goadrich
  roc[i,2] <- roc.curve(scores.class0 = pred2[test2$death == 1,2], scores.class1 = pred2[test2$death == 0,2])$auc
  pr[i,2] <- pr.curve(scores.class0 = pred2[test2$death == 1,2], scores.class1 = pred2[test2$death == 0,2])$auc.davis.goadrich
  roc[i,3] <- roc.curve(scores.class0 = pred3[test3$death == 1,2], scores.class1 = pred3[test3$death == 0,2])$auc
  pr[i,3] <- pr.curve(scores.class0 = pred3[test3$death == 1,2], scores.class1 = pred3[test3$death == 0,2])$auc.davis.goadrich
  print(i)
}

rocpr = cbind(roc,pr)
write.csv(rocpr,"rocprCovidPosM2.csv") 

# model 1 vs. model 2 ROC AUC
wilcox.test(rocpr[,1],rocpr[,2],paired = T)
# model 2 vs. model 3 ROC AUC
wilcox.test(rocpr[,2],rocpr[,3],paired = T)
# model 1 vs. model 3 ROC AUC
wilcox.test(rocpr[,1],rocpr[,3],paired = T)
# model 1 vs. model 2 PR AUC
wilcox.test(rocpr[,4],rocpr[,5],paired = T)
# model 2 vs. model 3 PR AUC
wilcox.test(rocpr[,5],rocpr[,6],paired = T)
# model 1 vs. model 3 PR AUC
wilcox.test(rocpr[,4],rocpr[,6],paired = T)


### Interpretability Using SHEP save M2 and M3 observations for Python
write.csv(observationM1,"M1Neg.csv") 
write.csv(observationM2,"M2Neg.csv") 
write.csv(observationM3,"M3Neg.csv") 


# Let's fine tune the fit
data(observationM1)
fineTune = function(mtry, observation){
  grid <-  expand.grid(mtry = mtry,splitrule = "gini", min.node.size=10)
  
  fitControl <- trainControl(method = "CV",
                             number = 5,
                             verboseIter = TRUE)
  
  observation$death <- factor(observation$death, levels=c(0,1), labels=c("Survivor", "Non-Survivor"))
  fit = train(
    x = observation[ , names(observation) != 'death'],
    y = observation[ , names(observation) == 'death'],
    method = 'ranger',
    num.trees = 500,
    tuneGrid = grid,
    trControl = fitControl
  )
  fit
}

## Supplementary Plots; hazard rate Covid Neg and Positive
oP = read.csv("M1Pos.csv")
oN = read.csv("M1Neg.csv")

# Plot histogram of length of stay in the hospital based on whether covid +/- and survived yes or no
oP$Death = factor(oP$death, labels = c("Survivor", "non-Survivor"))
ggplot(oP, aes(x=LOS, fill=Death)) +
  geom_histogram(binwidth = 1,aes(fill = Death )) +
  scale_fill_manual(values=c("palegreen1","red")) +
  theme_minimal() +
  labs(fill = "") + 
  theme(text = element_text(size=15))+
  xlab("Length of Stay") +
  ylab("Frequency")  +
  scale_x_continuous(breaks = c(seq(5,20,5),seq(29,59,10)),labels = c(seq(5,20,5),seq(30,60,10)))

## Some Discussion Variables
wilcox.test(noMultLab$TEMPERATURE[noMultLab$cc == "Detected" & noMultLab$death == 1],noMultLab$TEMPERATURE[noMultLab$cc == "Detected"& noMultLab$death == 0], alternative = "less")

# Ablation
observationM2 = covidPos7[, -c(grep("\\_1$", colnames(covidPos7)),grep("\\_2$", colnames(covidPos7)),grep("\\_3$", colnames(covidPos7)),grep("\\_0$", colnames(covidPos7)),grep("\\_4$", colnames(covidPos7)),grep("\\_5$", colnames(covidPos7)),grep("_first$", colnames(covidPos7)),grep("\\last$", colnames(covidPos7)))]
Ob2f = (observationM2)[,-c(grep("_first_recorded", colnames(observationM2)))] 
Ob2c = (observationM2)[,-c(grep("\\count$", colnames(observationM2)))] 
Ob2mm = (observationM2)[,-c(grep("\\min$", colnames(observationM2)),grep("\\max$", colnames(observationM2)))] 
Ob2m = (observationM2)[,-c(grep("\\median$", colnames(observationM2)))] 
Ob2l = (observationM2)[,-c(grep("_last_recorded$", colnames(observationM2)))] 
rf1 <- ranger(death ~ ., data = observationM2, write.forest = TRUE,probability = TRUE)
rf2 <- ranger(death ~ ., data = Ob2f, write.forest = TRUE,probability = TRUE, mtry = 50)
rf3 <- ranger(death ~ ., data = Ob2c, write.forest = TRUE,probability = TRUE, mtry = 50)
rf4 <- ranger(death ~ ., data = Ob2mm, write.forest = TRUE,probability = TRUE, mtry = 50)
rf5 <- ranger(death ~ ., data = Ob2m, write.forest = TRUE,probability = TRUE, mtry = 50)
rf6 <- ranger(death ~ ., data = Ob2l, write.forest = TRUE,probability = TRUE, mtry = 50)

obj1 = roc(observationM2$death~rf1$predictions[,2], data = observationM2, plot = TRUE, main = "ROC CURVE", col= "blue", ci=TRUE)
obj2 = roc(Ob2f$death~rf2$predictions[,2], data = Ob2f, plot = TRUE, main = "ROC CURVE", col= "blue", ci=TRUE)
obj3 = roc(Ob2c$death~rf3$predictions[,2], data = Ob2c, plot = TRUE, main = "ROC CURVE", col= "blue", ci=TRUE)
obj4 = roc(Ob2mm$death~rf4$predictions[,2], data = Ob2mm, plot = TRUE, main = "ROC CURVE", col= "blue", ci=TRUE)
obj5 = roc(Ob2m$death~rf5$predictions[,2], data = Ob2m, plot = TRUE, main = "ROC CURVE", col= "blue", ci=TRUE)
obj6 = roc(Ob2l$death~rf6$predictions[,2], data = Ob2l, plot = TRUE, main = "ROC CURVE", col= "blue", ci=TRUE)

# use bootstrap to sample 2000 times and create 2000 ROC curves
ciobj1 = ci.se(obj1, specificities=seq(0, 1, l=25))
ciobj2 = ci.se(obj2, specificities=seq(0, 1, l=25))
ciobj3 = ci.se(obj3, specificities=seq(0, 1, l=25))
ciobj4 = ci.se(obj4, specificities=seq(0, 1, l=25))
ciobj5 = ci.se(obj5, specificities=seq(0, 1, l=25))
ciobj6 = ci.se(obj6, specificities=seq(0, 1, l=25))

# reformat the result
dat.ci1 = data.frame(x = as.numeric(rownames(ciobj1)),
                     lower = ciobj1[, 1],
                     med = ciobj1[, 2],
                     upper = ciobj1[, 3])
dat.ci2 = data.frame(x = as.numeric(rownames(ciobj2)),
                     lower = ciobj2[, 1],
                     med = ciobj2[, 2],
                     upper = ciobj2[, 3])
dat.ci3 = data.frame(x = as.numeric(rownames(ciobj3)),
                     lower = ciobj3[, 1],
                     med = ciobj3[, 2],
                     upper = ciobj3[, 3])
dat.ci4 = data.frame(x = as.numeric(rownames(ciobj4)),
                     lower = ciobj4[, 1],
                     med = ciobj4[, 2],
                     upper = ciobj4[, 3])
dat.ci5 = data.frame(x = as.numeric(rownames(ciobj5)),
                     lower = ciobj5[, 1],
                     med = ciobj5[, 2],
                     upper = ciobj5[, 3])
dat.ci6 = data.frame(x = as.numeric(rownames(ciobj6)),
                     lower = ciobj6[, 1],
                     med = ciobj6[, 2],
                     upper = ciobj6[, 3])

# plot the CI ROCs 
ggplot(data=dat.ci1, aes(x=(1-x), y=med)) +
  geom_line(color = "limegreen") +
  theme_minimal() + 
  geom_abline(slope=1, intercept = 0, linetype = "dashed", 
              alpha=0.7, color = "red") + # the line for a random predictor
  coord_equal() + # square plot
  #geom_ribbon(data = dat.ci1, aes(x = (1-x), ymin = lower, ymax = upper), 
              #fill = "limegreen", alpha= 0.2) + # the ribbon representing CI
  #ggtitle("95% Confidence Interval ROC curve") + 
  #geom_ribbon(data = dat.ci3, aes(x = (1-x), ymin = lower, ymax = upper), 
             # fill = "dodgerblue", alpha= 0.2)+
  #geom_ribbon(data = dat.ci2, aes(x = (1-x), ymin = lower, ymax = upper), 
             # fill = "firebrick1", alpha= 0.2)+  
  #geom_ribbon(data = dat.ci4, aes(x = (1-x), ymin = lower, ymax = upper), 
             # fill = "lightpink", alpha= 0.2)+
  #geom_ribbon(data = dat.ci5, aes(x = (1-x), ymin = lower, ymax = upper), 
             # fill = "violet", alpha= 0.2)+
  #geom_ribbon(data = dat.ci6, aes(x = (1-x), ymin = lower, ymax = upper), 
             # fill = "orange", alpha= 0.2)+
  geom_line(data = dat.ci3, aes(x = (1-x), y = med), color = "dodgerblue")+
  geom_line(data = dat.ci4, aes(x = (1-x), y = med), color = "lightpink")+
  geom_line(data = dat.ci5, aes(x = (1-x), y = med), color = "violet")+
  geom_line(data = dat.ci6, aes(x = (1-x), y = med), color = "orange")+
  geom_line(data = dat.ci2, aes(x = (1-x), y = med), color = "firebrick1")+
  theme(text = element_text(size=15), legend.position = "bottomright")+
  ylab("Sensitivity")+
  xlab("Specificity")# font size

rp = data.frame(NA, nrow = 20, ncol = 12)
for(i in (1:20)){
  # obtain training/test subsets
  #train.idx <- sample(x = dim(covidPos7[covidPos7$LOS >= 5, ])[1], size = ceiling(dim(covidPos7[covidPos7$LOS >= 5, ])[1]*0.7))
  train.idx <- sample(x = dim(covidPos7)[1], size = ceiling(dim(covidPos7)[1]*0.7))
  train1All = (observationM2)[train.idx,] 
  test1All = (observationM2)[-train.idx,] 
  train2f = (observationM2)[train.idx,-c(grep("_first_recorded", colnames(observationM2)))] 
  test2f = (observationM2)[-train.idx,-c(grep("_first_recorded", colnames(observationM2)))] 
  train2c = (observationM2)[train.idx,-c(grep("\\count$", colnames(observationM2)))] 
  test2c = (observationM2)[-train.idx,-c(grep("\\count$", colnames(observationM2)))] 
  train2mm = (observationM2)[train.idx,-c(grep("\\min$", colnames(observationM2)),grep("\\max$", colnames(observationM2)))] 
  test2mm = (observationM2)[-train.idx,-c(grep("\\min$", colnames(observationM2)),grep("\\max$", colnames(observationM2)))] 
  train2m = (observationM2)[train.idx,-c(grep("\\median$", colnames(observationM2)))] 
  test2m = (observationM2)[-train.idx,-c(grep("\\median$", colnames(observationM2)))]
  train2l = (observationM2)[train.idx,-c(grep("_last_recorded$", colnames(observationM2)))] 
  test2l = (observationM2)[-train.idx,-c(grep("_last_recorded$", colnames(observationM2)))]
  
  # fit ranger random forests 
  rf1 <- ranger(death ~ ., data = train1All, write.forest = TRUE,probability = TRUE)
  rf2 <- ranger(death ~ ., data = train2f, write.forest = TRUE,probability = TRUE, mtry = 50)
  rf3 <- ranger(death ~ ., data = train2c, write.forest = TRUE,probability = TRUE, mtry = 50)
  rf4 <- ranger(death ~ ., data = train2mm, write.forest = TRUE,probability = TRUE, mtry = 50)
  rf5 <- ranger(death ~ ., data = train2m, write.forest = TRUE,probability = TRUE, mtry = 50)
  rf6 <- ranger(death ~ ., data = train2l, write.forest = TRUE,probability = TRUE, mtry = 50)
  pred1 <- predict(rf1, data = test1All)$predictions
  pred2 <- predict(rf2, data = test2f)$predictions
  pred3 <- predict(rf3, data = test2c)$predictions
  pred4 <- predict(rf4, data = test2mm)$predictions
  pred5 <- predict(rf5, data = test2m)$predictions
  pred6 <- predict(rf6, data = test2l)$predictions
  
  # Obtain ROC, PR AUC values
  rp[i,1] <- roc.curve(scores.class0 = pred1[test1All$death == 1,2], scores.class1 = pred1[test1All$death == 0,2])$auc
  rp[i,2] <- pr.curve(scores.class0 = pred1[test1All$death == 1,2], scores.class1 = pred1[test1All$death == 0,2])$auc.davis.goadrich
  rp[i,3] <- roc.curve(scores.class0 = pred2[test2f$death == 1,2], scores.class1 = pred2[test2f$death == 0,2])$auc
  rp[i,4] <- pr.curve(scores.class0 = pred2[test2f$death == 1,2], scores.class1 = pred2[test2f$death == 0,2])$auc.davis.goadrich
  rp[i,5] <- roc.curve(scores.class0 = pred3[test2c$death == 1,2], scores.class1 = pred3[test2c$death == 0,2])$auc
  rp[i,6] <- pr.curve(scores.class0 = pred3[test2c$death == 1,2], scores.class1 = pred3[test2c$death == 0,2])$auc.davis.goadrich
  rp[i,7] <- roc.curve(scores.class0 = pred4[test2mm$death == 1,2], scores.class1 = pred4[test2mm$death == 0,2])$auc
  rp[i,8] <- pr.curve(scores.class0 = pred4[test2mm$death == 1,2], scores.class1 = pred4[test2mm$death == 0,2])$auc.davis.goadrich
  rp[i,9] <- roc.curve(scores.class0 = pred5[test2m$death == 1,2], scores.class1 = pred5[test2m$death == 0,2])$auc
  rp[i,10] <- pr.curve(scores.class0 = pred5[test2m$death == 1,2], scores.class1 = pred5[test2m$death == 0,2])$auc.davis.goadrich
  rp[i,11] <- roc.curve(scores.class0 = pred6[test2l$death == 1,2], scores.class1 = pred5[test2l$death == 0,2])$auc
  rp[i,12] <- pr.curve(scores.class0 = pred6[test2l$death == 1,2], scores.class1 = pred5[test2l$death == 0,2])$auc.davis.goadrich
  print(i)
}

#### Precision Recall Plot
pr1 = pr.curve(scores.class0 = rf1$predictions[observationM2$death == 1,2], scores.class1 = rf1$predictions[observationM2$death == 0,2],curve = T)
pr2 = pr.curve(scores.class0 = rf2$predictions[Ob2f$death == 1,2], scores.class1 = rf2$predictions[Ob2f$death == 0,2],curve = T)
pr3 = pr.curve(scores.class0 = rf3$predictions[Ob2c$death == 1,2], scores.class1 = rf3$predictions[Ob2c$death == 0,2],curve = T)
pr4 = pr.curve(scores.class0 = rf4$predictions[Ob2mm$death == 1,2], scores.class1 = rf4$predictions[Ob2mm$death == 0,2],curve = T)
pr5 = pr.curve(scores.class0 = rf5$predictions[Ob2m$death == 1,2], scores.class1 = rf5$predictions[Ob2m$death == 0,2],curve = T)
pr6 = pr.curve(scores.class0 = rf6$predictions[Ob2l$death == 1,2], scores.class1 = rf6$predictions[Ob2l$death == 0,2],curve = T)

plot(smooth.spline(pr2$curve[,1],pr2$curve[,2]),type = "l", xlab = "Recall", ylab="Precision", lwd = 2,col = "firebrick1", cex.lab = 1.25,las = 1, axes = T, ylim = c(0,0.6))
lines(smooth.spline(pr1$curve[,1],pr1$curve[,2]),type = "l",col="limegreen", lwd = 2)
lines(smooth.spline(pr3$curve[,1],pr3$curve[,2]),type = "l",col="dodgerblue1", lwd = 2)
lines(smooth.spline(pr4$curve[,1],pr4$curve[,2]),type = "l",col="lightpink", lwd = 2)
lines(smooth.spline(pr5$curve[,1],pr5$curve[,2]),type = "l",col="violet", lwd = 2)
lines(smooth.spline(pr6$curve[,1],pr6$curve[,2]),type = "l",col="orange", lwd = 2)
legend("topright", legend = c("All","first","counts","max/min","median","last"), col = c("limegreen","firebrick1","dodgerblue1","lightpink","violet","orange"),lwd = 2, border = F)
ticks = c(0,0.1, 0.2, 0.3, 0.4,0.5,0.6)
axis(side = 2, at = ticks, las = 1, cex.lab = 0.8)
axis(side = 1)

