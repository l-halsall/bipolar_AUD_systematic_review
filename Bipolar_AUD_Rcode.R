## Note – the code is identical for depressive and manic symptoms, so the code is provided for depressive symptomatology only. An identical code was then ran for manic symptomatology. 
## All code begins after the dataset has been uploaded. 

## 1. Depressive symptomatology
install.packages("metafor")
install.packages("meta")
install.packages("moments")
library(metafor)
library(meta)
library(moments)
SMD <- escalc(measure = "SMD", m1i = data_depression$int_mean, m2i = data_depression$comp_mean, sd1i = data_depression$int_SD, sd2i = data_depression$comp_SD, n1i = data_depression$int_n, n2i = data_depression$comp_n)
depression <- cbind(data_depression, SMD)
sei <- sqrt(depression$vi/(depression$int_n + depression$comp_n))
depression <- cbind(depression, sei)
skewness(depression$yi) 
skewness(depression$yi[depression$Alc_category == "Included"])
skewness(depression$yi[depression$Alc_category == "Excluded"])
metaanalysis <- metagen(TE = depression$yi, seTE = depression$sei, studlab = depression$Citation, sm = "SMD", method.tau = "DL", byvar = depression$Alc_category)
## Sensitivity analyses
measures_SMD <- escalc(measure = "SMD", m1i = data_depression$otherint_mean, m2i = data_depression$othercomp_mean, sd1i = data_depression$otherint_SD, sd2i = data_depression$othercomp_SD, n1i = data_depression$int_n, n2i = data_depression$comp_n)
measures_depression <- cbind(depression, measures_SMD)
colnames(measures_depression) <- make.unique(colnames(measures_depression))
measures_sei <- sqrt(measures_depression$vi.1/(measures_depression$int_n + measures_depression$comp_n))
measures_depression <- cbind(measures_depression, measures_sei)
measures_MA <- metagen(TE = measures_depression$yi.1, seTE = measures_depression$sei, byvar = measures_depression$Alc_category,  method.tau = "DL")
duration_SMD <- measures_SMD <- escalc(measure = "SMD", m1i = data_depression$earlyint_mean, m2i = data_depression$earlycomp_mean, sd1i = data_depression$earlyint_SD, sd2i = data_depression$earlycomp_SD, n1i = data_depression$int_n, n2i = data_depression$comp_n)
duration_depression <- cbind(data_depression, duration_SMD)
duration_sei <- sqrt(duration_depression$vi/(duration_depression$int_n + duration_depression$comp_n))
duration_depression <- cbind(duration_depression, duration_sei)
duration_MA <- metagen(TE = duration_depression$yi, seTE = duration_depression$duration_sei, byvar = duration_depression$Alc_category, method.tau = "DL")
metagen_overall <- metagen(TE = depression$yi, seTE = depression$sei, studlab = depression$Citation, method.tau = "DL")
egger <- metabias(metagen_overall, method.bias = "Egger")
influence <- metainf(metaanalysis, pooled = "random")
plot(influence)
influence$seTE
influence_depression <- depression[-6,]
influence_metagen <- metagen(TE = influence_depression$yi, seTE = influence_depression$sei, studlab = influence_depression$Citation, method.tau = "DL", byvar = influence_depression$Alc_category)
substance_SMD <- escalc(measure = "SMD", m1i = substance_depression$int_mean, m2i = substance_depression$comp_mean, sd1i = substance_depression$int_SD, sd2i = substance_depression$comp_SD, n1i = substance_depression$int_n, n2i = substance_depression$comp_n)
substance_depression <- cbind(substance_depression, substance_SMD)
substance_sei <- sqrt(substance_depression$vi/(substance_depression$int_n + substance_depression$comp_n))
substance_depression <- cbind(substance_depression, substance_sei)
substance_metaanalysis <- metagen(TE = substance_depression$yi, seTE = substance_depression$substance_sei, studlab = substance_depression$Citation, sm = "SMD", method.tau = "DL", byvar = substance_depression$Alc_category)

## 3. Relapse rates. 
data_relapse$Alc_category <- factor(data_relapse$Alc_category)
RR <- escalc(measure = "RR", ai = data_relapse$int_relapse, bi = data_relapse$int_no_relapse, ci = data_relapse$comp_relapse, di = data_relapse$comp_no_relapse, data = data_relapse)
relapse <- cbind(data_relapse, RR)
relapse_sei <- sqrt(relapse$vi/(relapse$int_n + relapse$comp_n))
relapse <- cbind(relapse, relapse_sei)
skewness(relapse$yi)
skewness(relapse$yi[relapse$Alc_category == “Included”])
skewness(relapse$yi[relapse$Alc_category == “Excluded”])
relapse_MA <- metagen(TE = relapse$yi, seTE = relapse$relapse_sei, byvar = relapse$Alc_category)

## 4. Dropout. 
overall_dropout <- metaprop(event = data_dropout$total_dropout, n = data_dropout$total_n, byvar = data_dropout$Alc_category, method = "inverse", method.tau = "DL")
RR <- escalc(measure = "RR", ai = data_dropout$int_dropout, bi = data_dropout$int_no_dropout, ci = data_dropout$comp_dropout, di = data_dropout$comp_no_dropout, data = data_dropout)
dropout <- cbind(data_dropout, RR)
skewness(dropout$yi)
skewness(dropout$yi[dropout$Alc_category == “Included”])
skewness(dropout$yi[dropout$Alc_category == “Excluded”])
differential_dropout <- metagen(TE = dropout$yi, seTE = dropout$vi, byvar = dropout$Alc_category, method.tau = "DL")

