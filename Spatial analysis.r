
## This script has been used for the paper `https://pubmed.ncbi.nlm.nih.gov/30989204/`
## The script has two sections: ANOVA, CART 

## ANOVA ## -----------

# Read data
my_data <- read.csv(file.choose("E:/PhD project file/10_04_2018/woreda_128/cases/Output/LISAGIS.csv"))
my_data <- as.data.frame(my_data)
my_data

# Show a random sample 
set.seed(1234)
dplyr::sample_n(my_data,10)

# Show the levels 
levels(my_data$LISA)
my_data$LISA <- ordered(my_data$LISA,
                        levels = c("HH", "outlier", "Notsig"))
levels(my_data$LISA)

# Compute summary statistics by group

library(dplyr)

group_by(my_data,LISA)%>%
  summarise(
    count = n(),
    mean = mean(propallcase,na.rm = TRUE),
    sd=sd(propallcase,na.rm = TRUE)
  )

# Visualize your data with ggpubr:
# Box plots
# Plot weight by group and color by group

library("ggpubr")
ggboxplot(my_data, x = "LISA", y = "propallcase", 
          color = "LISA", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("HH", "outlier", "Notsig"),
          ylab = "propallcase", xlab = "LISA")

# Mean plots
# Plot weight by group
# Add error bars: mean_se
# (other values include: mean_sd, mean_ci, median_iqr, ....)

library("ggpubr") 

ggline(my_data, x = "LISA", y = "propallcase", 
       add = c("mean_se", "jitter"), 
       order = c("HH", "outlier", "Notsig"),
       ylab = "propallcase", xlab = "LISA")

# Box plot
boxplot(weight ~ group, data = my_data,
        xlab = "Treatment", ylab = "Weight",
        frame = FALSE, col = c("#00AFBB", "#E7B800", "#FC4E07"))
# plot means

library("gplots")
plotmeans(propallcase ~ LISA, data = my_data, frame = FALSE,
          xlab = "LISA", ylab = "propallcase",
          main="Mean Plot with 95% CI") 

# Compute the analysis of variance

res.aov <- aov(propallcase ~ LISA, data = my_data)

# Summary of the analysis

summary(res.aov)

# Tukey multiple pairwise-comparisons

TukeyHSD(res.aov)

# Multiple comparision 

library(multcomp)
summary(glht(res.aov, linfct = mcp(LISA = "Tukey")))

# Pairwise t-test.Here, the p-values have been adjusted by the Benjamini-Hochberg method.

pairwise.t.test(my_data$propallcase, my_data$LISA,
                p.adjust.method = "BH")

# Check ANOVA assumptions: test validity?

# 1. Homogeneity of variances

plot(res.aov, 1)

# Levene's test using car package

library(car)
leveneTest(propallcase ~ LISA, data = my_data)

# ANOVA test with no assumption of equal variances. 
# Welch one-way test), that does not require that assumption have been implemented in the function oneway.test().

oneway.test(propallcase ~ LISA, data = my_data)

# Pairwise t-tests with no assumption of equal variances

pairwise.t.test(my_data$propallcase, my_data$LISA,
                p.adjust.method = "BH", pool.sd = FALSE)
# Check the normality assumptions
# 2. Normality

plot(res.aov, 2)

# To support the graphical presentation 
# Extract the residuals

aov_residuals <- residuals(object = res.aov )

# Run Shapiro-Wilk test

shapiro.test(x = aov_residuals )


## CART ## --------

##Regression tree

library(rpart)
# tree
Tuberculosis
fit <- rpart(propallcase~Immigrant5yr+popdensity+urbanResi+roomscrowding+Housecrowding+Popcharcoal+Popwood+Popdung+noteconomicactive+unemployment+noschooling+PLHIV+PLHIVART+GraduteHEP+COverageHC_Ho+Percentmale+aveelva, method = "anova",data = Tuberculosis)
printcp(fit)# display the result
plotcp(fit)# visualize cross-validation results 
summary(fit)
# plot tree
plot(fit, uniform=TRUE,main="Regression Tree for TB")
# create additional plots
par(mfrow=c(1,2)) # two plots on one page
rsq.rpart(fit) # visualize cross-validation results

text(fit, use.n=TRUE, all=TRUE, cex=.8)
# create attractive postscript plot of tree
post(fit, file = "C:/Users/s4415062/Documents/TB_document/Final data_analysis/Spatial analysis _P1/tree.ps",title = "Regression Tree for TB")
# prune the tree
pfit<- prune(fit, cp=   fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])

# plot the pruned tree
plot(pfit, uniform=TRUE,main="Pruned Classification Tree for TB")
text(pfit, use.n=TRUE, all=TRUE, cex=.8)
post(pfit, file = "C:/Users/s4415062/Documents/TB_document/Final data_analysis/Spatial analysis _P1/ptree.ps",title = "Pruned Classification Tree for Kyphosis")

# Conditional Inference Tree for Kyphosis
library(party)
fit2 <- ctree(propallcase~Immigrant5yr+popdensity+urbanResi+roomscrowding+Housecrowding+Popcharcoal+Popwood+Popdung+noteconomicactive+unemployment+noschooling+PLHIV+PLHIVART+GraduteHEP+COverageHC_Ho+Percentmale+aveelva,data = Tuberculosis)
plot(fit2, main="Conditional Inference Tree for TB")

# Random Forest prediction of Kyphosis data
library(randomForest)
fit3 <- randomForest(propallcase~Immigrant5yr+popdensity+urbanResi+roomscrowding+Housecrowding+Popcharcoal+Popwood+Popdung+noteconomicactive+unemployment+noschooling+PLHIV+PLHIVART+GraduteHEP+COverageHC_Ho+Percentmale+aveelva,data = Tuberculosis)
print(fit3) # view results
importance(fit3) # importance of each predictor 

