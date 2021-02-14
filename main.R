# Author: Sharder Islam
# Document: A Survival Analysis of Time to Wean
# Course: Applied Survival Analysis Spring 2020
# Institution: NYU School of Global Public Health

# ---Prepare dataset---------------------------------------------------------- #
library(survival)
library(survminer)

df <- read.csv("breastfeed.csv")
str(df)

df[,c(3:6, 10)] <- lapply(df[,c(3:6, 10)], factor)
df$age_c <-cut(df$age, breaks=2)
df$birthyr_c <-cut(df$birthyr, breaks=2)
df$educ_c <-cut(df$educ, breaks=2)

# ---Examine ties and univariates--------------------------------------------- #
summary(df)
table(df$length, df$complete)

# ---Create Kaplan-Meier plot------------------------------------------------- #
km <- survfit(Surv(length, complete) ~1, data=df, conf.type="log-log")
summary(km)
ggsurvplot(km,
           legend="none",
           xlab="# Weeks",
           ylab="Survival probability",
           conf.int = TRUE,
           surv.median.line = "hv")

# ---Kaplan Meier by covariate------------------------------------------------ #
par(mfrow=c(4, 2))
plot(survfit(Surv(df$length, df$complete) ~ df$race), main="race")
plot(survfit(Surv(df$length, df$complete) ~ df$poverty), main="poverty")
plot(survfit(Surv(df$length, df$complete) ~ df$smoke), main="smoke")
plot(survfit(Surv(df$length, df$complete) ~ df$alcohol), main="alcohol")
plot(survfit(Surv(df$length, df$complete) ~ df$prenatal3), main="prenatal3")
plot(survfit(Surv(df$length, df$complete) ~ df$age_c), main="age_c")
plot(survfit(Surv(df$length, df$complete) ~ df$birthyr_c), main="birthyr_c")
plot(survfit(Surv(df$length, df$complete) ~ df$educ_c), main="educ_c")

# ---Checking distribution for time to wean----------------------------------- #
par(mfrow=c(1, 2))
plot(km,
     fun='cumhaz',
     main='Cumulative hazard over time',
     xlab='Time',
     ylab='Cumulative hazard')
plot(km,
     fun='cloglog',
     main="Log - Log Hazard over log time",
     xlab='log time',
     ylab = 'Log cumulative hazard')

# ---Performing two-group comparisons----------------------------------------- #
for (i in colnames(df)[c(3:13)]){
    print(i)
    col <- df[, i]
    print(survdiff(formula = Surv(df$length, df$complete) ~ col, rho=0))
    print(survdiff(formula = Surv(df$length, df$complete) ~ col, rho=1))
    cat("\n\n")
}

# ---Inclusion into base model------------------------------------------------ #
for (i in colnames(df)[c(4, 6, 7, 10)]){
    var <- df[, i]
    mod <- coxph(Surv(df$length, df$complete)
        ~ df$race + df$smoke + df$birthyr + df$educ + var, ties="exact")
    if (coef(summary(mod))[30] <= 0.2){
        print(i)
        print(mod)
        cat("\n")
    }
}

# ---Create models------------------------------------------------------------ #
mod0 <- coxph(Surv(length, complete) ~ 1, data = df, ties="exact")
mod1 <- coxph(Surv(length, complete)
     ~ race + smoke + birthyr + educ, data = df, ties="exact")
mod2 <- coxph(Surv(length, complete)
     ~ race + smoke + birthyr + educ + poverty, data = df, ties="exact")
mod3 <- coxph(Surv(length, complete)
     ~ race + smoke + birthyr + educ * poverty,  data = df, ties="exact")
mod4 <- coxph(Surv(length, complete)
     ~ race + smoke + birthyr + strata(educ), data = df, ties="exact")

# ---Fit statistics---------------------------------------------------------- #
-2*logLik(mod0)[1]
-2*logLik(mod1)[1]
-2*logLik(mod2)[1]
-2*logLik(mod3)[1]
-2*logLik(mod4)[1]
AIC(mod0)
AIC(mod1)
AIC(mod2)
AIC(mod3)
AIC(mod4)

# ---Test proportional hazards assumption------------------------------------- #
cox.zph(mod3, terms=TRUE, global=FALSE)
cox.zph(mod4, terms=TRUE, global=FALSE)

# ---Testing linearity of birthyr variable-------------------------------------#
modby1 <- coxph(Surv(length, complete)
       ~ smoke + birthyr, data = df, ties="exact")
modby2 <- coxph(Surv(length, complete)
       ~ smoke + birthyr_c , data = df, ties="exact")
modby3 <- coxph(Surv(length, complete)
       ~ smoke + + birthyr + birthyr_c , data = df, ties="exact")

anova(modby1, modby3)
anova(modby2, modby3)

# ---Check for influential observations----------------------------------------#
ggcoxdiagnostics(mod4, type = 'deviance', linear.predictions = FALSE)