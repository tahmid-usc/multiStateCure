# Simulate multistate cure survival data
# Reference: Beesley, L. J., & Taylor, J. M. (2019). EM algorithms for fitting multistate cure models. Biostatistics, 20(3), 416-432.
# Situation 01: no covariate missingness or unequal follow-up


# Simulation 01 

library(mvtnorm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survminer)
library(survival)
library(ggfortify)
library(gridExtra)
set.seed(2021)


#devtools::install_github("lbeesleyBIOSTAT/MultiCure", build = TRUE, build_opts = c()) 
library(MultiCure) 
#utils::browseVignettes('MultiCure')


n <- 2000
id <- 1:2000

# X1, X2 from multivariate normal var = 1, corr = .5
X <- rmvnorm(n, rep(0, 2), cbind(c(1,.5), c(.5,1)))

#cure status: expit(P(Gi = 1jX)) = 0.5 + 0.5X1 + 0.5X2

prob <- 1 / (1 + exp(-(0.5 + 0.5 * X[,1] + 0.5 * X[,2])))
G <- rbinom(n, 1, prob = prob)

#Survival times

#For cured subjects, death time using a proportional hazards model with H_0_24(t) = 0.002t^1.4

a <- 0.002; b <- 1.4
T_24 <- (-log(runif(n)) * exp(-0.5*(rowSums(X))) / a)^(1/b) 

T_13 <- 10^6
T_14 <- 10^6
T_34 <- 10^6

C <- runif(n, 10, 80) # censoring time

# 1 --> 4; Cured => G = 0

dt <- data.frame(id = id, T_13 = T_13, T_14 = T_14, T_24 = T_24, T_34 = T_34, C = C, x1 = X[,1], x2 = X[,2], G = G)



# For non-cured subjects
# Recurrence              H_0_13(t) = 0.005t^2 
# Death from other causes H_0_14(t) = 0.002t^1.4
# Death after reccurernce H_0_34(t) = 0.08t^1.9

a <- 0.005; b <- 2
T_13 <- (-log(runif(n)) * exp(-0.5*(rowSums(X))) / a)^(1/b) 

a <- 0.002; b <- 1.4
T_14 <- (-log(runif(n)) * exp(-0.5*(rowSums(X))) / a)^(1/b) 

T_24 <- 10^6

a <- 0.08; b <- 1.9
T_34 <- (-log(runif(n)) * exp(-0.5*(rowSums(X))) / a)^(1/b) 


#replace times in the dataframe
dt[G == 1,]$T_13 <- T_13[G == 1]
dt[G == 1,]$T_14 <- T_14[G == 1]
dt[G == 1,]$T_34 <- T_34[G == 1]
dt[G == 1,]$T_24 <- T_24[G == 1]



#create variables
#for cured
dt <- dt %>% mutate(Y_d = pmin(T_24, C), delta_d = as.numeric(T_24 < C))

# for uncured; death before reccurence 
dt <- dt %>% mutate(Y_d = if_else(G == 1 & T_14 < T_13, pmin(T_14, C), Y_d),
                    Y_r = Y_d,
                    delta_d = if_else(G == 1 & T_14 < T_13, as.numeric(T_14 < C), delta_d),
                    delta_r = 0)

# uncured, death after recurrence
dt <- dt %>% mutate(
  Y_r =     if_else(G == 1 & T_14 > T_13, pmin(T_13, C), Y_r),
  Y_d =     if_else(G == 1 & T_14 > T_13, pmin(T_13 + T_34, C), Y_d),
  delta_r = if_else(G == 1 & T_14 > T_13, as.numeric(T_13 < C), delta_r),
  delta_d = if_else(G == 1 & T_14 > T_13, as.numeric(T_13 + T_34 < C), delta_d)
  )                 


none <- dt %>%  select(id, Y_d, Y_r, C, delta_d, delta_r, x1, x2, G)

NONE = SimulateMultiCure(type = 'NoMissingness') 

summary(NONE); summary(none)


#------------------------------------------

vis_none <- none %>% 
  mutate(event = ifelse(delta_r == 1, 'Dead', 'Alive')) %>%
  mutate(event = ifelse(delta_d == 1, 'Recurrence', event)) %>% 
  mutate(event = ifelse(delta_r == 0 & delta_d == 0, 'Censored', event))

vis_none %>% 
  ggplot(aes(x = id, y = Y_d)) +
  geom_linerange(aes(ymin = 0, ymax = Y_d)) +
  geom_point(aes(shape = event, color = event), stroke = 1, cex = 2) +
  scale_shape_manual(values = c(1, 3, 4)) +
  labs(y = "Time (years)", x = "Subject ID") + coord_flip() + theme_classic()

  
# Basic survival curves

par(mfrow = c(2,2))

fit_km <- survfit(Surv(Y_d, delta_d) ~ 1, data = none)
p1 <- ggsurvplot(fit_km, risk.table = TRUE, xlab = "Time (years)", censor = T, title = 'Overall survival')


fit_km <- survfit(Surv(Y_D, delta_D) ~ 1, data = NONE)
p2 <- ggsurvplot(fit_km, risk.table = TRUE, xlab = "Time (years)", censor = T, title = 'Overall survival (multicure)')

arrange_ggsurvplots(list(p1, p2))
