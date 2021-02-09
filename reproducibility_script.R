## ----setup, include=FALSE---------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)


## ---------------------------------------------------------------------------------------------------------
library(survival)
library(coxme)
library(frailtyEM)
library(tidyverse)


## ---------------------------------------------------------------------------------------------------------
dat <- read.csv("EORTC_10854.csv", stringsAsFactors = FALSE)

dat$periop <- 
  factor(dat$periop, 
         levels=c("no periop chemo","periop chemo"))
dat$surgery <- 
  factor(dat$surgery, 
         levels=c("mastectomy with RT","mastectomy without RT","breast conserving"))
dat$tusi <- 
  factor(dat$tusi, 
         levels=c("<2 cm","2-5 cm",">5 cm"))
dat$nodal <- 
  factor(dat$nodal, 
         levels=c("node negative","node positive"))
dat$age50 <- 
  factor(dat$age50, 
         levels=c("<=50",">50"))
dat$adjchem <- 
  factor(dat$adjchem, 
         levels=c("no adj chemo","adj chemo"))
dat$tam <- 
  factor(dat$tam, 
         levels=c("no tam","tam"))

dat$hospno <- factor(dat$hospno)



## ---------------------------------------------------------------------------------------------------------
dat %>% 
  group_by(hospno) %>% 
  summarize(npat = n()) %>% 
  ggplot(aes(x = npat)) + geom_histogram() + 
  labs(x = "Center size") +
  # scale_x_continuous(breaks = c(0, 50, 100, 182, 304, 602, 880)) +
  scale_x_continuous(breaks = seq(0, 900, by=50)) +
  theme_bw()


## ---------------------------------------------------------------------------------------------------------
# data set with Kaplan-Meier curves per center
sf <- survfit(Surv(survyrs, survstat) ~ strata(hospno), dat)
dd <- data.frame(time = sf$time, surv = sf$surv, center = as.factor(rep(1:length(sf$strata), sf$strata))) 
startpoint <- data.frame(center = unique(dd$center), time = 0, surv = 1)
dd <- rbind(dd, startpoint)

# data set with overall Kaplan-Meier curve
sf2 <- survfit(Surv(survyrs, survstat) ~ 1, dat)
dd2 <- data.frame(time = sf2$time, surv = sf2$surv) 

dd %>% 
  ggplot(aes(x = time, y = surv)) + 
  geom_step(aes(group = center),  show.legend = FALSE, alpha = 0.2) + 
  geom_step(data = dd2) + 
  theme_classic() + 
  labs(x = "Years since randomisation", y = "Survival")


## ---------------------------------------------------------------------------------------------------------
m_cph <- coxph(Surv(survyrs, survstat) ~ surgery + tusi + nodal + age50 + adjchem +
        tam + periop + cluster(hospno), dat, ties = "breslow")
m_cph

m_fixed <- coxph(Surv(survyrs, survstat) ~ surgery + tusi + nodal + age50 +
                 adjchem + tam + periop + hospno, dat, ties = "breslow")

m_fixed
m_fr_cov <- coxph(Surv(survyrs, survstat) ~ frailty(hospno) + surgery + tusi + nodal + age50 +
                    adjchem + tam + periop, data = dat, ties = "breslow")
m_fr_cov


m_emfrail <- emfrail(Surv(survyrs, survstat) ~ cluster(hospno) + surgery + tusi + nodal + age50 +
                       adjchem + tam + periop,
                     data = dat)
summary(m_emfrail)
# -5450.928 likelihood

m_fr_cov_lnorm <- coxme(Surv(survyrs, survstat) ~ (1|hospno) + surgery + tusi + nodal + age50 +
                          adjchem + tam + periop, data = dat, ties = "breslow")

m_fr_cov_lnorm
# -5449.445 


## ---------------------------------------------------------------------------------------------------------
# Obtain frailty estimates
sm <- summary(m_emfrail)
frailties <- sm$frail

# Obtain center estimates and calculate standard error
centers <- grep("hospno", names(m_fixed$coefficients))
nc <- length(centers)
var_centers <- m_fixed$var[centers, centers]

varmean <- sum(var_centers)/ nc^2

newse <- sqrt(c(varmean, diag(var_centers) + varmean - 2/nc * apply(var_centers, 1, sum)))

# Arrange frailty estimates and fixed effects estimates
logz <- frailties %>% 
  mutate_if(is.numeric, log) %>% 
  rename(estimate = z, ymin = lower_q, ymax = upper_q) %>% 
  mutate(type = "frailty")
  
fe <- data.frame(beta = c(hospnoA = 0, m_fixed$coefficients[centers]), sd = newse) %>% 
  rownames_to_column("Center") %>% 
  mutate(Center = substr(Center, 7, 9)) %>% 
  mutate(sbeta = sum(beta) / 15) %>% 
  mutate(beta = beta - sbeta) %>% 
  arrange(beta) %>%
  mutate(ord = 1:n()) %>% 
  mutate(Center = as.character(Center)) %>% 
  mutate(ymin = beta - 1.96 * sd, ymax = beta + 1.96 * sd) %>% 
  rename(id = Center, estimate = beta) %>% 
  select(-sd, -sbeta, -ord) %>% 
  mutate(type = "fixed effects") %>% 
  mutate(ord = 1:n()) 

ord_fe <- data.frame(id = fe$id, ord = fe$ord)
logz <- logz %>% 
  right_join(ord_fe)

datt <- bind_rows(logz, fe) 

# Finally, the plot
datt %>% 
 mutate_at(vars(estimate, ymin, ymax), .funs = exp) %>% # make it exp()
  ggplot(aes(x = ord, y = estimate)) + 
  geom_point(aes(colour = type), 
             position = position_dodge(.9)) + 
  geom_errorbar(aes(ymin = ymin, ymax = ymax, colour = type, linetype = type),
                position = position_dodge()) + 
  scale_x_continuous(labels = as.character(datt$id)[1:15], 
                     breaks = seq_along(datt$estimate)[1:15]) + 
  scale_y_continuous(trans = "log", 
                     breaks = seq(from = 0.5, to = 2.5, by = .5)) + 
  labs(x = "Center", 
       y="Center effect (hazard ratio)") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +  
  guides(colour = guide_legend(title = element_blank()), linetype = guide_legend(title = element_blank())) + 
  scale_colour_brewer(palette = "Set1")


## ---------------------------------------------------------------------------------------------------------
# checking proportional hazards assumption

# in Cox model
cox.zph(m_cph)

# in Cox model with log-frailties as offset
logz_long <- log(m_emfrail$frail)[dat$hospno]

m_cph_offset_gamma <- coxph(Surv(survyrs, survstat) ~ surgery + tusi + nodal + age50 + adjchem +
                 tam + periop + offset(logz_long), dat, ties = "breslow")
cox.zph(m_cph_offset_gamma)



## ---------------------------------------------------------------------------------------------------------
data(cgd)
# head(cgd)

# rebrand id
cgd <- cgd %>% 
  mutate(id = as.numeric(as.factor(id))) 


## ---------------------------------------------------------------------------------------------------------
cgd %>% 
  ggplot() +
  geom_segment(aes(x = tstart, xend = tstop, y = id, yend = id)) + 
  geom_point(data = filter(cgd, status == 1), aes(x = tstop, y = id)) + 
  theme_bw() + 
  labs(x = "time", y = "patient")




## ---------------------------------------------------------------------------------------------------------
cgd %>% 
  group_by(id) %>% 
  summarize(n_events = sum(status)) %>% 
  ggplot(aes(x = n_events)) + 
  geom_bar() + 
  labs(x= "Number of events / individual") + 
  theme_bw()


## ---------------------------------------------------------------------------------------------------------
mvals <- seq(from = 0.1, to = 2, by = 0.2)
models <- lapply(mvals, function(x) emfrail(formula = Surv(tstart, tstop, status) ~ sex +
    treat + age + propylac + inherit + steroids + cluster(id), data = cgd, 
    distribution = emfrail_dist(dist = "pvf", pvfm = x)) )

likelihoods <- sapply(models, function(x) x$loglik[2])

mvals[which.max(likelihoods)]



## ---------------------------------------------------------------------------------------------------------
m_cph <- coxph(Surv(tstart, tstop, status) ~ treat + age +inherit + steroids, cgd)
# summary(m_cph)


m_gamma_emf <- emfrail(Surv(tstart, tstop, status) ~ treat + sex + age +inherit + steroids + cluster(id), cgd)
# summary(m_gamma_emf)

mod_ig <- emfrail(formula = Surv(tstart, tstop, status) ~  treat + sex +
                    age + inherit + steroids + cluster(id), 
                  distribution = emfrail_dist(dist = "pvf"),
                  data = cgd)

mod_cp <- emfrail(formula = Surv(tstart, tstop, status) ~  treat + sex +
                    age + inherit + steroids + cluster(id), 
                  distribution = emfrail_dist(dist = "pvf", pvfm = 0.5),
                  data = cgd)

mod_cp_11 <- emfrail(formula = Surv(tstart, tstop, status) ~  treat + sex +
                    age + inherit + steroids + cluster(id), 
                  distribution = emfrail_dist(dist = "pvf", pvfm = 1.1),
                  data = cgd)

mod_stab <- emfrail(formula = Surv(tstart, tstop, status) ~ treat + sex + 
                      age + inherit + steroids + cluster(id), 
                    distribution = emfrail_dist(dist = "stable"),
                    data = cgd)



## ---------------------------------------------------------------------------------------------------------
autoplot(m_gamma_emf, type = "hist") + 
  theme_bw()


## ---------------------------------------------------------------------------------------------------------
sessionInfo()

