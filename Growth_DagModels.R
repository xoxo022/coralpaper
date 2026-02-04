##----------STATISTICAL DATA ANALYSES---------------
rm(list=ls())

library(tidyverse)
library(rstanarm)
library(brms)
library(tidybayes)
library(parallel)
library(here)

##---------------------------------------------
##Model Growth rates probabilities
##---------------------------------------------
myplots <- paste(here(), "Plots", sep="/")

load("growth_data.Rdata")


num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

rawdata <- list(enviro.vars, seasonal.data)

df <- seasonal.data

names(df)[names(df) == "isize"] <- "size"
names(df)[names(df) == "Max.Diameter"] <- "growth"
names(df)[names(df) == "Macroalgae"] <- "macroalgae"
names(df)[names(df) == "Coral"] <- "coralcov"
names(df)[names(df) == "Algae overgrown"] <- "overgrowth"
names(df)[names(df) == "BGratio"] <- "Turbidity"
names(df)[names(df) == "NoonPAR"] <- "Light"
names(df)[names(df) == "maxTemp"] <- "Temp"
names(df)[grep("WaveHeight", names(df))] <- "waves"
names(df)[grep("CurSpeed", names(df))] <- "currents"
names(df)[grep("WindSpeed", names(df))] <- "wind"
names(df)[grep("tide", names(df))] <- "Tides"

df %<>% na.omit()

## Transform variables
transformed.variables <- df[,!names(df) %in% c("Survey","SiteQ","Site","UniqueID","size","Date","Ctype")]
environ.vars <- c("Light","Temp","waves","currents","wind","Tides")

transformed.variables <- transformed.variables %>%
  mutate(across(all_of(environ.vars), scale)) %>%
  mutate(
    growth = log(growth),
    overgrowth = as.integer(overgrowth)
  )

## Mutate waves to categorical
transformed.variables$waves[transformed.variables$waves < 0] <- "low"
transformed.variables$waves[!transformed.variables$waves=="low"] <- "High"
transformed.variables$waves <- as.integer(factor(transformed.variables$waves, levels=c("low","High")))

##-----------------------------------------------------------------------------
## FINAL MODELS
##-----------------------------------------------------------------------------
final.structure <- df[, c("Site","SiteQ","UniqueID","size","Days","growth","Ctype")]
final.structure$size <- log(final.structure$size)
final.structure$Days <- scale(final.structure$Days, center=T, scale=T)

finalmodel <- cbind(final.structure, transformed.variables[, c("macroalgae","overgrowth","Turf","Light","Temp","currents","waves","wind","Tides")])

source("coral_biophysical_dag.R")
mydag <- growth_dag

## Extract DAG-informed covariates
temp.covars <- c("Temp", adjustmentSets(mydag,"Temp","growth")[[1]])
Light.covars <- c("Light", adjustmentSets(mydag,"Light","growth")[[1]])
current.covars <- c("currents", adjustmentSets(mydag,"currents","growth")[[1]])
tides.covars <- c("Tides", adjustmentSets(mydag,"Tides","growth")[[1]])
Turf.covars <- c("Turf", adjustmentSets(mydag,"Turf","growth")[[1]])
overgrowth.covars <- c("overgrowth", adjustmentSets(mydag,"overgrowth","growth")[[2]])

mymodels <- list(
  te = temp.covars,
  li = Light.covars,
  cu = current.covars,
  tu = Turf.covars,
  ov = overgrowth.covars,
  ti = tides.covars
)

##----------------------
## Fit bayesian models
##---------------------- 
rstan_options(auto_write = TRUE)
options(mc.cores = detectCores())

model_results <- list()
posterior_summaries <- list()

for (model_name in names(mymodels)) {
  predictors <- mymodels[[model_name]]
  formula_str <- paste("growth ~", paste(predictors, collapse = " + "), "+ size + Days + (1 |SiteQ/UniqueID)")
  model_formula <- as.formula(formula_str)
  
  fit <- brm(
    formula = model_formula,
    data = finalmodel,
    family = Gamma(link="log"),
    prior = c(
      prior(normal(0, 1.5), class = "b"),
      prior(normal(0, 1), class = "sd"), 
      prior(gamma(4, 0.5), class = "shape")
    ), 
    sample_prior = TRUE,
    chains = 4, cores = 4, iter = 5000, warmup = 1000,
    control = list(adapt_delta = 0.999, max_treedepth = 20),
    silent = TRUE, refresh = 0
  )
  
  model_results[[model_name]] <- fit
  
  draws <- as_draws_df(fit)
  mypred <- names(draws)[2]
  
  coef.summary <- tibble(
    term = mypred,
    estimate = mean(draws[[mypred]]),
    lower_95 = quantile(draws[[mypred]], 0.025),
    upper_95 = quantile(draws[[mypred]], 0.975),
    lower_50 = quantile(draws[[mypred]], 0.25),
    upper_50 = quantile(draws[[mypred]], 0.75),
    prob.dir = mean(draws[[mypred]] < 0)
  ) %>% mutate(model = model_name)
  
  posterior_summaries[[paste(model_name, mypred, sep="_")]] <- coef.summary
}

## Combine results
all_estimates <- bind_rows(posterior_summaries)
all_estimates$term <- factor(tools::toTitleCase(gsub("^b_","",all_estimates$term)), 
                             levels=c("Overgrowth","Turf","Temp","Light","Currents","Tides"))
all_estimates <- all_estimates %>%
  mutate(
    effect_type = case_when(
      prob.dir > 0.975 ~ "negative",
      prob.dir < 0.025 ~ "positive",
      TRUE ~ "uncertain"
    ),
    percentchange = (exp(estimate) - 1) * 100
  )

save(rawdata, conditional.test, mymodels,finalmodel,model_results, 
     posterior_summaries, all_estimates, 
     file = paste(df_name,"Growth_models.RData", sep="_"))

##------------------------
gby<-ggplot(all_estimates, aes(y = term, x=estimate, colour=prob.dir)) +
  geom_pointrange(aes(xmin = lower_50, xmax = upper_50), linewidth = 0.8, color = "gray20") +
  geom_pointrange(aes(xmin = lower_95, xmax = upper_50), linewidth = 0.2, color = "gray40") +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradient2(
    low = "#0072B2",     # blue (low P < 0 = positive)
    mid = "gray70",      # uncertain
    high = "#D90E00",    # red (high P < 0 = negative)
    midpoint = 0.5,
    limits = c(0, 1)
  ) +
  labs(x = "Posterior mean", y = "")+ggtitle("Growth")+
  theme_minimal()+theme(legend.position = "none", title = element_text(size=8))+
  theme(
    panel.grid.major = element_line(colour = "grey98"),
    panel.grid.minor = element_line(colour = "grey98")
  )


gby<-gby + labs(tag = "b") 

ggsave(gby,filename =paste(myplots,"/",df_name,"_growth_model_estimates.png",sep = ""), dpi=300, width = 3, height = 3.5)

