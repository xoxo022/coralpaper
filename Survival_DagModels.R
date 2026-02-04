##===========================================================
##   CORAL SURVIVAL ANALYSIS — DAG-INFORMED BAYESIAN MODELS
##===========================================================

rm(list = ls())

##------------------
## Load libraries
##------------------
library(tidyverse)   # Data manipulation, pipes, plotting
library(here)        # File path handling
library(dagitty)     # DAG creation and analysis
library(ggdag)       # DAG plotting
library(brms)        # Bayesian modeling
library(parallel)    # Detect CPU cores
library(rstan)       # Backend for brms
library(ggplot2)     # Plotting and ggsave
library(tools)       # toTitleCase()

##------------------
## Setup
##------------------
myplots <- file.path(here(), "Plots")

load("survival_data.Rdata")

num_cores <- detectCores() - 1


df_name <- "survival.data"

rawdata <- list(enviro.vars, survival.data)
df <- survival.data

## Rename variables to DAG equivalents
df <- df %>%
  rename(
    size = isize,
    macroalgae = Macroalgae,
    overgrowth = `Algae overgrown`,
    Turbidity = BGratio,
    Light = NoonPAR,
    Temp = maxTemp
  )

names(df)[grep("WaveHeight", names(df))] <- "waves"
names(df)[grep("CurSpeed", names(df))]   <- "currents"
names(df)[grep("WindSpeed", names(df))]  <- "wind"
names(df)[grep("tide", names(df))]       <- "Tides"

df <- na.omit(df)

##--------------------------
## Scale and transform data
##--------------------------
apply_transformation <- function(data, cols, transformation) {
  data %>% mutate(across(all_of(cols), transformation))
}

environ.vars <- c("Light", "Temp", "waves", "currents", "wind", "Tides")

transformed.variables <- df %>%ungroup() %>%
  select(-Survey, -Site, -SiteQ, -Ctype, -UniqueID, -size)

## Scale environmental variables
transformed.variables <- apply_transformation(transformed.variables, environ.vars, scale)

## Convert binary to integer
transformed.variables$overgrowth <- as.integer(transformed.variables$overgrowth)
transformed.variables$survival   <- as.integer(transformed.variables$survival)

##---------------------------------------------
## Test conditional independencies (diagnostic)
##---------------------------------------------

source("coral_biophysical_dag.R")

conditional.test <- localTests(survival_dag, transformed.variables, max.conditioning.variables = 4)
conditional.test$p.value <- p.adjust(conditional.test$p.value)

png(filename = file.path(myplots, paste0(df_name, "_dag_survival_conditionals.png")),
    width = 900, height = 700)
plotLocalTestResults(localTests(survival_dag, transformed.variables, abbreviate.names = FALSE))
title("Survival model")
dev.off()

##---------------------------------------------
## "Categorical" waves variable
##---------------------------------------------
transformed.variables$waves[transformed.variables$waves < 0] <- "low"
transformed.variables$waves[!transformed.variables$waves == "low"] <- "High"
transformed.variables$waves <- ordered(transformed.variables$waves, levels = c("low", "High"))
transformed.variables$waves <- as.integer(transformed.variables$waves)

##---------------------------------------------
## Final model dataset
##---------------------------------------------
final.structure<-df[, c("Site","SiteQ","UniqueID","size","Days","Ctype", "survival")]
final.structure$Days<-as.numeric(scale(final.structure$Days, center=T, scale=T))
final.structure$size<-log(final.structure$size)


finalmodel <- cbind(
  final.structure,
  transformed.variables[, c("macroalgae", "overgrowth", "Turf", "Light", "Temp", "currents", "waves", "wind", "Tides")]
)

##---------------------------------------------
## DAG-informed covariates
##---------------------------------------------
mydag <- survival_dag

temp.covars      <- c("Temp", adjustmentSets(mydag, "Temp", "survival")[[1]])
Light.covars     <- c("Light", adjustmentSets(mydag, "Light", "survival")[[1]])
tides.covars     <- c("Tides", adjustmentSets(mydag, "Tides", "survival")[[1]])
current.covars   <- c("currents", adjustmentSets(mydag, "currents", "survival")[[1]])
Turf.covars      <- c("Turf", adjustmentSets(mydag, "Turf", "survival")[[1]])
overgrowth.covars <- c("overgrowth", adjustmentSets(mydag, "overgrowth", "survival")[[2]])

mymodels <- list(
  te = temp.covars,
  li = Light.covars,
  cu = current.covars,
  tu = Turf.covars,
  ov = overgrowth.covars,
  ti = tides.covars
)

save(finalmodel, mymodels, file = "Survival.RData")

##====================================================
##   BAYESIAN MODELS — WITH 50% & 95% CREDIBLE INTERVALS
##====================================================

ncores <- detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = ncores)

model_results <- list()
posterior_summaries <- list()

for (model_name in names(mymodels)) {
  predictors <- mymodels[[model_name]]
  formula_str <- paste(
    "survival ~", paste(predictors, collapse = " + "),
    "+ size + Days + (1 |SiteQ/UniqueID)"
  )
  
  fit <- brm(
    formula = as.formula(formula_str),
    data = finalmodel,
    family = bernoulli(link = "logit"),
    prior = c(
      prior(normal(0, 2.5), class = "b"),
      prior(normal(0, 1), class = "sd"),
      prior(normal(2.28, 1.5), class = "Intercept")
    ),
    chains = 4, cores = 4, iter = 4000, warmup = 2000,
    control = list(adapt_delta = 0.999),
    silent = TRUE, refresh = 0
  )
  
  model_results[[model_name]] <- fit
  
  draws <- as_draws_df(fit)
  mypred <- names(draws)[2]
  
  coef.summary <- tibble(
    term = tools::toTitleCase(gsub("^b_", "", mypred)),
    estimate = exp(mean(draws[[mypred]])),
    lower_95 = exp(quantile(draws[[mypred]], 0.025)),
    upper_95 = exp(quantile(draws[[mypred]], 0.975)),
    lower_50 = exp(quantile(draws[[mypred]], 0.25)),
    upper_50 = exp(quantile(draws[[mypred]], 0.75)),
    prob.dir = mean(draws[[mypred]] < 0),
    model = model_name
  )
  
  posterior_summaries[[model_name]] <- coef.summary
}

all_estimates <- bind_rows(posterior_summaries)
all_estimates$term <- factor(tools::toTitleCase(all_estimates$term),
                             levels = c("Overgrowth", "Turf", "Temp", "Light", "Currents", "Tides"))

all_estimates <- all_estimates %>%
  mutate(
    effect_type = case_when(
      prob.dir > 0.975 ~ "negative",
      prob.dir < 0.025 ~ "positive",
      TRUE ~ "uncertain"
    ),
    percentchange = (estimate - 1) * 100
  )

save(rawdata, conditional.test, mymodels, finalmodel, model_results,
     posterior_summaries, all_estimates,
     file = paste0(df_name, "_Survival_models.RData"))

##---------------------------------------------
## PLOT RESULTS
##---------------------------------------------
sby <- ggplot(all_estimates, aes(y = term, x = estimate, colour = prob.dir)) +
  geom_pointrange(aes(xmin = lower_50, xmax = upper_50), linewidth = 0.8, color = "gray20") +
  geom_pointrange(aes(xmin = lower_95, xmax = upper_95), linewidth = 0.2, color = "gray40") +
  geom_point(size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_gradient2(
    low = "#0072B2", mid = "gray70", high = "#D90E00",
    midpoint = 0.5, limits = c(0, 1)
  ) +
  ggtitle("Survival") +
  labs(x = "Posterior Odds ratio", y = "") +
  theme_minimal() +
  theme(
    legend.position = "none",
    title = element_text(size = 8),
    panel.grid.major = element_line(colour = "grey98"),
    panel.grid.minor = element_line(colour = "grey98")
  ) +
  labs(tag = "a")

ggsave(
  filename = file.path(myplots, paste0(df_name, "_survival_model_estimates.png")),
  plot = sby, dpi = 300, width = 3, height = 3.5
)
