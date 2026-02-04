# ============================================================
# Script: coral_biophysical_dag.R
# Purpose: Build and visualize DAGs for biophysical drivers of coral demography
# Date: 13-09-2025
# ============================================================

# ---- Load required packages ----
library(dagitty)
library(ggdag)
library(ggplot2)
library(here)
library(tools)

# ---- Set output directory ----
plot_dir <- here("Plots")

# Create directory if it doesn’t exist
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ============================================================
# 1. Define DAGs (Created in Dagitty)
# ============================================================

# --- DAG for coral growth ---
growth_dag <- dagitty('dag {
bb="-2.928,-3.545,4.676,3.504"
CCA [pos="0.475,-0.390"]
Light [pos="0.319,-1.760"]
Temp [pos="-0.625,-0.613"]
Tides [pos="-0.442,-2.968"]
Turf [pos="1.230,-2.014"]
currents [pos="-0.636,-1.767"]
grazing [latent,pos="2.152,0.310"]
macroalgae [pos="1.041,1.642"]
overgrowth [exposure,pos="0.038,2.394"]
growth [outcome,pos="-1.198,2.671"]
waves [pos="-1.569,-0.567"]
wind [pos="-1.542,-2.121"]

CCA -> growth
Light -> CCA
Light -> Temp
Light -> Turf
Light -> macroalgae
Light -> growth
Temp -> CCA
Temp -> macroalgae
Temp -> growth
Tides -> Light
Tides -> currents
Tides -> waves
Tides -> Temp
Turf -> CCA
Turf -> overgrowth
currents -> Light
currents -> Temp
currents -> Turf
currents -> macroalgae
currents -> growth
grazing -> Turf
grazing -> macroalgae
macroalgae -> CCA
macroalgae -> Turf
macroalgae -> overgrowth
overgrowth -> growth
waves -> CCA
waves -> Light
waves -> Temp
waves -> currents
waves -> grazing
waves -> growth
wind -> Light
wind -> Temp
wind -> Turf
wind -> currents
wind -> waves
}')

## Create survival dag by replacing growth response variable

survival_dag <- as.character(growth_dag)
survival_dag <- gsub("\\bgrowth\\b", "survival", survival_dag)
survival_dag <- dagitty(survival_dag)

# ============================================================
# 2. Prepare node metadata
# ============================================================

latent_vars <- c("grazing")
env_vars    <- c("wind", "Tides", "waves", "currents", "Temp", "Light")
biotic_vars <- c("Turf", "CCA", "macroalgae", "overgrowth")
response    <- c("growth")

node_info <- data.frame(
  name = c(latent_vars, env_vars, biotic_vars, response),
  status = c(
    "latent",
    rep("environment", length(env_vars)),
    rep("biotic", length(biotic_vars)),
    "response"
  )
)

# ============================================================
# 3. Tidy and relabel for plotting
# ============================================================

mydag <- tidy_dagitty(growth_dag)
mydag <- merge(mydag, node_info, by = "name", all.x = TRUE)

# Rename key nodes
mydag$name[mydag$name == "growth"]   <- "Growth/Survival"
mydag$name[mydag$name == "grazing"]  <- "Herbivores"

# Title case all node names
mydag$name <- tools::toTitleCase(mydag$name)

# ============================================================
# 4. Plot DAG
# ============================================================

p <- ggdag(growth_dag, text_size = 2, text_col = "black") +
  theme_dag_blank() +
  geom_dag_edges(edge_colour = "grey60", edge_width = 0.5) +
  geom_dag_node(aes(color = status), data = mydag, alpha = 0.9) +
  geom_dag_text(data = mydag, color = "black", size = 2) +
  scale_color_manual(values = c(
    latent = "grey90",
    environment = "lightblue",
    biotic = "lightgreen",
    response = "pink"
  )) +
  guides(color = "none") +
  scale_y_reverse() +
  labs(title = "Biophysical Drivers of Coral Demography")

# ============================================================
# 5. Save plot
# ============================================================

# ggsave(
#   filename = file.path(plot_dir, "coral_biophysical_dag.png"),
#   plot = p,
#   dpi = 400,
#   width = 5,
#   height = 5
# )
