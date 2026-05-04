library(tidyverse)
library(cowplot)

# 1. Load your data
df <- read_csv("script_05_optim_recross_0.2.csv")

# 2. Reshape the data for plotting
plot_data <- df %>%
  group_by(Start_sk) %>%
  mutate(ID = row_number()) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(rpt = row_number()) %>%
  ungroup() %>%
  mutate(ID = factor(ID)) %>%
  pivot_longer(
    cols = starts_with(c("Start", "Best")),
    names_to = c(".value", "Parameter"),
    names_pattern = "^(Start|Best)_(.*)$" ) %>%
  mutate(Parameter = factor(Parameter, levels = c("sk", "ku", "RR", "Rampl", "bpm1", "bpm2"))) %>%
  pivot_longer(c(Start, Best), names_to = "Optimization_State", values_to = "value") %>%
  mutate(Optimization_State = factor(Optimization_State, levels = c("Start", "Best")))

# a plot with facets for each parameter
# x has 2 values: initial and best
# y is the value of the parameter

limits <- list(sk = c(0, 2), ku = c(0, 10), RR = c(0, 1), 
               Rampl = c(0, 1), bpm1 = c(0, 50), bpm2 = c(50, 130))

plots <- map(names(limits), ~{
  plot_data %>%
    filter(Parameter == .x) %>%
    ggplot(aes(x = Optimization_State, y = value, group = interaction(ID, rpt), color = ID)) +
    geom_line(alpha = 0.5) +
    geom_point(color = "grey20", alpha=0.5) +
    scale_y_continuous(limits = limits[[.x]]) +
    labs(title = .x, x = "Parameter State", y = "Parameter Value") +
    theme_minimal() +
    theme(text = element_text(size = 12)) +
    theme(legend.position = "none")
})

ggsave(filename = "script_05_optimization_plots.png", plot = plot_grid(plotlist = plots, ncol = 3))