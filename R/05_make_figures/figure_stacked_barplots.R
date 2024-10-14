library(ggplot2)
library(dplyr)
library(here)

scenarios <- c(
  "Low Mean, Ovarian",
  "Low Mean, Lobular Breast",
  "Low Mean, Ductal Breast",
  "Low Mean, Subtype Breast",
  "Small Lengthscale, Ovarian",
  "Small Lengthscale, Lobular Breast",
  "Small Lengthscale, Ductal Breast",
  "Small Lengthscale, Subtype Breast",
  "High Mean, Ovarian"
)

scenarios_row1 <- list(
  c(A = 5, "A&B" = 2),  # 1: 7 genes in List 1, 2 overlaps
  c(A = 20, "A&B" = 11), # 2: 31 genes in List 1, 11 overlaps
  c(A = 18, "A&B" = 8),   # 3: 26 genes in List 1, 8 overlaps
  c(A = 1, "A&B" = 0)      # 4: 1 gene in List 1, 0 overlaps
)

scenarios_row2 <- list(
  c(A = 69, "A&B" = 8),   # 1: 77 genes in List 1, 8 overlaps
  c(A = 65, "A&B" = 19),  # 2: 84 genes in List 1, 19 overlaps
  c(A = 59, "A&B" = 19),  # 3: 78 genes in List 1, 19 overlaps
  c(A = 43, "A&B" = 16)   # 4: 59 genes in List 1, 16 overlaps
)

create_df <- function(scenarios_data, scenario_names) {
  data <- do.call(rbind, lapply(seq_along(scenarios_data), function(i) {
    scenario <- scenarios_data[[i]]
    data.frame(
      Scenario = scenario_names[i],
      Component = names(scenario),
      Count = as.numeric(scenario)
    )
  }))
  return(data)
}

df_row1 <- create_df(scenarios_row1, scenarios[1:4])
df_row2 <- create_df(scenarios_row2, scenarios[5:9])

df_combined <- bind_rows(df_row1, df_row2)

p <- ggplot(df_combined, aes(x = Scenario, y = Count, fill = Component)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c("A" = "#1f77b4", "A&B" = "#ff7f0e"),
    labels = c("A" = "All Genes", "A&B" = "Cancer Genes")
  ) +
  labs(
    title = "",
    x = "Dataset",
    y = "Count",
    fill = "Component"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_flip()

ggsave(here("plots", "main", "cancer_barplots.png"),
       plot = p,
       width = 15, height = 5)
