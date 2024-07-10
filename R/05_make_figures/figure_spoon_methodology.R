library(ggplot2)
library(ggformula)
library(patchwork)
library(here)

fn <- here("outputs", "results", "weights_output_list_humanBreast.rds")
output_list <- readRDS(file=fn)

s_g <- output_list[["s_g"]]
r_tilda <- output_list[["r_tilda"]]

dat <- data.frame(
  y = sqrt(s_g),
  x = r_tilda
) 

p1 <- ggplot(dat, aes(x=x, y=y)) +
  geom_point(size=0.5) +
  ylim(c(0.15, 1)) + xlim(c(-1,6)) + 
  labs(
    x = "average log2(count size)",
    y = "sqrt(std dev)"
  ) +
  theme_bw()

p2 <- ggplot(dat, aes(x=x, y=y)) +
  geom_point(size=0.5, alpha = 0.2) +
  stat_spline(nknots=10, linewidth=1.5, color="purple") +
  ylim(c(0.15, 1)) + xlim(c(-1,6)) + 
  labs(
    x = "average log2(count size)",
    y = "sqrt(std dev)"
  )  +
  annotate("text", x=3.9, y=0.91, label= "gene-wise mean-var trend") +
  theme_bw() 

p3 <- ggplot(dat, aes(x=x, y=y)) +
  stat_spline(nknots=10, linewidth=1.5, color="purple") +
  ylim(c(0.15, 1)) + xlim(c(-1,6)) + 
  labs(
    x = "fitted log2(count size)",
    y = "sqrt(std dev)"
  )  +
  annotate("text", x=3.9, y=0.91, label= "gene-wise mean-var trend") +
  theme_bw() 

ggsave(here("plots", "main", "spoon_methodology.png"),
       wrap_plots(p1, p2, p3, ncol=3, nrow=1),
       width=15, height=4)
