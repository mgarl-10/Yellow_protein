Figure S8. Number of genes with a log₂ fold change (FC) < –1 or > 1 across different treatments under low-humidity conditions compared to control conditions 


# Bar plots for genes with |log2FC|

library(ggplot2)

deg_counts <- data.frame(
  Condition = rep(c("CTL", "GFP", "YELL"), each = 2),
  Direction = rep(c("Log2FC<–1", "Log2FC>1"), times = 3),
  Count = c(17, 36, 17, 20, 47, 31)
)


# Plotting control

ctl<-deg_counts %>% filter(Condition == "CTL")

ctl_bar <- ggplot(ctl, aes(x = Condition, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c(
      "Log2FC<–1" = "#9FC6A9",
      "Log2FC>1"  = "#058F4D"
    )
  ) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(
    y = "Number of genes",
    fill = "Direction"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_blank()
  )


ggsave("FigS8_barplot_ctl.pdf", p, width = 5, height = 5)
ggsave("FigS8_barplot_ctl.svg", p, width = 5, height = 5)


# Plotting GFP

gfp<-deg_counts %>% filter(Condition == "GFP")

ctl_bar <- ggplot(gfp, aes(x = Condition, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c(
      "Log2FC<–1" = "#9FC6A9",
      "Log2FC>1"  = "#058F4D"
    )
  ) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(
    y = "Number of genes",
    fill = "Direction"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_blank()
  )


ggsave("FigS8_barplot_gfp.pdf", p, width = 5, height = 5)
ggsave("FigS8_barplot_gfp.svg", p, width = 5, height = 5)



# Plotting Yellow


yell<-deg_counts %>% filter(Condition == "YELL")

ctl_bar <- ggplot(yell, aes(x = Condition, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c(
      "Log2FC<–1" = "#9FC6A9",
      "Log2FC>1"  = "#058F4D"
    )
  ) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(
    y = "Number of genes",
    fill = "Direction"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_blank()
  )


ggsave("FigS8_barplot_yell.pdf", p, width = 5, height = 5)
ggsave("FigS8_barplot_yell.svg", p, width = 5, height = 5)


# Fisher exact test 


table_combined <- matrix(c(47, 31, 34, 56), nrow = 2,
                         dimnames = list(Direction = c("Down", "Up"),
                                         Group = c("YELL", "REF")))
fisher.test(table_combined)


