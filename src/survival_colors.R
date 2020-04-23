surv_colors = list(
  overall_survival = c("1" = "black", "0" = "white"),
  distant_recurrence = c("1" = "darkblue", "0" = "white")
)

saveRDS(object = surv_colors,
        file = here::here("data/Rds/survival_colors.Rds"))

#color_grid(surv_colors$distant_recurrence)