surv_colors = list(
  overall_survival = c("1" = "black", "0" = "white"),
  distant_recurrence = c("1" = "darkblue", "0" = "white"),
  involution_duration = c("<= 6 months" = "#3538F2", "> 6 months" = "#FF831E"),
  breastfeeding_duration = c("<= 1 month" = "#264653", "> 1 month" = "#E9C46A")
)


saveRDS(object = surv_colors,
        file = here::here("data/Rds/survival_colors.Rds"))

#scales::show_col(surv_colors$involution_duration)
#scales::show_col(surv_colors$breastfeeding_duration)