ageAdjustmentPlot <- function(df, longevity){
  longevity$lty <- factor(2, levels = 2, labels = "Longevity")
  ggplot(df, aes(x = age, fill = processed)) +
    geom_histogram(alpha = 0.5, position = position_identity()) +
    geom_vline(data = longevity, aes(xintercept = longevity, linetype = lty)) +
    facet_wrap( ~ speciesCode, nrow = 2, scales = "free") +
    theme_bw() +
    scale_linetype_manual(name = NULL, values = "dashed")+
    labs(y = "Number of cohorts", x = "Age", fill = NULL)
}
