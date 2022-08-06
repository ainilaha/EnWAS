forest_plot <- function(xwas_result) {
    ggplot(xwas_result,aes(
      x = term,
      y = estimate,
      colour = estimate
    )) +
    geom_point(size = 2, position = position_dodge(width = 1)) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 1), width=0.5,cex=1) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Exposures") + ylab("Estimate (95% CI)") +
    theme_bw()  # use a white background
}




forest_plot_mult <- function(xwas_result_list) {
  xwas_result <- do.call("rbind", xwas_result_list)
  xwas_result$EnWAS <-
    rep(names(xwas_result_list), each = nrow(xwas_result_list[[1]]))

    ggplot(xwas_result, aes(
      x = term,
      y = estimate,
      colour = EnWAS
    ))  +
    geom_point(size = 2, position = position_dodge(width = 1)) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 1), width=0.5,cex=1) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Exposures") + ylab("Estimate (95% CI)") +
    theme_bw()  # use a white background
}
