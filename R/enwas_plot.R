
#' Forest Plot for single model
#'
#' @param xwas_result the EnWAS result data frame
#'
#' @return plot result
#' @export
#'
#' @examples forest_plot(linear_res$enwas_res)

forest_plot <- function(xwas_result) {
  xwas_result$col <- as.numeric(rownames(xwas_result)) %% 2
  n <- nrow(xwas_result)
  xwas_result |> mutate(xmin = seq(0.5, n - 0.5, by = 1),
                        xmax = seq(1.5, n + 0.5, by = 1)) |>
    ggplot(aes(x = term,
               y = estimate,
               colour = estimate)) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.5,
      cex = 1)
    ) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = +Inf,
        fill = factor(col)
      ),
      color = 'black',
      alpha = 0.3
    ) +
    scale_fill_manual(values = c('white', 'grey78'), guide = 'none') +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Exposures") + ylab("Estimate (95% CI)") +
    theme_bw()  # use a white background
}






#' Title
#'
#' @param xwas_result_list list of the EnWAS result data frames
#'
#' @return plot multiple model results in a single image
#' @export
#'
#' @examples forest_plot_mult(list(linear = linear_res$enwas_res,
#'                                ns = ns_res$enwas_res))
forest_plot_mult <- function(xwas_result_list) {
  xwas_result <- do.call("rbind", xwas_result_list)
  xwas_result$EnWAS <-
    rep(names(xwas_result_list), each = nrow(xwas_result_list[[1]]))

  tem_df <- xwas_result_list[[1]]
  n <- nrow(tem_df)
  tem_df$col <- as.numeric(rownames(tem_df)) %% 2
  tem_df <- tem_df|> mutate(xmin = seq(0.5, n - 0.5, by = 1),
                            xmax = seq(1.5, n + 0.5, by = 1))
  xwas_result |>  ggplot(aes(
    x = term,
    y = estimate,
    colour = EnWAS
  ))  +
    geom_point(size = 2,position = position_dodge(width = 1)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width=0.5,
                  position = position_dodge(width = 1),cex=1) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_rect(data=tem_df,
              aes(
                xmin = xmin,
                xmax = xmax,
                ymin = -Inf,
                ymax = +Inf,
                fill = factor(col)
              ),
              color = 'black',
              alpha = 0.3
    ) +
    scale_fill_manual(values = c('white', 'grey78'), guide = 'none') +
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Exposures") + ylab("Estimate (95% CI)") +
    theme_bw()  # use a white background
}
