#' Calculates the ROC for intervals
#'
#' `roci()` calculates sensitivity, specificity and the Youden-Index for tests,
#' where two cut-offs (interval) are used for binary classification.
#' All possible intervals are tested. A full description of the algorithm is
#' available in Krannich 2023 ...
#'
#' @param data A numerical vector, which should be classified into a binary.
#' outcome.
#' @param classes A numerical binary vector with the knwon classes of the data.
#'
#' @return A data.frame containing the cut-off values, sensitivity, specificity
#' and the Youden-Index for all possible intervals.
#' @export
#'
#' @examples
#' set.seed(123)
#' ph_values <- rnorm(20, 7, 2)
#' med_class <- ifelse(ph_values >= quantile(ph_values, 0.25) &
#'   ph_values >= quantile(ph_values, 0.75),
#' sample(c(0, 1), prob = c(0.3, 0.7)),
#' sample(c(0, 1), prob = c(0.8, 0.2))
#' )
#' roci(ph_values, med_class)
roci <- function(data, classes) {
  data_sorted <- data.frame(
    data = data,
    classes = classes
  )
  data_sorted <- data_sorted[order(data_sorted$data), ]
  n <- nrow(data_sorted)

  # Initialize results
  # Number of rows sum_{i=1}^n (i)
  n_rows <- sum(seq(1, n-1))
  results <- data.frame(
    "width" = numeric(n_rows),
    "rp" = numeric(n_rows),
    "rn" = numeric(n_rows),
    "fp" = numeric(n_rows),
    "fn" = numeric(n_rows),
    "sensitivity" = numeric(n_rows),
    "specificity" = numeric(n_rows),
    "J" = numeric(n_rows),
    "l.cutoff" = numeric(n_rows),
    "u.cutoff" = numeric(n_rows)
  )
  row_idx <- 1
  # Laufe über Breite des Intervalls
  for (int_width in 1:(n - 1)) {
    s <- n - int_width
    for (int_pos in 0:s) {
      # Erstelle Vektor mit "vorhergesagten" Klassen
      pred_classes_tmp <- c(
        rep(0, int_pos),
        rep(1, int_width),
        rep(0, (n - int_pos - int_width))
      )
      # Anzahl richtig positiver/negativer und falsch positiver/negativer
      rp <- sum(ifelse(data_sorted$classes == 1 & pred_classes_tmp == 1, 1, 0))
      rn <- sum(ifelse(data_sorted$classes == 0 & pred_classes_tmp == 0, 1, 0))
      fp <- sum(ifelse(data_sorted$classes == 0 & pred_classes_tmp == 1, 1, 0))
      fn <- sum(ifelse(data_sorted$classes == 1 & pred_classes_tmp == 0, 1, 0))
      l_cut <- ifelse(
        int_width != 0,
        min(subset(data_sorted$data, pred_classes_tmp == 1)), 0
      )
      u_cut <- ifelse(
        int_width != 0,
        max(subset(data_sorted$data, pred_classes_tmp == 1)), 1
      )
      sens <- rp / (rp + fn)
      spec <- rn / (rn + fp)
      you_idx <- sens + spec - 1
      results_tmp <- data.frame(
        "width" = int_width,
        "rp" = rp,
        "rn" = rn,
        "fp" = fp,
        "fn" = fn,
        "sensitivity" = sens,
        "specificity" = spec,
        "J" = you_idx,
        "l_cutoff" = l_cut,
        "u_cutoff" = u_cut
      )
      results[row_idx, ] <- results_tmp
      row_idx <- row_idx + 1
    }
  }
  return(results)
}



#' Plot results of ROCI algorithm
#'
#' Plot a ROC curve for intervals as calculated by `roci()`.
#'
#' @param roci_results A data frame as returned by `roci()`.
#' @param add_width_plot Logical if width of intervals should be added to plot.
#' @param add_yi_plot Logical if maximal Youden-Index should be written to plot
#' and if the corresponding interval should be highlighted
#' @param add_cutoff_plot Logical if corresponding interval of maximal
#' Youden-Index should be written to plot.
#' @param ... Additional arguments passed to ggplot.
#'
#' @return A ggplot.
#' @export
#'
#' @examples
#' set.seed(123)
#' ph_values <- rnorm(20, 7, 2)
#' med_class <- ifelse(ph_values >= quantile(ph_values, 0.25) &
#'   ph_values >= quantile(ph_values, 0.75),
#' sample(c(0, 1), prob = c(0.3, 0.7)),
#' sample(c(0, 1), prob = c(0.8, 0.2))
#' )
#' roci_results <- roci(ph_values, med_class)
#' plot_roci(roci_results)
plot_roci <- function(roci_results,
                      add_width_plot = TRUE,
                      add_yi_plot = TRUE,
                      add_cutoff_plot = TRUE,
                      ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }


  if (sum(roci_results$width == 0) == 0) {
    roci_results <- roci_results %>%
      dplyr::add_row(
        width = 0,
        sensitivity = 0,
        specificity = 1,
        J = 0
      )
  }

  if (sum(roci_results$width == max(roci_results$width) + 1) == 0) {
    roci_results <- roci_results %>%
      dplyr::add_row(
        width = max(roci_results$width) + 1,
        sensitivity = 1,
        specificity = 0,
        J = 0
      )
  }

  roci_results_gg <- roci_results %>%
    dplyr::arrange(width) %>%
    dplyr::group_by(width) %>%
    dplyr::mutate(max_j = ifelse(J == max(J), 1, 0)) %>%
    dplyr::ungroup()

  roci_plot <- ggplot2::ggplot(roci_results_gg) +
    ggplot2::geom_point(ggplot2::aes(x = 1 - specificity, y = sensitivity)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = 1, y = 0, yend = 1),
      linetype = "dotted"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = 1 - specificity, y = sensitivity, group = width)
    ) +
    ggplot2::geom_line(
      data = subset(roci_results_gg, max_j == 1),
      ggplot2::aes(x = 1 - specificity, y = sensitivity)
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1)
    )

  if (add_width_plot) {
    roci_plot <- roci_plot +
      ggplot2::geom_text(
        data = dplyr::distinct(
          subset(roci_results_gg, max_j == 1),
          sensitivity, specificity, J,
          .keep_all = TRUE
        ),
        ggplot2::aes(x = 1 - specificity, y = sensitivity, label = width),
        vjust = -0.5, hjust = 1.5
      )
  }

  if (add_yi_plot) {
    roci_plot <- roci_plot +
      ggplot2::geom_point(
        data = subset(roci_results_gg, J == max(J)),
        ggplot2::aes(x = 1 - specificity, y = sensitivity), color = "red"
      ) +
      ggplot2::annotate(
        geom = "text",
        x = 0.65, y = 0.2, label = paste0(
          "J = ",
          round(max(roci_results_gg$J), 3),
          "\n"
        ),
        hjust = 0
      )
  }
  if (add_cutoff_plot) {
    roci_plot <- roci_plot +
      ggplot2::annotate(
        geom = "text",
        x = 0.65, y = 0.1,
        label = paste0(
          "lower bound = ",
          round(min(roci_results_gg[
            roci_results_gg$J ==
              max(roci_results_gg$J),
            "l.cutoff"
          ]), 3), "\n"
        ),
        hjust = 0
      ) +
      ggplot2::annotate(
        geom = "text",
        x = 0.65, y = 0,
        label = paste0(
          "upper bound = ",
          round(max(roci_results_gg[
            roci_results_gg$J ==
              max(roci_results_gg$J),
            "u.cutoff"
          ]), 3), "\n"
        ),
        hjust = 0
      )
  }
  roci_plot
}
