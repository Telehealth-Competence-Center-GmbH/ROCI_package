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
#' @param event One of "inside" or "outside" indicating whether class 1 should be
#' predicted inside or outside of the interval.
#'
#' @return A data.frame containing the cut-off values, sensitivity, specificity
#' and the Youden-Index for all possible intervals.
#' @export
#'
#' @examples
#' set.seed(123)
#' ph_values <- rnorm(20, 7, 2)
#' med_class <- ifelse(ph_values >= quantile(ph_values, 0.25) &
#'   ph_values <= quantile(ph_values, 0.75),
#' sample(c(0, 1), prob = c(0.3, 0.7)),
#' sample(c(0, 1), prob = c(0.8, 0.2))
#' )
#' roci(ph_values, med_class)
roci <- function(data, classes, event = "inside") {

  if(any(duplicated(data))){
    warning("Data contains ties. ROCI curves may not be an appropriate tool. Consider making values unique by adding some minimal noise.")
  }

  if (event == "outside") {
    # switch classes 0 and 1 and correct at the end of the function
    classes = 1 - classes
  } else if (event != "inside") {
    stop("Only event = \"inside\" or event = \"outside\" allowed.")
  }

  data_df <- data.frame(data = data, classes = classes) %>%
    summarize(events = sum(classes),
              non_events = sum(1 - classes),
              .by = data) %>%
    arrange(data)

  ord <- order(data)
  data_sorted <- data_df$data
  events_cumsum <- cumsum(data_df$events)
  non_events_cumsum <- cumsum(data_df$non_events)
  n <- length(data_sorted)

  # Initialize results
  n_rows <- sum(seq(1, n-1))
  width_ <- numeric(n_rows)
  rp_ <- numeric(n_rows)
  rn_ <- numeric(n_rows)
  fp_ <- numeric(n_rows)
  fn_ <- numeric(n_rows)
  l_cut_ <- numeric(n_rows)
  u_cut_ <- numeric(n_rows)
  row_idx <- 1
  # Laufe über Breite des Intervalls
  for (int_width in 1:(n - 1)) {
    if (int_width %% 100 == 0) print(int_width)
    s <- n - int_width
    for (int_pos in 0:s) {
      l_cut_[row_idx] <- data_sorted[int_pos + 1]
      u_cut_[row_idx] <- data_sorted[int_pos + int_width]
      # Anzahl richtig positiver/negativer und falsch positiver/negativer
      if (int_pos > 0) {
        rp_[row_idx] <- events_cumsum[int_pos + int_width] - events_cumsum[int_pos]
        fp_[row_idx] <- non_events_cumsum[int_pos + int_width] - non_events_cumsum[int_pos]
        fn_[row_idx] <- events_cumsum[int_pos] + events_cumsum[n] - events_cumsum[int_pos + int_width]
        rn_[row_idx] <- non_events_cumsum[int_pos] + non_events_cumsum[n] - non_events_cumsum[int_pos + int_width]
      } else {
        rp_[row_idx] <- events_cumsum[int_pos + int_width]
        fp_[row_idx] <- non_events_cumsum[int_pos + int_width]
        fn_[row_idx] <- events_cumsum[n] - events_cumsum[int_pos + int_width]
        rn_[row_idx] <- non_events_cumsum[n] - non_events_cumsum[int_pos + int_width]
      }
      width_[row_idx] = rp_[row_idx] + fp_[row_idx]
      row_idx <- row_idx + 1
    }
  }

  sensitivity_ <- rp_ / (rp_ + fn_)
  specificity_ <- rn_ / (rn_ + fp_)
  J_ <- sensitivity_ + specificity_ - 1

  results <- data.frame(
    width = width_,
    rp = rp_,
    rn = rn_,
    fp = fp_,
    fn = fn_,
    sensitivity = sensitivity_,
    specificity = specificity_,
    J = J_,
    l.cutoff = l_cut_,
    u.cutoff = u_cut_
  )

  results <- results[order(results$width, results$J),]

  if (event == "outside") {
    # switch sensitivity and specificity etc.
    colnames(results) <- c("width", "rn", "rp", "fn", "fp", "specificity", "sensitivity", "J", "l.cutoff", "u.cutoff")
    results <- results[c("width", "rp", "rn",  "fp", "fn", "sensitivity", "specificity", "J", "l.cutoff", "u.cutoff")]
  }
  attr(results, "event") <- event
  return(results)
}

#' Compute area under the curve (AUC) for curve provided by ROCI algorithm.
#'
#' Compute area under the curve (AUC) for curve computed via `roci()` function using the maximal sensitivity corresponding to each specificity.
#'
#' @return A numeric giving the AUC value.
#' @export
#'
#' @examples
#' set.seed(123)
#' ph_values <- rnorm(20, 7, 2)
#' med_class <- ifelse(ph_values >= quantile(ph_values, 0.25) &
#'   ph_values <= quantile(ph_values, 0.75),
#' sample(c(0, 1), prob = c(0.3, 0.7)),
#' sample(c(0, 1), prob = c(0.8, 0.2))
#' )
#' roci_results <- roci(ph_values, med_class)
#' auc_roci(roci_results)

auc_roci <- function(roci_results) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("magrittr", quietly = TRUE)) {
    stop(
      "Package \"magrittr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # Add bottom left and top right value for roci curve
  event <- attr(roci_results, "event")

  if (sum(roci_results$width == 0) == 0) {
    roci_results <- roci_results %>%
      dplyr::add_row(
        width = 0,
        sensitivity = ifelse(event == "inside", 0, 1),
        specificity = ifelse(event == "inside", 1, 0),
        J = 0
      )
  }

  if (sum(roci_results$width == max(roci_results$width) + 1) == 0) {
    roci_results <- roci_results %>%
      dplyr::add_row(
        width = max(roci_results$width) + 1,
        sensitivity = ifelse(event == "inside", 1, 0),
        specificity = ifelse(event == "inside", 0, 1),
        J = 0
      )
  }

  # Compute AUC via Riemann sum
  roci_comp_auc <- roci_results %>%
    dplyr::mutate(fpr = 1 - specificity) %>%
    dplyr::select(fpr, sensitivity) %>%
    dplyr::slice_max(sensitivity, by = fpr, with_ties = FALSE) %>%
    dplyr::arrange(fpr) %>%
    dplyr::mutate(fpr_next = dplyr::lead(fpr),
                  sensitivity_next = dplyr::lead(sensitivity),
                  df = (sensitivity + sensitivity_next)/2 * (fpr_next - fpr))

  auc <- sum(roci_comp_auc$df, na.rm = TRUE)

  return(auc)
}




#' Plot results of ROCI algorithm
#'
#' Plot a ROC curve for intervals as calculated by `roci()`.
#'
#' @param roci_results A data frame as returned by `roci()`.
#' @param hide_whisker Logical if suboptimal sensitivity-specificity pairs for
#' fixed interval lengths should be hidden in the plot.
#' @param add_width_plot Logical if width of intervals should be added to plot.
#' @param add_yi_plot Logical if maximal Youden-Index should be written to plot
#' and if the corresponding interval should be highlighted.
#' @param add_cutoff_plot Logical if corresponding interval of maximal
#' Youden-Index should be written to plot.
#' @param add_auc_plot Logical if corresponding AUC should be written to plot.
#' @param ... Additional arguments passed to ggplot.
#'
#' @return A ggplot.
#' @export
#'
#' @examples
#' set.seed(123)
#' ph_values <- rnorm(20, 7, 2)
#' med_class <- ifelse(ph_values >= quantile(ph_values, 0.25) &
#'   ph_values <= quantile(ph_values, 0.75),
#' sample(c(0, 1), prob = c(0.3, 0.7)),
#' sample(c(0, 1), prob = c(0.8, 0.2))
#' )
#' roci_results <- roci(ph_values, med_class)
#' plot_roci(roci_results)
plot_roci <- function(roci_results,
                      hide_whiskers = FALSE,
                      add_width_plot = TRUE,
                      add_yi_plot = TRUE,
                      add_cutoff_plot = TRUE,
                      add_auc_plot = TRUE,
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

  if(add_auc_plot) {
    auc <- auc_roci(roci_results)
  }

  # Add bottom left and top right value for roci curve
  event <- attr(roci_results, "event")

  if (sum(roci_results$width == 0) == 0) {
    roci_results <- roci_results %>%
      dplyr::add_row(
        width = 0,
        sensitivity = ifelse(event == "inside", 0, 1),
        specificity = ifelse(event == "inside", 1, 0),
        J = 0
      )
  }

  if (sum(roci_results$width == max(roci_results$width) + 1) == 0) {
    roci_results <- roci_results %>%
      dplyr::add_row(
        width = max(roci_results$width) + 1,
        sensitivity = ifelse(event == "inside", 1, 0),
        specificity = ifelse(event == "inside", 0, 1),
        J = 0
      )
  }

  roci_results_gg <- roci_results %>%
    dplyr::group_by(width) %>%
    dplyr::mutate(max_j = ifelse(J == max(J), 1, 0)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(width)

  roci_plot <- ggplot2::ggplot(roci_results_gg) +
    ggplot2::geom_point(ggplot2::aes(x = 1 - specificity, y = sensitivity)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = 1, y = 0, yend = 1),
      linetype = "dotted"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = 1 - specificity, y = sensitivity, group = width)
    ) +
    ggplot2::geom_path(
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
  if (add_auc_plot) {
    roci_plot <- roci_plot +
      ggplot2::annotate(
        geom = "text",
        x = 0.65, y = 0.3, label = paste0(
          "AUC = ",
          round(auc, 3),
          "\n"
        ),
        hjust = 0
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
