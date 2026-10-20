#' Extract the coefficient matrix from a fit (internal)
#'
#' Accepts a matrix, the output of \code{run_FCR_ISTA()} (a list with
#' \code{$beta}), or a cross-validation result holding the final fit in
#' \code{$sparse_fit$beta}.
#'
#' @param fit A beta matrix or a fit object.
#' @return The beta matrix (genes x transitions, possibly with an intercept row).
#' @keywords internal
#' @noRd
get_beta_matrix <- function(fit) {
  if (is.matrix(fit)) {
    return(fit)
  }
  if (is.list(fit) && is.matrix(fit$beta)) {
    return(fit$beta)
  }
  if (is.list(fit) && is.list(fit$sparse_fit) && is.matrix(fit$sparse_fit$beta)) {
    return(fit$sparse_fit$beta)
  }
  stop("Could not find a beta matrix: pass a matrix, a run_FCR_ISTA() fit, or an object with $sparse_fit$beta.")
}

#' Apply a magnitude threshold to the gene coefficients (internal)
#'
#' Sets gene coefficients with absolute value at or below \code{threshold} to
#' zero. The intercept row, if present, is never thresholded.
#'
#' @param beta Beta matrix.
#' @param threshold Non-negative magnitude threshold.
#' @return The thresholded beta matrix.
#' @keywords internal
#' @noRd
threshold_beta <- function(beta, threshold) {
  gene_rows <- rownames(beta) != "Intercept"
  if (is.null(rownames(beta))) {
    gene_rows <- rep(TRUE, nrow(beta))
  }
  sub <- beta[gene_rows, , drop = FALSE]
  sub[abs(sub) <= threshold] <- 0
  beta[gene_rows, ] <- sub
  beta
}

#' Select genes with a nonzero effect at any transition
#'
#' Thresholds the estimated coefficients and keeps the genes whose coefficient
#' exceeds \code{threshold} in absolute value at one or more transitions.
#' Works for any fitted FCR model, not only simulated data.
#'
#' @param fit A beta matrix, the output of \code{run_FCR_ISTA()}, or an object
#'   holding the final fit in \code{$sparse_fit}.
#' @param threshold Non-negative number. Coefficients with
#'   \code{abs(beta) <= threshold} are treated as zero. Default 0.
#' @param rank_by How to rank the selected genes: \code{"max"} (largest absolute
#'   coefficient across transitions, the default) or \code{"norm"} (Euclidean norm
#'   across transitions).
#'
#' @return A matrix of the selected genes (rows) by transitions (columns), with
#'   sub-threshold entries set to zero and rows ordered from strongest to
#'   weakest effect. The intercept row is dropped. The ranking score is stored
#'   in \code{attr(, "score")}.
#'
#' @examples
#' sim <- generate_sim_data(N = 500, P = 20)
#' fit <- run_FCR_ISTA(sim$X, sim$Y, lambda1 = 0.05, lambda2 = 0.05)
#' select_genes(fit, threshold = 0.1)
#'
#' @export
select_genes <- function(fit, threshold = 0, rank_by = c("max", "norm")) {
  rank_by <- match.arg(rank_by)
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
    stop("threshold must be a single non-negative number.")
  }
  beta <- get_beta_matrix(fit)
  if (is.null(rownames(beta))) {
    rownames(beta) <- paste0("Gene_", seq_len(nrow(beta)))
  }
  if (is.null(colnames(beta))) {
    colnames(beta) <- paste0("Time_", seq_len(ncol(beta)))
  }
  beta <- beta[rownames(beta) != "Intercept", , drop = FALSE]
  beta <- threshold_beta(beta, threshold)
  keep <- apply(beta != 0, 1, any)
  selected <- beta[keep, , drop = FALSE]
  if (rank_by == "max") {
    score <- apply(abs(selected), 1, max)
  } else {
    score <- sqrt(rowSums(selected^2))
  }
  ord <- order(score, decreasing = TRUE)
  selected <- selected[ord, , drop = FALSE]
  attr(selected, "score") <- score[ord]
  selected
}

#' Lollipop plots of transition-specific coefficients for the top genes
#'
#' Selects genes with a nonzero coefficient at any transition (see
#' \code{\link{select_genes}}), keeps the \code{top_n} with the largest effect,
#' and draws one lollipop panel per gene showing its coefficient at each
#' transition. Panels are laid out with at most \code{ncol} per row, so
#' \code{top_n = 10} with \code{ncol = 6} gives two rows. When \code{file} is
#' given, the saved figure grows in height with the number of rows.
#'
#' @inheritParams select_genes
#' @param top_n Integer number of genes to plot. Default 6. If fewer genes are
#'   selected, all of them are plotted.
#' @param ncol Maximum number of panels per row. Default 6.
#' @param transition_labels Optional character vector of x-axis labels, one per
#'   transition (e.g. \code{c("D0-D3", "D3-D6", "D6-D8")}). Defaults to
#'   \code{1, 2, ...}.
#' @param scales Facet y-axis scaling: \code{"fixed"} (shared axis, the default,
#'   so effect sizes are comparable across genes) or \code{"free_y"}.
#' @param colors Named vector of two colors for positive and negative effects.
#' @param base_size Base font size. Default 9.
#' @param file Optional path (e.g. \code{"lollipop.pdf"} or \code{".png"}). If
#'   supplied, the plot is saved with \code{ggplot2::ggsave()}.
#' @param width_per_col,height_per_row Panel size in inches used when saving.
#'   Total width is \code{width_per_col * (panels per row)} and total height is
#'   \code{height_per_row * (number of rows) + 0.6}.
#' @param dpi Resolution for raster formats. Default 300.
#'
#' @return The \code{ggplot} object, invisibly when \code{file} is supplied.
#'
#' @examples
#' sim <- generate_sim_data(N = 500, P = 20)
#' fit <- run_FCR_ISTA(sim$X, sim$Y, lambda1 = 0.05, lambda2 = 0.05)
#' plot_lollipop(fit, threshold = 0.1, top_n = 3)
#'
#' @importFrom ggplot2 .data
#' @export
plot_lollipop <- function(fit, threshold = 0, top_n = 6, ncol = 6, rank_by = c("max", "norm"), transition_labels = NULL, scales = c("fixed", "free_y"), colors = c(Positive = "#B2182B", Negative = "#2166AC"), base_size = 9, file = NULL, width_per_col = 1.6, height_per_row = 1.8, dpi = 300) {
  scales <- match.arg(scales)
  selected <- select_genes(fit, threshold = threshold, rank_by = rank_by)
  if (nrow(selected) == 0) {
    stop("No genes have |beta| > ", threshold, " at any transition. Try a smaller threshold.")
  }
  if (nrow(selected) < top_n) {
    message("Only ", nrow(selected), " genes pass the threshold; plotting all of them.")
    top_n <- nrow(selected)
  }
  top <- selected[seq_len(top_n), , drop = FALSE]
  TT <- ncol(top)
  if (is.null(transition_labels)) {
    transition_labels <- as.character(seq_len(TT))
  }
  if (length(transition_labels) != TT) {
    stop("transition_labels must have one label per transition (", TT, ").")
  }
  plot_df <- data.frame(gene = factor(rep(rownames(top), times = TT), levels = rownames(top)), transition = factor(rep(transition_labels, each = top_n), levels = transition_labels), beta = as.vector(top))
  plot_df$direction <- ifelse(plot_df$beta > 0, "Positive", ifelse(plot_df$beta < 0, "Negative", "Zero"))
  n_col <- min(ncol, top_n)
  n_row <- ceiling(top_n / n_col)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$transition, y = .data$beta)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$transition, y = 0, yend = .data$beta), color = "grey40", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = .data$direction), size = 2) +
    ggplot2::scale_color_manual(values = c(colors, Zero = "grey60"), breaks = c("Positive", "Negative"), name = "Effect") +
    ggplot2::facet_wrap(~ gene, ncol = n_col, scales = scales) +
    ggplot2::labs(x = "Transition", y = expression(hat(beta))) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"), panel.border = ggplot2::element_rect(color = "grey80", fill = NA, linewidth = 0.3))
  if (!is.null(file)) {
    ggplot2::ggsave(file, p, width = width_per_col * n_col, height = height_per_row * n_row + 0.6, dpi = dpi)
    return(invisible(p))
  }
  p
}

#' Expected pseudotime from a fitted FCR model
#'
#' Predicts stage probabilities for new cells, optionally after thresholding
#' small gene coefficients to zero, and returns the expected stage
#' \eqn{\sum_k k\,P(Y = k)}, rescaled to [0, 1] by default.
#'
#' @inheritParams select_genes
#' @param X_new Matrix of new observations (cells x genes), on the same scale
#'   as the data used to fit the model, without an intercept column.
#' @param y_levels Numeric values attached to the stages, used for the
#'   expectation (e.g. \code{c(0, 3, 6, 8)} for days). Defaults to
#'   \code{1, 2, ..., T}.
#' @param rescale Logical; if \code{TRUE} (default), rescale to [0, 1] using
#'   \code{min(y_levels)} and \code{max(y_levels)}.
#'
#' @return A list with \code{pseudotime} (numeric vector, one per cell),
#'   \code{prob_mat} (cells x stages probability matrix), and
#'   \code{pred_class} (most likely stage index).
#'
#' @examples
#' sim <- generate_sim_data(N = 500, P = 20)
#' fit <- run_FCR_ISTA(sim$X, sim$Y, lambda1 = 0.05, lambda2 = 0.05)
#' pt <- pseudotime_FCR(fit, sim$X, threshold = 0.1)
#' summary(pt$pseudotime)
#'
#' @export
pseudotime_FCR <- function(fit, X_new, threshold = 0, y_levels = NULL, rescale = TRUE) {
  beta <- threshold_beta(get_beta_matrix(fit), threshold)
  pred <- predict_FCR(X_new = X_new, beta_matrix = beta)
  prob_mat <- pred$probabilities
  K <- ncol(prob_mat)
  if (is.null(y_levels)) {
    y_levels <- seq_len(K)
  }
  if (length(y_levels) != K) {
    stop("y_levels must have one value per stage (", K, ").")
  }
  colnames(prob_mat) <- paste0("P(Y=", y_levels, ")")
  pseudotime <- as.vector(prob_mat %*% y_levels)
  if (rescale) {
    pseudotime <- (pseudotime - min(y_levels)) / (max(y_levels) - min(y_levels))
  }
  list(pseudotime = pseudotime, prob_mat = prob_mat, pred_class = pred$predictions)
}

#' Density plot of predicted pseudotime by true stage
#'
#' Draws overlaid density curves of predicted pseudotime, one per true stage,
#' to show how well the predicted ordering separates the stages. Takes the
#' pseudotime from \code{\link{pseudotime_FCR}} (or from any other method, for
#' comparison).
#'
#' @param pseudotime Numeric vector of predicted pseudotime, or the list
#'   returned by \code{pseudotime_FCR()}.
#' @param true_stage Vector of true stages, same length as \code{pseudotime}.
#'   Converted to a factor; the factor level order sets the legend order.
#' @param colors Optional named vector of colors, one per stage (names matching
#'   the stage labels). Defaults to the viridis palette.
#' @param bw Density bandwidth. Default \code{"nrd0"} (R's default rule); for
#'   pseudotime on [0, 1], small values such as 0.008 show more detail.
#' @param xlim Range of the x axis. Default \code{c(0, 1)}; use \code{NULL} for
#'   unscaled pseudotime.
#' @param x_lab,y_lab,legend_title Axis and legend titles.
#' @param title Optional plot title.
#' @param base_size Base font size. Default 9.
#' @param file Optional path to save the plot with \code{ggplot2::ggsave()}.
#' @param width,height Size in inches when saving. Default 3.5 x 2.8.
#' @param dpi Resolution for raster formats. Default 300.
#'
#' @return The \code{ggplot} object, invisibly when \code{file} is supplied.
#'
#' @examples
#' sim <- generate_sim_data(N = 500, P = 20)
#' fit <- run_FCR_ISTA(sim$X, sim$Y, lambda1 = 0.05, lambda2 = 0.05)
#' pt <- pseudotime_FCR(fit, sim$X, threshold = 0.1)
#' plot_pseudotime_density(pt, true_stage = sim$Y)
#'
#' @importFrom ggplot2 .data
#' @export
plot_pseudotime_density <- function(pseudotime, true_stage, colors = NULL, bw = "nrd0", xlim = c(0, 1), x_lab = "Expected pseudotime (scaled, 0–1)", y_lab = "Density", legend_title = "True stage", title = NULL, base_size = 9, file = NULL, width = 3.5, height = 2.8, dpi = 300) {
  if (is.list(pseudotime)) {
    pseudotime <- pseudotime$pseudotime
  }
  if (length(pseudotime) != length(true_stage)) {
    stop("pseudotime and true_stage must have the same length.")
  }
  true_stage <- as.factor(true_stage)
  df <- data.frame(pseudotime = pseudotime, true_stage = true_stage)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$pseudotime, fill = .data$true_stage, color = .data$true_stage)) +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.6, bw = bw)
  if (is.null(colors)) {
    p <- p + ggplot2::scale_fill_viridis_d() + ggplot2::scale_color_viridis_d()
  } else {
    p <- p + ggplot2::scale_fill_manual(values = colors) + ggplot2::scale_color_manual(values = colors)
  }
  if (!is.null(xlim)) {
    p <- p + ggplot2::scale_x_continuous(limits = xlim)
  }
  p <- p +
    ggplot2::labs(title = title, x = x_lab, y = y_lab, fill = legend_title, color = legend_title) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1), color = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = base_size - 3), legend.title = ggplot2::element_text(size = base_size - 2), legend.key.size = ggplot2::unit(0.4, "cm"))
  if (!is.null(file)) {
    ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi)
    return(invisible(p))
  }
  p
}
