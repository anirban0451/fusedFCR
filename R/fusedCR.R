#' Fit FusedCR with Cross-Validated Lambda Selection
#'
#' Runs a cross-validated grid search over fused-lasso tuning parameters for the
#' forward continuation ratio (FCR) model, then refits the final model(s) on the
#' full data at the selected lambda(s).
#'
#' @param X Numeric feature matrix (samples x genes), no intercept column.
#' @param Y Numeric or integer vector of ordinal outcomes.
#' @param lam1_grid Numeric vector of candidate lambda1 (sparsity penalty) values.
#' @param lam2_grid Numeric vector of candidate lambda2 (fusion penalty) values.
#' @param cv_fold Integer number of cross-validation folds. Default 5.
#' @param cr_method Character; continuation ratio method. Currently only "FCR" is supported.
#' @param max_iter Integer maximum ISTA iterations. Default 10000.
#' @param tol Numeric convergence tolerance for ISTA. Default 1e-4.
#' @param T_max Integer total number of ordinal stages (TT = T_max - 1 transitions modeled). Default 6.
#' @param cores Integer number of parallel workers for the CV grid search. Default 10.
#' @param fit Character; which final model(s) to fit and return. One of
#'   \code{"min"} (minimum CV log-loss), \code{"sparse"} (1SE sparsity-selected model),
#'   or \code{"both"}. Default \code{"both"}.
#'
#' @return A list with:
#' \describe{
#'   \item{final_fit}{If \code{fit = "min"} or \code{"sparse"}, the fitted model list
#'     from \code{run_CR_ISTA}. If \code{fit = "both"}, a named list with elements
#'     \code{$min} and \code{$sparse}, each a fitted model list.}
#'   \item{final_fit_lambdas}{If \code{fit = "min"} or \code{"sparse"}, a one-row
#'     data frame of the selected lambda1/lambda2 (and CV metrics). If \code{fit = "both"},
#'     a named list with elements \code{$min} and \code{$sparse}.}
#' }
#'
#' @export
FusedCR = function(X, Y, lam1_grid, lam2_grid,
                    cv_fold = 5, cr_method = "FCR",
                    max_iter = 10000, tol = 1e-4, T_max = 6, cores = 10,
                    fit = "both"){
  if(cr_method == "FCR"){
    run_CR_ISTA = run_FCR_ISTA
    predict_CR = predict_FCR
  }
  else{
    stop("cr_method should be populated with 'FCR'." ) #either 'FCR' or 'BCR'."
  }
  # else if(cr_method == "BCR"){
  #   run_CR_ISTA = run_BCR_ISTA
  #   predict_CR = predict_BCR
  # }

  if(!fit %in% c("min", "sparse", "both")){
    stop("fit should be one of 'min', 'sparse', or 'both'.")
  }

  TT = T_max - 1
  cat("Executing Micro-Grid Search...\n")
  cv_output <- cv_ISTA_parallel(X=X, Y=Y, lambda1_grid = lam1_grid, lambda2_grid = lam2_grid, K_folds = cv_fold,
                                max_iter = max_iter, tol = tol, cr_method = cr_method, ncores = cores)
  # 1. Extract the Base Results
  cv_df <- cv_output$cv_results
  # 2. Extract the Target Model(s) lambdas
  min_pt <- NULL
  sparse_pt <- NULL
  if(fit %in% c("min", "both")){
    min_pt <- cv_df[which.min(cv_df$log_loss), ]
  }
  if(fit %in% c("sparse", "both")){
    sparse_pt <- extract_optimal_lambdas(cv_df, criterion = "sparsity")
  }

  # 3. Fit the Final Model(s)
  if(fit == "min"){
    final_fit <- run_CR_ISTA(X=X, Y=Y,
                             lambda1=min_pt$lambda1, lambda2=min_pt$lambda2,
                             standardize = TRUE,
                             max_iter = max_iter,
                             tol = tol)
    final_fit_lambdas <- min_pt
  }
  else if(fit == "sparse"){
    final_fit <- run_CR_ISTA(X=X, Y=Y,
                             lambda1=sparse_pt$lambda1, lambda2=sparse_pt$lambda2,
                             standardize = TRUE,
                             max_iter = max_iter,
                             tol = tol)
    final_fit_lambdas <- sparse_pt
  }
  else if(fit == "both"){
    final_fit_min_cv <- run_CR_ISTA(X=X, Y=Y,
                                    lambda1=min_pt$lambda1, lambda2=min_pt$lambda2,
                                    standardize = TRUE,
                                    max_iter = max_iter,
                                    tol = tol)
    final_fit_sp_cv <- run_CR_ISTA(X=X, Y=Y,
                                   lambda1=sparse_pt$lambda1, lambda2=sparse_pt$lambda2,
                                   standardize = TRUE,
                                   max_iter = max_iter,
                                   tol = tol)
    final_fit <- list(min = final_fit_min_cv, sparse = final_fit_sp_cv)
    final_fit_lambdas <- list(min = min_pt, sparse = sparse_pt)
  }

  return(list(final_fit = final_fit, final_fit_lambdas = final_fit_lambdas))
}

#' Parallel Cross-Validation for FCR ISTA Models
#'
#' Performs stratified K-fold cross-validation in parallel over a grid of
#' lambda1/lambda2 values for the forward continuation ratio (FCR) ISTA model,
#' returning mean log-loss, misclassification rate, and log-loss standard
#' error for each lambda pair.
#'
#' @param X Numeric feature matrix (samples x genes), no intercept column.
#' @param Y Numeric or integer vector of ordinal outcomes.
#' @param cr_method Character; continuation ratio method. Only \code{"FCR"} is
#'   currently supported (BCR support is disabled).
#' @param lambda1_grid Numeric vector of candidate lambda1 (sparsity penalty) values.
#' @param lambda2_grid Numeric vector of candidate lambda2 (fusion penalty) values.
#' @param K_folds Integer number of cross-validation folds. Default 5.
#' @param standardize Logical; whether to standardize predictors within \code{run_CR_ISTA}. Default TRUE.
#' @param max_iter Integer maximum ISTA iterations. Default 1000.
#' @param tol Numeric convergence tolerance for ISTA. Default 1e-4.
#' @param ncores Integer number of parallel workers to request (capped at
#'   \code{parallel::detectCores() - 2}). Default 10.
#' @param seed Integer random seed for fold assignment. Default 123.
#' @return A list with a single element \code{cv_results}, a data frame with
#'   columns \code{lambda1}, \code{lambda2}, \code{log_loss}, \code{misclass_rate},
#'   and \code{log_loss_se}, one row per lambda1/lambda2 combination.
#'
#' @importFrom foreach %dopar% foreach
#' @export
cv_ISTA_parallel <- function(X, Y, cr_method = c("FCR"), lambda1_grid, lambda2_grid, K_folds = 5
                             , standardize = TRUE, max_iter = 1000, tol = 1e-4,
                             ncores = 10,
                             seed = 123
) {

  N <- nrow(X)
  cr_method <- match.arg(cr_method)

  if(cr_method == "FCR"){
    run_CR_ISTA = run_FCR_ISTA
    predict_CR = predict_FCR
  }
  # else if(cr_method == "BCR"){
  #   run_CR_ISTA = run_BCR_ISTA
  #   predict_CR = predict_BCR
  # }
  else{
    stop("cr_method should be populated with 'FCR'.") # either 'FCR' or 'BCR'."
  }

  # Setup stratified data folds
  set.seed(seed)
  Y_fac <- as.factor(Y)

  class_counts <- table(Y_fac)
  if (any(class_counts < K_folds)) {
    warning("Some outcome classes have fewer observations than K_folds.")
  }

  fold_indices <- integer(length(Y_fac))

  for (cls in levels(Y_fac)) {
    cls_idx <- which(Y_fac == cls)
    cls_idx <- sample(cls_idx, length(cls_idx), replace = FALSE)
    fold_indices[cls_idx] <- rep(1:K_folds, length.out = length(cls_idx))
  }

  # Flatten the Task Grid
  task_grid <- expand.grid(lambda1 = lambda1_grid,
                           lambda2 = lambda2_grid,
                           fold = 1:K_folds)

  cat(sprintf("Starting parallel CV for %d total models...\n", nrow(task_grid)))


  n_workers <- max(1, min(ncores, parallel::detectCores() - 2))
  cl <- parallel::makeCluster(n_workers)

  doParallel::registerDoParallel(cl)

  parallel::clusterExport(
    cl,
    varlist = c(
      "task_grid", "X", "Y", "fold_indices",
      "run_CR_ISTA", "predict_CR",
      "standardize", "max_iter", "tol", "K_folds", "post_process_beta"
    ),
    envir = environment()
  )

  parallel::clusterEvalQ(cl, {
    library(fusedFCR)
    NULL
  })

  on.exit({
    parallel::stopCluster(cl)
  }, add = TRUE)


  #  Execute the Parallel Loop (%dopar%)
  parallel_results <- foreach(i = 1:nrow(task_grid), .combine = rbind) %dopar% {

    lam1 <- task_grid$lambda1[i]
    lam2 <- task_grid$lambda2[i]
    k    <- task_grid$fold[i]

    # Split the data for this fold
    test_idx  <- which(fold_indices == k)
    X_train   <- X[-test_idx, , drop = FALSE]
    Y_train   <- Y[-test_idx]
    X_test    <- X[test_idx, , drop = FALSE]
    Y_test    <- as.numeric(Y[test_idx])

    # Fit the model using our C++ engine
    fit <- run_CR_ISTA(X = X_train, Y = Y_train,
                       lambda1 = lam1, lambda2 = lam2,
                       standardize = standardize,
                       max_iter = max_iter, tol = tol)

    # Predict on the test set
    preds <- predict_CR(X_new = X_test, beta_matrix = fit$beta)

    # Calculate Log-Loss
    predicted_probs <- preds$probabilities
    true_class_probs <- predicted_probs[cbind(1:length(Y_test), Y_test)]
    true_class_probs <- pmax(pmin(true_class_probs, 1 - 1e-15), 1e-15)
    fold_log_loss <- -mean(log(true_class_probs))

    # Calculate Misclassification Rate
    fold_misclass <- mean(preds$predictions != Y_test)

    # Return metrics
    data.frame(lambda1 = lam1, lambda2 = lam2, fold = k,
               log_loss = fold_log_loss, misclass_rate = fold_misclass)
  }

  # Aggregate the results across the K folds

  # A. Calculate the Mean Log-Loss and Misclassification Rate
  avg_results <- aggregate(cbind(log_loss, misclass_rate) ~ lambda1 + lambda2,
                           data = parallel_results,
                           FUN = mean)

  # B. Calculate the Standard Error of the Log-Loss
  # SE = standard deviation / sqrt(K)
  se_results <- aggregate(log_loss ~ lambda1 + lambda2,
                          data = parallel_results,
                          FUN = function(x) sd(x) / sqrt(K_folds))

  # Rename the column so it doesn't conflict during the merge
  names(se_results)[names(se_results) == "log_loss"] <- "log_loss_se"

  # C. Merge the Means and Standard Errors together
  final_results <- merge(avg_results, se_results, by = c("lambda1", "lambda2"))

  cat("Cross-Validation Complete.\n")

  return(list(
    cv_results = final_results
  ))
}

#' Extract Optimal Lambdas Using the 1SE Rule
#'
#' Given a cross-validation results data frame, identifies the minimum
#' log-loss model, then selects a more parsimonious model within one
#' standard error of that minimum, using either a sparsity-first or
#' balanced-penalty criterion.
#'
#' @param cv_df Data frame of CV results with columns \code{lambda1}, \code{lambda2},
#'   \code{log_loss}, and \code{log_loss_se} (e.g. as returned by
#'   \code{cv_ISTA_parallel}'s \code{$cv_results}).
#' @param criterion Character; selection rule applied among models within 1SE
#'   of the minimum. One of \code{"sparsity"} (prioritize the highest lambda1,
#'   breaking ties with the highest lambda2) or \code{"combination"} (maximize
#'   the min-max normalized sum of lambda1 and lambda2). Default \code{"sparsity"}.
#'
#' @return A one-row data frame giving the selected lambda1/lambda2 (and
#'   associated CV metrics, plus \code{total_penalty} if \code{criterion = "combination"}).
#'
#' @export
extract_optimal_lambdas <- function(cv_df, criterion = "sparsity") {

  # 1. Identify the absolute minimum model
  idx_min <- which.min(cv_df$log_loss)
  min_point <- cv_df[idx_min, ]

  # 2. Calculate the 1SE Threshold
  threshold_1se <- min_point$log_loss + min_point$log_loss_se

  # 3. Filter for all models within the 1SE boundary
  candidates <- cv_df[cv_df$log_loss <= threshold_1se, ]

  # If the threshold only leaves 1 model (the minimum itself), return it safely
  if (nrow(candidates) == 1) {
    cat("Only one model within 1SE. Returning absolute minimum.\n")
    return(candidates[1, ])
  }

  # 4. Apply the chosen selection criterion
  if (criterion == "combination") {

    # COMBINATION STYLE:
    # We use Min-Max normalization (0 to 1) instead of standard scale().
    # This prevents mathematical NaNs if all candidates share the same L1 or L2 value,
    # while ensuring a grid range of 0-6 doesn't mathematically crush a range of 0-0.5.
    range_l1 <- max(candidates$lambda1) - min(candidates$lambda1) + 1e-8
    range_l2 <- max(candidates$lambda2) - min(candidates$lambda2) + 1e-8

    norm_l1 <- (candidates$lambda1 - min(candidates$lambda1)) / range_l1
    norm_l2 <- (candidates$lambda2 - min(candidates$lambda2)) / range_l2

    candidates$total_penalty <- norm_l1 + norm_l2

    # Select the model that maximizes the combined normalized penalty
    optimal_point <- candidates[which.max(candidates$total_penalty), ]

    cat("1SE Selection: Combination (Balanced Penalties)\n")

  } else if (criterion == "sparsity") {

    # SPARSITY STYLE:
    # We order the dataframe descending by L1 first, then descending by L2.
    # This forces the algorithm to pick the model that drops the maximum number
    # of genes possible. If multiple models share that max L1, it picks the
    # one with the most temporal smoothing (highest L2) among those ties.
    sorted_candidates <- candidates[order(-candidates$lambda1, -candidates$lambda2), ]

    optimal_point <- sorted_candidates[1, ]

    cat("1SE Selection: Sparsity (Prioritizing High L1)\n")

  } else {
    stop("Invalid criterion. Choose 'combination' or 'sparsity'.")
  }

  # Print the final shift to the console
  cat(sprintf("Min Loss Point -> L1: %.4f | L2: %.4f | Loss: %.4f\n",
              min_point$lambda1, min_point$lambda2, min_point$log_loss))
  cat(sprintf("Selected 1SE   -> L1: %.4f | L2: %.4f | Loss: %.4f\n\n",
              optimal_point$lambda1, optimal_point$lambda2, optimal_point$log_loss))

  return(optimal_point)
}

#' Post-process Beta Coefficient Matrix
#'
#' Cleans a fitted beta coefficient matrix by hard-thresholding small values to
#' zero (removing floating-point noise) and rounding remaining values to a fixed
#' number of decimals so that fused/equal blocks snap together exactly.
#'
#' @param beta_matrix Numeric matrix of fitted coefficients (genes x time transitions).
#' @param zero_threshold Numeric; absolute values below this are set to 0. Default 1e-2.
#' @param fusion_decimals Integer number of decimal places to round to. Default 3.
#'
#' @return A cleaned numeric matrix of the same dimensions as \code{beta_matrix}.
#'
#' @export
post_process_beta <- function(beta_matrix, zero_threshold = 1e-2, fusion_decimals = 3) {
  clean_beta <- beta_matrix

  # 1. Hard thresholding for sparsity (crush floating-point noise to 0)
  clean_beta[abs(clean_beta) < zero_threshold] <- 0.0

  # 2. Precision rounding for fusion (snap roughly equal blocks together)
  clean_beta <- round(clean_beta, digits = fusion_decimals)

  return(clean_beta)
}

#' Calculate Log-Likelihood for FCR
#'
#' Computes the (positive) log-likelihood of a fitted forward continuation
#' ratio (FCR) model via the compiled \code{obj_g} C routine.
#'
#' @param X Numeric feature matrix (samples x genes), no intercept column.
#' @param Y Numeric or integer vector of ordinal outcomes/survival times.
#' @param beta_matrix Numeric matrix of fitted coefficients, including intercept row
#'   (P rows x TT columns), as returned by \code{run_FCR_ISTA}.
#' @param standardize Logical; whether the model was fit with standardized predictors. Default FALSE.
#'
#' @return Numeric scalar log-likelihood value.
#'
#' @useDynLib fusedFCR, .registration = TRUE
#' @export
calculate_loglik_FCR <- function(X, Y, beta_matrix, standardize = FALSE) {
  N <- as.integer(nrow(X))
  TT <- as.integer(ncol(beta_matrix))

  # THE FIX: Add the intercept column (column of 1s) so it matches beta_matrix
  X_with_int <- cbind(1, X)
  P_with_int <- as.integer(ncol(X_with_int))

  Nt <- as.integer(sapply(1:TT, function(t) sum(Y >= t)))
  std_flag <- as.integer(ifelse(standardize, 1, 0))

  result <- .C("obj_g", PACKAGE = "fusedFCR",
               YY = as.double(Y),
               XX = as.double(X_with_int), # Pass the padded matrix
               beta = as.double(beta_matrix),
               pNN = N,
               pPP = P_with_int,           # Pass the padded P
               pTT = TT,
               Nt = Nt,
               pStandardize = std_flag,
               pGG = double(1))

  return(-result$pGG)
}

#' Fit FCR Model via ISTA
#'
#' Fits a forward continuation ratio (FCR) model with fused-lasso and lasso
#' penalties via ISTA (Iterative Shrinkage-Thresholding Algorithm), using the
#' compiled \code{GGD} C routine as the optimization engine.
#'
#' @param X Numeric feature matrix (samples x genes), no intercept column.
#' @param Y Numeric or integer vector of ordinal outcomes.
#' @param lambda1 Numeric sparsity (lasso) penalty.
#' @param lambda2 Numeric fusion (fused-lasso) penalty.
#' @param standardize Logical; whether to standardize predictors before fitting. Default TRUE.
#' @param max_iter Integer maximum ISTA iterations. Default 1000.
#' @param tol Numeric convergence tolerance. Default 1e-5.
#' @param zero_threshold Numeric; absolute beta values below this are hard-thresholded to zero
#'   in post-processing. Default 1e-2.
#' @param fusion_decimals Integer number of decimals to round beta to in post-processing,
#'   so fused blocks snap together. Default 4.
#'
#' @return A list with:
#' \describe{
#'   \item{beta}{Fitted, post-processed coefficient matrix (P rows including intercept x TT columns).}
#'   \item{iterations}{Integer number of ISTA iterations actually run.}
#'   \item{objective}{Numeric vector of the objective function trace.}
#'   \item{tau_trace}{Numeric vector of the step-size trace.}
#'   \item{convergence}{Numeric vector of the convergence criterion trace.}
#' }
#'
#' @export
run_FCR_ISTA <- function(X, Y, lambda1, lambda2,
                         standardize = TRUE,
                         max_iter = 1000,
                         tol = 1e-5,
                         zero_threshold = 1e-2,
                         fusion_decimals = 4) {

  # --- DATA PREP ---
  N <- nrow(X)

  # Ensure X has an intercept! (Assuming X doesn't already have one)
  # If X already has an intercept, skip this cbind.
  X_mat <- cbind(1, as.matrix(X))
  P <- ncol(X_mat) # P now includes the intercept

  # Ensure Y is numeric
  Y_num <- as.numeric(Y)

  # In discrete survival, if max(Y) is T_max, we model transitions up to T_max - 1
  # (Because everyone left at the final time point has nowhere else to transition to)
  TT <- max(Y_num) - 1

  # Compute the risk set size (N_t) for each time point t
  Nt <- sapply(1:TT, function(t) sum(Y_num >= t))

  # Initialize the parameter matrix (Beta) with zeros
  # Dimensions: P rows (genes + intercept) by TT columns (time points)
  init_beta <- matrix(0, nrow = P, ncol = TT)

  # --- EXECUTE C++ CODE ---
  # Every argument must be strictly coerced for the .C interface
  fit <- .C("GGD", PACKAGE = "fusedFCR",
            YY          = as.double(Y_num),
            XX          = as.double(X_mat),
            beta        = as.double(init_beta),
            lambda1     = as.double(lambda1),
            lambda2     = as.double(lambda2),
            pNN         = as.integer(N),
            pPP         = as.integer(P),
            pTT         = as.integer(TT),
            pNt         = as.integer(Nt),
            pStandardize= as.integer(ifelse(standardize, 1, 0)),
            pNiter      = as.integer(max_iter),
            pStop       = as.integer(1), # 1 = stop based on objective function
            pEps        = as.double(tol),
            pTau        = as.double(1.0), # Initial step size for line search
            pFactor     = as.double(0.8), # Backtracking shrinkage factor
            new_Beta    = as.double(init_beta),
            objectives  = as.double(numeric(max_iter)),
            taus        = as.double(numeric(max_iter)),
            criterions  = as.double(numeric(max_iter)),
            actualNiter = as.integer(0)
  )

  # --- FORMAT OUTPUT ---
  actual_iters <- fit$actualNiter

  # Reconstruct the flat C++ array back into an R matrix
  final_beta_matrix <- matrix(fit$new_Beta, nrow = P, ncol = TT, byrow = FALSE)
  final_beta_matrix <- post_process_beta(final_beta_matrix, zero_threshold = zero_threshold, fusion_decimals = fusion_decimals)
  rownames(final_beta_matrix) <- c("Intercept", colnames(X))
  colnames(final_beta_matrix) <- paste0("Time_", 1:TT)

  # Return a clean list of results
  return(list(
    beta       = final_beta_matrix,
    iterations = actual_iters,
    objective  = fit$objectives[1:actual_iters],
    tau_trace  = fit$taus[1:actual_iters],
    convergence= fit$criterions[1:actual_iters]
  ))
}

# #' Fit BCR Model via ISTA
# #'
# #' Fits a backward continuation ratio (BCR) model by reversing the ordinal
# #' labels and reusing the FCR ISTA C++ engine, then reversing the resulting
# #' beta columns back into chronological order.
# #'
# #' @param X Numeric feature matrix (samples x genes), no intercept column.
# #' @param Y Numeric or integer vector of ordinal outcomes.
# #' @param lambda1 Numeric sparsity (lasso) penalty.
# #' @param lambda2 Numeric fusion (fused-lasso) penalty.
# #' @param standardize Logical; whether to standardize predictors before fitting. Default TRUE.
# #' @param max_iter Integer maximum ISTA iterations. Default 1000.
# #' @param tol Numeric convergence tolerance. Default 1e-5.
# #' @param zero_threshold Numeric; absolute beta values below this are hard-thresholded to zero
# #'   in post-processing. Default 1e-2.
# #' @param fusion_decimals Integer number of decimals to round beta to in post-processing. Default 4.
# #'
# #' @return A list with \code{beta}, \code{iterations}, \code{objective},
# #'   \code{tau_trace}, and \code{convergence}, analogous to \code{run_FCR_ISTA}.
# #'
# #' @export
# run_BCR_ISTA <- function(X, Y, lambda1, lambda2,
#                          standardize = TRUE,
#                          max_iter = 1000,
#                          tol = 1e-5,
#                          zero_threshold = 1e-2,
#                          fusion_decimals = 4) {
#
#   # --- 1. DATA PREP ---
#   N <- nrow(X)
#   X_mat <- cbind(1, as.matrix(X))
#   P <- ncol(X_mat)
#
#   Y_num <- as.numeric(Y)
#   T_max <- max(Y_num)
#   TT <- T_max - 1
#
#   # --- 2. THE REVERSAL HACK ---
#   # Reverse the Y labels.
#   # If T_max = 5: A true Y of 5 becomes Y_rev = 1. A true Y of 2 becomes Y_rev = 4.
#   Y_rev <- T_max - Y_num + 1
#
#   # Compute risk sets using the reversed labels, exactly as the FCR C++ expects
#   Nt <- sapply(1:TT, function(t) sum(Y_rev >= t))
#
#   init_beta <- matrix(0, nrow = P, ncol = TT)
#
#   # --- 3. CALL UNMODIFIED C++ FCR BACKEND ---
#   # Notice we pass Y_rev instead of Y_num
#   fit <- .C("GGD",
#             YY          = as.double(Y_rev),
#             XX          = as.double(X_mat),
#             beta        = as.double(init_beta),
#             lambda1     = as.double(lambda1),
#             lambda2     = as.double(lambda2),
#             pNN         = as.integer(N),
#             pPP         = as.integer(P),
#             pTT         = as.integer(TT),
#             pNt         = as.integer(Nt),
#             pStandardize= as.integer(ifelse(standardize, 1, 0)),
#             pNiter      = as.integer(max_iter),
#             pStop       = as.integer(1),
#             pEps        = as.double(tol),
#             pTau        = as.double(1.0),
#             pFactor     = as.double(0.8),
#             new_Beta    = as.double(init_beta),
#             objectives  = as.double(numeric(max_iter)),
#             taus        = as.double(numeric(max_iter)),
#             criterions  = as.double(numeric(max_iter)),
#             actualNiter = as.integer(0)
#   )
#
#   actual_iters <- fit$actualNiter
#
#   # --- 4. FORMAT OUTPUT & REVERSE BETA ---
#   # The C++ outputs the matrix where Column 1 = T_max, and Column TT = Time 2
#   final_beta_matrix <- matrix(fit$new_Beta, nrow = P, ncol = TT, byrow = FALSE)
#   final_beta_matrix <- post_process_beta(final_beta_matrix, zero_threshold = zero_threshold, fusion_decimals = fusion_decimals)
#
#   # Reverse the columns so they represent Time 2 -> Time T_max chronologically
#   final_beta_matrix <- final_beta_matrix[, TT:1, drop = FALSE]
#
#   rownames(final_beta_matrix) <- c("Intercept", colnames(X))
#   colnames(final_beta_matrix) <- paste0("Time_", 2:T_max)
#
#   return(list(
#     beta       = final_beta_matrix,
#     iterations = actual_iters,
#     objective  = fit$objectives[1:actual_iters],
#     tau_trace  = fit$taus[1:actual_iters],
#     convergence= fit$criterions[1:actual_iters]
#   ))
# }

#' Predict Probabilities for FCR
#'
#' Computes per-class predicted probabilities and hard classifications from
#' a fitted forward continuation ratio (FCR) model.
#'
#' @param X_new Numeric feature matrix of new observations, with or without an
#'   intercept column (an intercept column is added automatically if missing).
#' @param beta_matrix Numeric matrix of fitted coefficients (P rows including
#'   intercept x TT columns), as returned by \code{run_FCR_ISTA}.
#'
#' @return A list with:
#' \describe{
#'   \item{probabilities}{Numeric matrix (N_new x (TT + 1)) of predicted class probabilities.}
#'   \item{predictions}{Integer vector of hard class predictions (argmax over \code{probabilities}).}
#' }
#'
#' @export
predict_FCR <- function(X_new, beta_matrix) {

  # 1. Ensure X_new has an intercept column
  # (Assuming beta_matrix has P rows, and X_new currently has P-1 columns)
  if (ncol(X_new) == nrow(beta_matrix) - 1) {
    X_new <- cbind(1, as.matrix(X_new))
  } else {
    X_new <- as.matrix(X_new)
  }

  N_new <- nrow(X_new)
  T_cols <- ncol(beta_matrix)

  # 2. Linear Predictor (N_new x T matrix)
  eta <- X_new %*% beta_matrix

  # 3. Discrete Hazards (N_new x T matrix)
  h <- 1 / (1 + exp(-eta))

  # 4. Unconditional Probabilities (N_new x (T + 1) matrix)
  P_mat <- matrix(0, nrow = N_new, ncol = T_cols + 1)

  # 'surv' tracks the cumulative probability of NOT transitioning
  # It starts at 1.0 (100%) for everyone at t = 1
  surv <- rep(1.0, N_new)

  for (t in 1:T_cols) {
    # P(Y = t) = P(Surviving to t) * P(Transitioning at t)
    P_mat[, t] <- surv * h[, t]

    # Update survival probability for the next time step
    surv <- surv * (1 - h[, t])
  }

  # P(Y = T + 1) is the probability of surviving past all modeled time points
  P_mat[, T_cols + 1] <- surv

  # 5. Classify based on highest probability
  # max.col efficiently finds the column index with the highest value for each row
  predicted_labels <- max.col(P_mat, ties.method = "first")

  # Return both the raw probabilities and the hard classifications
  return(list(
    probabilities = P_mat,
    predictions   = predicted_labels
  ))
}

# #' Predict Probabilities for BCR
# #'
# #' Computes per-class predicted probabilities and hard classifications from
# #' a fitted backward continuation ratio (BCR) model, working top-down from
# #' the highest ordinal stage.
# #'
# #' @param X_new Numeric feature matrix of new observations, with or without an
# #'   intercept column (an intercept column is added automatically if missing).
# #' @param beta_matrix Numeric matrix of fitted BCR coefficients (P rows including
# #'   intercept x T_cols columns).
# #'
# #' @return A list with \code{probabilities} (N_new x T_max matrix) and
# #'   \code{predictions} (integer vector of hard class predictions).
# #'
# #' @export
# predict_BCR <- function(X_new, beta_matrix) {
#
#   # 1. Ensure X_new has an intercept column
#   if (ncol(X_new) == nrow(beta_matrix) - 1) {
#     X_new <- cbind(1, as.matrix(X_new))
#   } else {
#     X_new <- as.matrix(X_new)
#   }
#
#   N_new <- nrow(X_new)
#   T_cols <- ncol(beta_matrix)
#
#   # Total number of stages is T_cols + 1
#   T_max <- T_cols + 1
#
#   # 2. Linear Predictor (N_new x T_cols matrix)
#   eta <- X_new %*% beta_matrix
#
#   # 3. Discrete Hazards (N_new x T_cols matrix)
#   # h[, j] represents P(Y = j + 1 | Y <= j + 1)
#   h <- 1 / (1 + exp(-eta))
#
#   # 4. Unconditional Probabilities (N_new x T_max matrix)
#   P_mat <- matrix(0, nrow = N_new, ncol = T_max)
#   colnames(P_mat) <- paste0("Prob_Stage_", 1:T_max)
#
#   # 'surv_from_top' tracks the cumulative probability of being in a LOWER stage
#   # It starts at 1.0 (100%) for everyone, representing P(Y <= T_max)
#   surv_from_top <- rep(1.0, N_new)
#
#   # 5. The Top-Down BCR Loop
#   # We iterate backward from the highest transition column down to the first
#   for (j in T_cols:1) {
#
#     # The actual stage class being evaluated is j + 1
#     current_stage <- j + 1
#
#     # P(Y = current_stage) = P(Y <= current_stage) * P(Transitioning here)
#     P_mat[, current_stage] <- surv_from_top * h[, j]
#
#     # Update the probability of being in an even lower stage for the next step
#     # P(Y <= current_stage - 1) = P(Y <= current_stage) * (1 - hazard)
#     surv_from_top <- surv_from_top * (1 - h[, j])
#   }
#
#   # 6. The Base Category
#   # Whatever probability is left over must belong to Stage 1
#   P_mat[, 1] <- surv_from_top
#
#   # 7. Classify based on highest probability
#   predicted_labels <- max.col(P_mat, ties.method = "first")
#
#   return(list(
#     probabilities = P_mat,
#     predictions   = predicted_labels
#   ))
# }


#' Calculate Fused Degrees of Freedom
#'
#' Computes the total degrees of freedom of a fitted fused coefficient matrix,
#' counting one baseline degree of freedom per intercept plus, for each gene,
#' the number of distinct non-zero fused blocks in its coefficient trajectory.
#'
#' @param beta_matrix Numeric coefficient matrix (P rows including intercept x
#'   TT columns of time transitions).
#'
#' @return A list with:
#' \describe{
#'   \item{total_df}{Numeric total degrees of freedom (baseline + gene DF).}
#'   \item{gene_df}{Numeric total penalized (gene-level) degrees of freedom.}
#'   \item{gene_dfs}{Numeric vector of degrees of freedom for each individual gene.}
#' }
#'
#' @export
calculate_fused_df <- function(beta_matrix) {
  TT <- ncol(beta_matrix)

  # 1. Baseline DF: Intercepts (Row 1)
  baseline_df <- TT

  penalized_df <- 0

  # NEW: Create a vector to store the exact DF for each individual gene
  gene_dfs <- numeric(nrow(beta_matrix) - 1)

  # 2. Penalized DF: Iterate through genes
  for (p in 2:nrow(beta_matrix)) {
    traj <- beta_matrix[p, ]

    # Check for true zeros
    if (all(traj == 0)) {
      gene_dfs[p - 1] <- 0
      next
    }

    # Count the unique structural blocks
    blocks <- rle(as.numeric(traj))$values
    df_local <- sum(blocks != 0)

    # Store local DF and add to total
    gene_dfs[p - 1] <- df_local
    penalized_df <- penalized_df + df_local
  }

  return(list(
    total_df = baseline_df + penalized_df,
    gene_df  = penalized_df,
    gene_dfs = gene_dfs # Return the vector of individual gene DFs
  ))
}

#' Sequential Cross-Validation for FCR ISTA Models
#'
#' Performs K-fold cross-validation sequentially (no parallel workers) over a
#' grid of lambda1/lambda2 values for the forward continuation ratio (FCR)
#' ISTA model, returning mean log-loss, misclassification rate, and log-loss
#' standard error for each lambda pair.
#'
#' @param X Numeric feature matrix (samples x genes), no intercept column.
#' @param Y Numeric or integer vector of ordinal outcomes.
#' @param cr_method Character; continuation ratio method. Only \code{"FCR"} is
#'   currently supported (BCR support is disabled).
#' @param lambda1_grid Numeric vector of candidate lambda1 (sparsity penalty) values.
#' @param lambda2_grid Numeric vector of candidate lambda2 (fusion penalty) values.
#' @param K_folds Integer number of cross-validation folds. Default 5.
#' @param standardize Logical; whether to standardize predictors within \code{run_CR_ISTA}. Default TRUE.
#' @param max_iter Integer maximum ISTA iterations. Default 1000.
#' @param tol Numeric convergence tolerance for ISTA. Default 1e-5.
#' @param seed Integer random seed for fold assignment. Default 123.
#'
#' @return A list with a single element \code{cv_results}, a data frame with
#'   columns \code{lambda1}, \code{lambda2}, \code{log_loss}, \code{misclass_rate},
#'   and \code{log_loss_se}, one row per lambda1/lambda2 combination.
#'
#' @importFrom foreach %do% foreach
#' @export
cv_ISTA_sequential <- function(X, Y, cr_method = c("FCR"), lambda1_grid, lambda2_grid, K_folds = 5,
                               standardize = TRUE, max_iter =
                                 1000, tol = 1e-5, seed = 123) {

  N <- nrow(X)
  cr_method <- match.arg(cr_method)

  if(cr_method == "FCR"){
    run_CR_ISTA = run_FCR_ISTA
    predict_CR = predict_FCR
  }
  # else if(cr_method == "BCR"){
  #   run_CR_ISTA = run_BCR_ISTA
  #   predict_CR = predict_BCR
  # }
  else{
    stop("cr_method should be populated with 'FCR'.") # either 'FCR' or 'BCR'."
  }

  #  Setup Data Folds
  set.seed(seed)
  fold_indices <- sample(rep(1:K_folds, length.out = N))

  #  Flatten the Task Grid
  task_grid <- expand.grid(lambda1 = lambda1_grid,
                           lambda2 = lambda2_grid,
                           fold = 1:K_folds)

  cat(sprintf("Starting sequential CV for %d total models...\n", nrow(task_grid)))

  #  Execute the Sequential Loop (%do%)
  sequential_results <- foreach(i = 1:nrow(task_grid), .combine = rbind) %do% {

    lam1 <- task_grid$lambda1[i]
    lam2 <- task_grid$lambda2[i]
    k    <- task_grid$fold[i]

    # Split the data for this fold
    test_idx  <- which(fold_indices == k)
    X_train   <- X[-test_idx, , drop = FALSE]
    Y_train   <- Y[-test_idx]
    X_test    <- X[test_idx, , drop = FALSE]
    Y_test    <- as.numeric(Y[test_idx])

    # Fit the model using our C++ engine
    fit <- run_CR_ISTA(X = X_train, Y = Y_train,
                       lambda1 = lam1, lambda2 = lam2,
                       standardize = standardize,
                       max_iter = max_iter, tol = tol)

    # Predict on the test set
    preds <- predict_CR(X_new = X_test, beta_matrix = fit$beta)

    # Calculate Log-Loss
    predicted_probs <- preds$probabilities
    true_class_probs <- predicted_probs[cbind(1:length(Y_test), Y_test)]
    true_class_probs <- pmax(pmin(true_class_probs, 1 - 1e-15), 1e-15)
    fold_log_loss <- -mean(log(true_class_probs))

    # Calculate Misclassification Rate
    fold_misclass <- mean(preds$predictions != Y_test)

    # Return metrics
    data.frame(lambda1 = lam1, lambda2 = lam2, fold = k,
               log_loss = fold_log_loss, misclass_rate = fold_misclass)
  }

  # 5. Aggregate the results across the K folds

  # A. Calculate the Mean Log-Loss and Misclassification Rate
  avg_results <- aggregate(cbind(log_loss, misclass_rate) ~ lambda1 + lambda2,
                           data = sequential_results,
                           FUN = mean)

  # B. Calculate the Standard Error of the Log-Loss
  # SE = standard deviation / sqrt(K)
  se_results <- aggregate(log_loss ~ lambda1 + lambda2,
                          data = sequential_results,
                          FUN = function(x) sd(x) / sqrt(K_folds))

  # Rename the column so it doesn't conflict during the merge
  names(se_results)[names(se_results) == "log_loss"] <- "log_loss_se"

  # C. Merge the Means and Standard Errors together
  final_results <- merge(avg_results, se_results, by = c("lambda1", "lambda2"))

  cat("Cross-Validation Complete.\n")

  return(list(
    cv_results = final_results
  ))
}

# #' Calculate Log-Likelihood for BCR
# #'
# #' Computes the (positive) log-likelihood of a fitted backward continuation
# #' ratio (BCR) model via the compiled \code{obj_g} C routine, using the label
# #' reversal trick so the FCR C++ engine can be reused.
# #'
# #' @param X Numeric feature matrix (samples x genes), no intercept column.
# #' @param Y Numeric or integer vector of ordinal outcomes/survival times.
# #' @param beta_matrix Numeric matrix of fitted BCR coefficients, including intercept row.
# #' @param standardize Logical; whether the model was fit with standardized predictors. Default FALSE.
# #'
# #' @return Numeric scalar log-likelihood value.
# #'
# #' @export
# calculate_loglik_BCR <- function(X, Y, beta_matrix, standardize = FALSE) {
#   N <- as.integer(nrow(X))
#   TT <- as.integer(ncol(beta_matrix))
#   T_max <- TT + 1
#
#   # 1. Add the intercept column
#   X_with_int <- cbind(1, X)
#   P_with_int <- as.integer(ncol(X_with_int))
#
#   # 2. THE REVERSAL HACK (Matching your ISTA wrapper)
#   Y_num <- as.numeric(Y)
#   Y_rev <- as.double(T_max - Y_num + 1)
#
#   # We must also temporarily reverse the columns of the beta matrix
#   # so it matches the reversed Y data when fed to the C++ engine
#   beta_rev <- beta_matrix[, TT:1, drop = FALSE]
#
#   # Calculate risk sets based on reversed labels
#   Nt <- as.integer(sapply(1:TT, function(t) sum(Y_rev >= t)))
#
#   std_flag <- as.integer(ifelse(standardize, 1, 0))
#
#   # 3. Call the unmodified FCR C++ likelihood function
#   result <- .C("obj_g",
#                YY = Y_rev,                      # Pass reversed Y
#                XX = as.double(X_with_int),
#                beta = as.double(beta_rev),      # Pass reversed Beta
#                pNN = N,
#                pPP = P_with_int,
#                pTT = TT,
#                Nt = Nt,
#                pStandardize = std_flag,
#                pGG = double(1))
#
#   # Return the true Log-Likelihood
#   return(-result$pGG)
# }
