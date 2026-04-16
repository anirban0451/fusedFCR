#include <R.h>
#include <R_ext/Print.h>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

// ------------------------------------------------------------------
// FORWARD DECLARATIONS (Table of Contents)
// ------------------------------------------------------------------

// Memory and Math Helpers
double ReadMat(double *XX, int *pNN, int *pDD, int i, int d);
void WriteMat(double *XX, int *pNN, int *pDD, int i, int d, double input);
double vecMlt(double *XX, double *YY, int NN);
double norm2(double *XX, int NN);

// Likelihood and Penalty Engines (Updated for Binary FCR)
extern "C" void obj_g(double *YY, double *XX, double *beta, int *pNN, int *pPP, int *pTT, int *Nt, int *pStandardize, double *pGG);
void delta(double *YY, double *XX, double *beta, int *pNN, int *pPP, int *pTT, int *Nt, int *pStandardize, double *resDelta);
void obj_h(int *pPP, int *pTT, double *beta, double *plambda1, double *plambda2, double* pHH);
void obj_fnc(double *YY, double *XX, double *beta, int *pNN, int *pPP, int *pTT, int *Nt, int *pStandardize, double *plambda1, double *plambda2, double *result);

// Tibshirani's 1D Fused Lasso Solvers
void DP_for_1DFL(int n, double *y, double lam, double *beta);
void soft_thresh(int n, double *y, double lam, double *beta);
void DP_for_S1DFL(int n, double *y, double lam1, double lam2, double *beta);

// Optimization Steps
void newBeta(double *beta, double* DELTA, double *tau, int *pNN, int *pPP, int *pTT, double *lambda1, double *lambda2, double *resBeta);
void backtracking(double *YY, double *XX, double* pTau, double* pFactor, double* beta, double* DELTA, double* lambda1, double* lambda2, int *pNN, int *pPP, int *pTT, int *pNt, int *pStandardize, double* resTau, double* new_Beta);

// ------------------------------------------------------------------

/* * ReadMat: Reads an element from a 2D matrix (passed from R as a 1D array)
 * * Memory Layout: Column-Major Order (R standard).
 * To access element (i, d) in an N x D matrix, we skip (d-1) full columns
 * of length N, and then step forward (i-1) rows.
 * * Variables:
 * XX   : Pointer to the 1D array representing the matrix
 * pNN  : Pointer to the number of rows (N)
 * pDD  : Pointer to the number of columns (D)
 * i    : Row index (1-based, 1 to N)
 * d    : Column index (1-based, 1 to D)
 */
double ReadMat(double*XX, int *pNN, int *pDD, int i, int d){
  double out=0;
  // Index calculation: (col - 1) * N + (row - 1)
  out=XX[(d-1)*(pNN[0])+i-1];
  return out;
}

/* * WriteMat: Writes an element to a 2D matrix (passed from R as a 1D array)
 * * Variables are identical to ReadMat.
 * 'input' is the scalar double value to be written to location (i, d).
 */
void WriteMat(double*XX, int *pNN, int *pDD, int i, int d, double input){
  // Index calculation: (col - 1) * N + (row - 1)
  XX[(d-1)*(pNN[0])+i-1]=input;
}

/* * ReadArr: Reads an element from a 3D array (passed from R as a 1D array)
 * * Memory Layout: 3D Column-Major Order.
 * Array dimensions are (N x D x K).
 * To access element (ii, jj, kk), we skip (kk-1) full 2D "slices" of size N*D,
 * then skip (jj-1) full columns of length N, and then step forward (ii-1) rows.
 * * Variables:
 * XX   : Pointer to the 1D array representing the 3D block
 * pNN  : Pointer to Dim 1 size (N - rows)
 * pDD  : Pointer to Dim 2 size (D - columns)
 * pKK  : Pointer to Dim 3 size (K - slices/depth)
 * ii   : Dim 1 index (1-based, 1 to N)
 * jj   : Dim 2 index (1-based, 1 to D)
 * kk   : Dim 3 index (1-based, 1 to K)
 */
double ReadArr(double*XX, int *pNN, int *pDD, int *pKK, int ii, int jj, int kk){
  double out=0;
  // Index calculation: (slice - 1) * (N * D) + (col - 1) * N + (row - 1)
  out=XX[(kk-1)*(pNN[0])*(pDD[0])+(jj-1)*(pNN[0])+ii-1];
  return out;
}

/* * WriteArr: Writes an element to a 3D array (passed from R as a 1D array)
 * * Variables are identical to ReadArr.
 * 'input' is the scalar double value to be written to location (ii, jj, kk).
 */
void WriteArr(double*XX, int *pNN, int *pDD, int *pKK, int ii, int jj, int kk, double input){
  // Index calculation: (slice - 1) * (N * D) + (col - 1) * N + (row - 1)
  XX[(kk-1)*(pNN[0])*(pDD[0])+(jj-1)*(pNN[0])+ii-1]=input;
}

/* * vecMlt: Computes the dot product (inner product) of two vectors
 * * Mathematical Formulation:
 * Result = Sum_{i=1}^{NN} (XX[i] * YY[i])
 * * Variables:
 * XX : Pointer to the first vector
 * YY : Pointer to the second vector
 * NN : Length of the vectors
 */
double vecMlt(double *XX, double *YY, int NN){
  double res=0;

  // Loop through elements 1 to NN.
  // We subtract 1 inside the loop to use native 0-based C++ indexing.
  for (int i=1; i<=NN; i++) {
    res += XX[i-1]*YY[i-1];
  }

  return res;
}

//norm2 returns the squared norm of a vector
double norm2(double *XX, int NN){
  double res=0;

  for (int i=1; i<=NN; i++) {
    res += XX[i-1]*XX[i-1];
  }

  return res;
}


extern "C" {
  /* * obj_g: Computes the Negative Log-Likelihood for the FCR Model
   * * Assumptions:
   * 1. Binary FCR (Independent Bernoulli trials per time point).
   * 2. No Censoring (Risk set at time t is simply everyone with YY >= t).
   * 3. Static Baseline Covariates (XX is a 2D matrix, not 3D).
   * * Mathematical Formulation:
   * We calculate the negative log-likelihood (binary cross-entropy)
   * over the risk set R(t) at each time point t:
   * * $$-\log L = - \sum_{t=1}^{T} \sum_{i: y_i \ge t} \left[ I(y_i = t) \eta_{it} - \log(1 + \exp(\eta_{it})) \right]$$
   * * * Where the linear predictor is:
   * $$\eta_{it} = \sum_{j=1}^{P} X_{ij} \beta_{jt}$$
   * * Inputs passed from R:
   * YY           : 1D Array (length N). YY[i-1] is the integer time subject i transitioned.
   * XX           : 2D Array (N * P). Baseline Predictors.
   * beta         : 2D Array (P * T). Time-varying parameters.
   * pNN          : Pointer to N (Total number of subjects).
   * pPP          : Pointer to P (Total number of predictors/genes).
   * pTT          : Pointer to T (Max time points to fit).
   * Nt           : 1D Array (length T). Nt[t-1] is the risk set size at time t.
   * pStandardize : Pointer to boolean flag (1 = divide by Nt, 0 = pure sum).
   * pGG          : Pointer to store the final scalar log-likelihood result.
   */
  void obj_g(double *YY, double *XX, double *beta, int *pNN, int *pPP, int *pTT, int *Nt, int *pStandardize, double *pGG){

    // Initialize the total objective value
    pGG[0] = 0.0;

    // Loop over each time point t = 1, ..., T
    for (int t=1; t<=pTT[0]; t++) {

      // Accumulator for the log-likelihood at the current time point t
      double t_sum = 0.0;

      // POINTER ALIASING (MEMORY OPTIMIZATION):
      // beta is stored in column-major order. To access the t-th column,
      // we skip (t-1) full columns of length P.
      // beta_t now points directly to the start of the parameters for time t.
      double* beta_t = &beta[(t-1) * pPP[0]];

      // Loop over each subject i = 1, ..., N
      for (int i=1; i<=pNN[0]; i++) {

        // RISK SET CONDITION: i \in R(t)
        if (YY[i-1] >= t) {

          // INDICATOR TARGET: I(y_i = t)
          double y_it = (YY[i-1] == t) ? 1.0 : 0.0;
          double eta = 0.0;

          // LINEAR PREDICTOR: \eta_{it} = X_i^\top \beta_t
          // Using the aliased pointer beta_t to bypass matrix indexing math
          for (int j=1; j<=pPP[0]; j++) {
            eta += ReadMat(XX, pNN, pPP, i, j) * beta_t[j-1];
          }

          // STABLE LOG-LIKELIHOOD CONTRIBUTION
          // Uses Log-Sum-Exp trick to prevent numerical overflow to +Inf
          if (eta > 0) {
            t_sum += (y_it * eta) - (eta + log(1.0 + exp(-eta)));
          } else {
            t_sum += (y_it * eta) - log(1.0 + exp(eta));
          }
        }
      }

      // CONDITIONAL STANDARDIZATION
      // If requested, calculate Mean Log-Loss. Otherwise, use Pure Joint Likelihood.
      if (pStandardize[0] == 1) {
        pGG[0] -= t_sum / Nt[t-1];
      } else {
        pGG[0] -= t_sum;
      }
    }
    return;
  }
}

/* * delta: Computes the Gradient of the Negative Log-Likelihood wrt Beta
 * * Mathematical Formulation:
 * We calculate the partial derivative for each parameter \beta_{jt}.
 * Using the chain rule, this is the sum of the error term multiplied by the covariate:
 * * $$\frac{\partial (-\log L)}{\partial \beta_{jt}} = \sum_{i: y_i \ge t} X_{ij} (p_{it} - I(y_i = t))$$
 * * * Where the predicted probability is:
 * $$p_{it} = \frac{\exp(\eta_{it})}{1 + \exp(\eta_{it})} = \frac{1}{1 + \exp(-\eta_{it})}$$
 * * Inputs passed from R:
 * resDelta : 2D Array (P * T) to store the computed gradients.
 * (All other inputs are identical to obj_g).
 */
void delta(double *YY, double *XX, double *beta, int *pNN, int *pPP, int *pTT, int *Nt, int *pStandardize, double *resDelta){

  // Step 1: Initialize the master gradient matrix (resDelta) to 0.0
  for(int j=1; j<=pPP[0]; j++){
    for(int t=1; t<=pTT[0]; t++){
      WriteMat(resDelta, pPP, pTT, j, t, 0.0);
    }
  }

  // Step 2: Compute the full gradient
  for (int t=1; t<=pTT[0]; t++) {

    // POINTER ALIASING (MEMORY OPTIMIZATION):
    // Point directly to the t-th columns of both the beta matrix and the gradient matrix.
    // This eliminates redundant ReadMat/WriteMat index calculations inside the loops.
    double* beta_t = &beta[(t-1) * pPP[0]];
    double* grad_t = &resDelta[(t-1) * pPP[0]];

    // Loop over subjects
    for (int i=1; i<=pNN[0]; i++) {

      // RISK SET CONDITION: i \in R(t)
      if (YY[i-1] >= t) {

        // INDICATOR TARGET: 1 if transitioned at t, 0 if survived
        double y_it = (YY[i-1] == t) ? 1.0 : 0.0;
        double eta = 0.0;

        // LINEAR PREDICTOR: \eta_{it} = X_i^\top \beta_t
        for (int p=1; p<=pPP[0]; p++) {
          eta += ReadMat(XX, pNN, pPP, i, p) * beta_t[p-1];
        }

        // BULLETPROOF PREDICTED PROBABILITY: p_{it}
        double prob = 0.0;
        if (eta > 0) {
          prob = 1.0 / (1.0 + exp(-eta));     // Prevents +Inf overflow
        } else {
          prob = exp(eta) / (1.0 + exp(eta)); // Prevents -Inf overflow
        }

        // ONE-SHOT ERROR TERM: (p_{it} - y_{it})
        // Pre-scale the error by 1/N_t here if standardization is requested
        double error_term = prob - y_it;
        if (pStandardize[0] == 1) {
          error_term /= Nt[t-1];
        }

        // COMPUTE GRADIENT CONTRIBUTION (Backward Pass)
        // Distribute the subject's error to all P genes simultaneously
        for(int j=1; j<=pPP[0]; j++){
          double x_ij = ReadMat(XX, pNN, pPP, i, j);

          // Directly accumulate into the memory address for \nabla \beta_{jt}
          grad_t[j-1] += x_ij * error_term;
        }
      }
    }
  }
  return;
}

/* * obj_h: Computes the Total Penalty Function for the Fused Lasso
 * * Mathematical Formulation:
 * Evaluates the L1 penalty (sparsity) and the Fused penalty (temporal smoothness).
 * $$h(\beta) = \lambda_1 \sum_{t=1}^{T}\sum_{j=2}^{P} |\beta_{jt}| + \lambda_2 \sum_{t=1}^{T-1}\sum_{j=2}^{P} |\beta_{j, t+1} - \beta_{j, t}|$$
 */
void obj_h(int *pPP, int *pTT, double *beta, double *plambda1, double *plambda2, double* pHH){

  pHH[0] = 0.0;

  // Part 1: Standard Lasso (L1) Loop
  for (int t=1; t<=pTT[0]; t++) {
    // POINTER ALIASING: Point directly to the t-th column
    double* beta_t = &beta[(t-1) * pPP[0]];

    for (int j=2; j<=pPP[0]; j++) {
      // Direct array access bypasses ReadMat
      pHH[0] += plambda1[0] * abs(beta_t[j-1]);
    }
  }

  // Part 2: Fused Lasso Loop
  for (int t=1; t<=(pTT[0]-1); t++) {
    // POINTER ALIASING: Point to column t, and column t+1
    double* beta_t = &beta[(t-1) * pPP[0]];
    double* beta_t_next = &beta[t * pPP[0]];

    for (int j=2; j<=pPP[0]; j++) {
      // Direct array access for the adjacent time difference
      pHH[0] += plambda2[0] * abs(beta_t_next[j-1] - beta_t[j-1]);
    }
  }

  return;
}

/* * obj_fnc: Computes the Global Objective Function
 * * Mathematical Formulation:
 * The total cost to minimize is the sum of the negative log-likelihood
 * and the regularization penalties:
 * * $$F(\beta) = -\log L(\beta) + h(\beta)$$
 * * Inputs:
 * Combines all inputs required for both obj_g (likelihood) and obj_h (penalty).
 * Includes the newly added pStandardize flag.
 * result : Pointer to store the final combined objective scalar.
 */
void obj_fnc(double *YY, double *XX, double *beta, int *pNN, int *pPP, int *pTT, int *Nt, int *pStandardize, double *plambda1, double *plambda2, double *result){

  // Create local variables to hold the separate likelihood and penalty values
  double GG;
  double HH;

  // Create pointers to these locals to pass into the sub-functions
  double* pGG = &GG;
  double* pHH = &HH;

  // Calculate the Penalty (h) and store it in HH
  obj_h(pPP, pTT, beta, plambda1, plambda2, pHH);

  // Calculate the Negative Log-Likelihood (g) and store it in GG
  // Passing the new pStandardize flag here
  obj_g(YY, XX, beta, pNN, pPP, pTT, Nt, pStandardize, pGG);

  // The final result is the sum of the likelihood cost and the penalty cost
  result[0] = pGG[0] + pHH[0];

  return;
}



// ------------------------------------------------------------------
// Ryan Tibshirani's 1D Fused Lasso Solvers
// These functions solve the exact proximal operators for L1 and TV penalties.
// ------------------------------------------------------------------

/* * DP_for_1DFL: Dynamic programming algorithm for the 1D Fused Lasso
 * * Mathematical Formulation:
 * Exact solver for the Total Variation (Fused) penalty.
 * Minimizes: $$ \frac{1}{2} \sum_{t=1}^{n} (y_t - \beta_t)^2 + \lambda \sum_{t=1}^{n-1} |\beta_{t+1} - \beta_t| $$
 * * Complexity: O(n)
 */
void DP_for_1DFL(int n, double *y, double lam, double *beta) {
  if (n==0) return;
  if (n==1 || lam==0) {
    for (int i=0; i<n; i++) beta[i] = y[i];
    return;
  }

  double afirst, alast, bfirst, blast;
  double *x = (double*)malloc(2*n*sizeof(double));
  double *a = (double*)malloc(2*n*sizeof(double));
  double *b = (double*)malloc(2*n*sizeof(double));
  int l,r;

  double *tm = (double*)malloc((n-1)*sizeof(double));
  double *tp = (double*)malloc((n-1)*sizeof(double));

  tm[0] = -lam+y[0];
  tp[0] = lam+y[0];
  l = n-1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = 1;
  b[l] = -y[0]+lam;
  a[r] = -1;
  b[r] = y[0]+lam;
  afirst = 1;
  bfirst = -lam-y[1];
  alast = -1;
  blast = -lam+y[1];

  int lo, hi;
  double alo, blo, ahi, bhi;
  for (int k=1; k<n-1; k++) {
    alo = afirst;
    blo = bfirst;
    for (lo=l; lo<=r; lo++) {
      if (alo*x[lo]+blo > -lam) break;
      alo += a[lo];
      blo += b[lo];
    }

    tm[k] = (-lam-blo)/alo;
    l = lo-1;
    x[l] = tm[k];

    ahi = alast;
    bhi = blast;
    for (hi=r; hi>=l; hi--) {
      if (-ahi*x[hi]-bhi < lam) break;
      ahi += a[hi];
      bhi += b[hi];
    }

    tp[k] = (lam+bhi)/(-ahi);
    r = hi+1;
    x[r] = tp[k];

    a[l] = alo;
    b[l] = blo+lam;
    a[r] = ahi;
    b[r] = bhi+lam;
    afirst = 1;
    bfirst = -lam-y[k+1];
    alast = -1;
    blast = -lam+y[k+1];
  }

  alo = afirst;
  blo = bfirst;
  for (lo=l; lo<=r; lo++) {
    if (alo*x[lo]+blo > 0) break;
    alo += a[lo];
    blo += b[lo];
  }
  beta[n-1] = -blo/alo;

  for (int k=n-2; k>=0; k--) {
    if (beta[k+1]>tp[k]) beta[k] = tp[k];
    else if (beta[k+1]<tm[k]) beta[k] = tm[k];
    else beta[k] = beta[k+1];
  }

  free(x); free(a); free(b); free(tm); free(tp);
}

/* * soft_thresh: Applies the standard L1 (Lasso) penalty
 * * Mathematical Formulation:
 * $$ \beta_i = \text{sign}(y_i) \max(|y_i| - \lambda, 0) $$
 */
void soft_thresh(int n, double *y, double lam, double *beta) {
  for (int i=0; i<n; i++) {
    if (y[i]>lam) beta[i] = y[i]-lam;
    else if (y[i]<-lam) beta[i] = y[i]+lam;
    else beta[i]=0;
  }
}

/* * DP_for_S1DFL: Solves the combined Sparse 1D Fused Lasso
 * * Applies the Fused penalty first, then applies the Lasso sparsity.
 */
void DP_for_S1DFL(int n, double *y, double lam1, double lam2, double *beta) {
  DP_for_1DFL(n,y,lam2,beta);
  soft_thresh(n,beta,lam1,beta);
}

/* * backtracking: Finds the optimal step size (tau) via ISTA line search
 * * Mathematical Formulation (Majorization-Minimization):
 * We evaluate a quadratic approximation of the smooth function g(beta) at the
 * current point. The step size tau (\tau) is valid if the actual likelihood at the
 * new proposed beta is less than or equal to this quadratic bound.
 * * Condition to satisfy: F_1 \le F_2
 * $$ F_1 = g(\beta^{(new)}) $$
 * $$ F_2 = g(\beta) + \langle \nabla g(\beta), \beta^{(new)} - \beta \rangle + \frac{1}{2\tau} \|\beta^{(new)} - \beta\|_2^2 $$
 * * If F_1 > F_2, we shrink \tau by a factor (e.g., \tau = 0.8 \times \tau) and try again.
 */
void backtracking(double *YY, double *XX, double* pTau, double* pFactor, double* beta, double* DELTA, double* lambda1, double* lambda2, int *pNN, int *pPP, int *pTT, int *pNt, int *pStandardize, double* resTau, double* new_Beta ){

  double lam1, lam2;
  double* pLam1 = &lam1;
  double* pLam2 = &lam2;

  // F1 evaluates the true objective, F2 evaluates the quadratic bound
  double F1 = 1.0;      double* pF1 = &F1;
  double F2 = 0.0;      double* pF2 = &F2;
  double firstF2 = 0.0; double* pFirstF2 = &firstF2;

  // Array to hold the difference vector: G = (\beta - \beta^{(new)}) / \tau
  double* GG = new double[pPP[0] * pTT[0]];

  // Initialize tau and scale the penalty tuning parameters accordingly
  resTau[0] = pTau[0];
  lam1 = resTau[0] * lambda1[0];
  lam2 = resTau[0] * lambda2[0];

  int b = 1;
  // Loop until the quadratic bound condition (F1 <= F2) is met, max 50 iterations
  while(pF1[0] > pF2[0] && b < 50){
    b = b + 1;

    // 1. Propose a new beta using the Proximal Gradient Step
    // \beta^{(new)} = \text{prox}_{\tau h} (\beta - \tau \nabla g(\beta))
    newBeta(beta, DELTA, resTau, pNN, pPP, pTT, pLam1, pLam2, new_Beta);

    // 2. Compute the scaled difference vector GG
    for (int i=1; i<=(pPP[0]*pTT[0]); i++) {
      GG[i-1] = (beta[i-1] - new_Beta[i-1]) / resTau[0];
    }

    // 3. Compute F1: The true negative log-likelihood at the proposed \beta^{(new)}
    obj_g(YY, XX, new_Beta, pNN, pPP, pTT, pNt, pStandardize, pF1);

    // 4. Compute F2: The quadratic approximation at the old \beta
    // First get g(\beta)
    obj_g(YY, XX, beta, pNN, pPP, pTT, pNt, pStandardize, pFirstF2);

    // F2 = g(\beta) - \tau * (\nabla g(\beta) \cdot GG) + (\tau / 2) * ||GG||_2^2
    // Note: vecMlt and norm2 operate on the flattened P*T arrays
    F2 = pFirstF2[0]
    - resTau[0] * vecMlt(DELTA, GG, pPP[0]*pTT[0])
      + (resTau[0] / 2.0) * norm2(GG, pPP[0]*pTT[0]);

      // 5. Shrink step size if condition fails
      resTau[0] = pFactor[0] * resTau[0];
      lam1 = resTau[0] * lambda1[0];
      lam2 = resTau[0] * lambda2[0];
  }

  // The loop aggressively shrinks tau. We revert the final multiplication
  // to return the optimal tau that satisfied the condition.
  resTau[0] = resTau[0] / pFactor[0];

  delete[] GG;
}


/* * newBeta: Computes the Proximal Gradient Step
 * * Mathematical Formulation:
 * $$ \beta^{(new)} = \text{prox}_{\tau h} \left( \beta - \tau \nabla g(\beta) \right) $$
 * * Memory Note: We extract "rows" across time for each gene j. Because data is column-major,
 * we jump memory addresses by P to get the next time point for the same gene.
 */
void newBeta(double *beta, double* DELTA, double *tau, int *pNN, int *pPP, int *pTT, double *lambda1, double *lambda2, double *resBeta){

  double *yy = new double[pTT[0]];
  double *tBeta = new double[pTT[0]];

  // Loop over genes (j starts at 2, intercept is separate)
  for (int j=2; j<=pPP[0]; j++) {

    // 1. Take the standard gradient descent step across all time points t
    for (int t=1; t<=(pTT[0]); t++) {
      // DIRECT 1D INDEXING: (col - 1) * N_rows + (row - 1)
      int flat_index = (t-1) * pPP[0] + (j-1);

      // Bypass ReadMat entirely
      yy[t-1] = beta[flat_index] - tau[0] * DELTA[flat_index];
    }

    // 2. Apply the Proximal Operator
    DP_for_S1DFL(pTT[0], yy, lambda1[0], lambda2[0], tBeta);

    // 3. Save the penalized parameters back to the result matrix
    for (int t=1; t<=(pTT[0]); t++) {
      int flat_index = (t-1) * pPP[0] + (j-1);

      // Bypass WriteMat entirely
      resBeta[flat_index] = tBeta[t-1];
    }
  }

  // For the intercept (j=1), manually update WITHOUT penalties
  for (int t=1; t<=(pTT[0]); t++) {
    int flat_index = (t-1) * pPP[0] + 0; // j=1 means row index 0

    resBeta[flat_index] = beta[flat_index] - tau[0] * DELTA[flat_index];
  }

  delete[] yy;
  delete[] tBeta;
}



// R wrapper

extern "C" {

  /* * GGD (Generalized Gradient Descent): Main ISTA Wrapper for Binary FCR
   * * * Memory / Storing Strategy Note:
   * beta, new_Beta, DELTA, and betaCriterion mathematically represent P x T matrices,
   * but because of R's column-major format, they are passed and stored as flat 1D
   * arrays of size P * T. We iterate over them from 0 to (P*T - 1) for fast
   * contiguous memory block operations, entirely bypassing ReadMat/WriteMat here.
   * * * Inputs passed from R via .C():
   * YY           : 1D Array (length N). Transition times.
   * XX           : 1D Array (N * P). Baseline Predictors.
   * beta         : 1D Array (P * T). Initial parameter estimates.
   * lambda1      : L1 (Lasso) penalty tuning parameter.
   * lambda2      : 1D Total Variation (Fused Lasso) penalty tuning parameter.
   * pNN          : Pointer to N (Number of subjects).
   * pPP          : Pointer to P (Number of predictors + 1 for intercept).
   * pTT          : Pointer to T (Max time points).
   * pNt          : 1D Array (length T). Risk set sizes at each time point.
   * pStandardize : Pointer to boolean flag (1 = Mean Log-Loss, 0 = Pure Sum).
   * pNiter       : Maximum number of iterations to run.
   * pStop        : Stopping criterion (1 = Objective change, 2 = Beta change).
   * pEps         : Tolerance threshold for the stopping criterion.
   * pTau         : Initial step size for the backtracking line search.
   * pFactor      : Shrinkage factor for backtracking (e.g., 0.8).
   * new_Beta     : 1D Array (P * T). Stores the updated parameters for the next step.
   * objectives   : 1D Array to track objective history for R output.
   * taus         : 1D Array to track step size history for R output.
   * criterions   : 1D Array to track stopping criterion history for R output.
   * actualNiter  : Pointer to return the exact number of iterations ran.
   */
  void GGD(double *YY, double *XX, double* beta, double* lambda1, double* lambda2,
           int *pNN, int *pPP, int *pTT, int *pNt, int *pStandardize,
           int* pNiter, int* pStop, double* pEps, double* pTau, double* pFactor,
           double* new_Beta, double* objectives, double* taus, double* criterions, int *actualNiter){

    // Total flat memory size of the P x T matrices
    int PT_size = pPP[0] * pTT[0];

    // Allocate local memory arrays for the Gradient and Convergence Tracker
    double* DELTA = new double[PT_size];
    double* betaCriterion = new double[PT_size];

    double stopping = pEps[0] + 1.0;
    double resTau = 0.0;
    double* presTau = &resTau;
    double beforeF = 0.0;
    double* pbeforeF = &beforeF;
    double afterF = 0.0;
    double* pafterF = &afterF;

    // Iteration step counter
    int ss = 1;

    // Initialize new_Beta with the starting beta values provided by R
    for (int i = 0; i < PT_size; i++) {
      new_Beta[i] = beta[i];
    }

    // ------------------------------------------------------------------
    // MAIN ISTA ALGORITHM LOOP
    // ------------------------------------------------------------------
    while (ss <= pNiter[0] && stopping > pEps[0]) {

      // 1. UPDATE STATE: The new beta becomes the old beta for the next step
      for (int i = 0; i < PT_size; i++) {
        beta[i] = new_Beta[i];
      }

      // 2. COMPUTE GRADIENT: Delta = nabla g(beta)
      // Evaluates the full gradient and stores it in the DELTA array
      delta(YY, XX, beta, pNN, pPP, pTT, pNt, pStandardize, DELTA);

      // 3. STORE OBJECTIVE BEFORE STEP: g(beta) + h(beta)
      if (ss == 1) {
        obj_fnc(YY, XX, beta, pNN, pPP, pTT, pNt, pStandardize, lambda1, lambda2, pbeforeF);
      } else {
        // Recycle the evaluated value from the end of the previous iteration
        pbeforeF[0] = pafterF[0];
      }

      // 4. BACKTRACKING LINE SEARCH & PROXIMAL STEP
      // Finds optimal tau and computes Beta_new = prox_{tau * h}(beta - tau * Delta)
      backtracking(YY, XX, pTau, pFactor, beta, DELTA, lambda1, lambda2, pNN, pPP, pTT, pNt, pStandardize, presTau, new_Beta);

      // 5. STORE OBJECTIVE AFTER STEP: g(Beta_new) + h(Beta_new)
      obj_fnc(YY, XX, new_Beta, pNN, pPP, pTT, pNt, pStandardize, lambda1, lambda2, pafterF);

      // 6. EVALUATE STOPPING CRITERION
      if (pStop[0] == 1) {
        // Criterion 1: Relative change in the Objective Function
        // stopping = |F(Beta_new) - F(beta)| / F(beta)
        stopping = abs(pafterF[0] - pbeforeF[0]) / pbeforeF[0];
      }
      else if (pStop[0] == 2) {
        // Criterion 2: Relative L2 Norm change in the Parameter Matrix
        // stopping = ||Beta_new - beta||_2 / ||beta||_2
        for (int i = 0; i < PT_size; i++) {
          betaCriterion[i] = new_Beta[i] - beta[i];
        }
        stopping = (sqrt(norm2(betaCriterion, PT_size))) / sqrt(norm2(beta, PT_size));
      }
      else {
        // Criterion 3: Force the algorithm to run for the maximum number of iterations
        stopping = pEps[0] + 1.0;
      }

      // 7. RECORD METRICS for R output analysis
      taus[ss-1] = presTau[0];
      objectives[ss-1] = pafterF[0];
      criterions[ss-1] = stopping;

      ss = ss + 1;
    }

    // Record the actual number of iterations completed before convergence
    actualNiter[0] = ss - 1;

    // Free allocated RAM to prevent memory leaks
    delete[] DELTA;
    delete[] betaCriterion;
  }
}
