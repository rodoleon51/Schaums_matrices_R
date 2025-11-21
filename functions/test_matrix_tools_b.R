library(pander)

source("functions/matrix_tools_b.R")

# Improve pander formatting
panderOptions('table.split.table', Inf)
panderOptions('round', 6)


# ------------------------------------------------------------
# Helper: Pretty print matrices using round_small() + pander()
# ------------------------------------------------------------
pander_clean <- function(X, tol = 1e-12, digits = 6) {
  Xc <- round_small(X, tol)
  Xc <- round(Xc, digits)
  pander(Xc)
}



# ------------------------------------------------------------
# 1. Test Householder QR
# ------------------------------------------------------------
test_householder_qr <- function() {
  cat("\n=== Test: Householder QR ===\n")
  
  A <- matrix(c(-4,2,2,
                3,-3,3,
                6,6,0), nrow = 3, byrow = TRUE)
  
  QR <- qr_decomposition(A, method = "householder")
  
  Q <- QR$Q
  R <- QR$R
  
  cat("\nMatrix A:\n"); pander_clean(A)
  
  cat("\nQ Matrix:\n"); pander_clean(Q)
  cat("\nR Matrix:\n"); pander_clean(R)
  
  cat("\nCheck QᵀQ = I:\n")
  pander_clean(t(Q) %*% Q)
  
  cat("\nCheck QR = A:\n")
  pander_clean(Q %*% R)
}



# ------------------------------------------------------------
# 2. Test Modified Gram–Schmidt QR
# ------------------------------------------------------------
test_mgs_qr <- function() {
  cat("\n=== Test: Modified Gram-Schmidt QR ===\n")
  
  A <- matrix(c(1,2,3,
                4,5,6,
                7,8,10), nrow = 3, byrow = TRUE)
  
  QR <- qr_decomposition(A, method = "mgs")
  
  Q <- QR$Q
  R <- QR$R
  
  cat("\nMatrix A:\n"); pander_clean(A)
  
  cat("\nQ Matrix:\n"); pander_clean(Q)
  cat("\nR Matrix:\n"); pander_clean(R)
  
  cat("\nCheck QᵀQ = I:\n")
  pander_clean(t(Q) %*% Q)
  
  cat("\nCheck QR = A:\n")
  pander_clean(Q %*% R)
}



# ------------------------------------------------------------
# 3. Test QR Algorithm (Eigenvalues)
# ------------------------------------------------------------
test_qr_algorithm <- function() {
  cat("\n=== Test: QR Algorithm (Eigenvalues) ===\n")
  
  A <- matrix(c(2,1,0,
                1,2,1,
                0,1,2), nrow = 3, byrow = TRUE)
  
  out <- qr_algorithm(A, diag = FALSE)
  
  cat("\nMatrix A:\n"); pander_clean(A)
  
  cat("\nEigenvalues (QR Algorithm):\n")
  pander_clean(out$values)
  
  cat("\nEigenvalues (base R eigen()):\n")
  pander_clean(eigen(A)$values)
  
  cat("\nFinal Matrix from QR Iterations:\n")
  pander_clean(out$matrix)
}

# ------------------------------------------------------------
# 4. Test Shifted QR Algorithm (Wilkinson shift)
# ------------------------------------------------------------
test_qr_algorithm_shifted <- function() {
  cat("\n=== Test: Shifted QR Algorithm (Wilkinson shift) ===\n")
  
  A <- matrix(c(2,1,0,
                1,2,1,
                0,1,2), nrow = 3, byrow = TRUE)
  
  out <- qr_algorithm_shifted(A, diag = FALSE)
  
  cat("\nMatrix A:\n"); pander_clean(A)
  
  cat("\nEigenvalues (Shifted QR Algorithm):\n")
  pander_clean(out$values)
  
  cat("\nEigenvalues (base R eigen()):\n")
  pander_clean(eigen(A)$values)
  
  cat("\nFinal Matrix from Shifted QR Iterations:\n")
  pander_clean(out$matrix)
}


# ------------------------------------------------------------
# Run All Tests
# ------------------------------------------------------------
run_all_tests <- function() {
  test_householder_qr()
  test_mgs_qr()
  test_qr_algorithm()
  test_qr_algorithm_shifted()
  
  cat("\nAll tests completed successfully.\n")
}

run_all_tests()
