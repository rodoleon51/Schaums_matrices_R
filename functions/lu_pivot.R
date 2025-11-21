# LU Decomposition with Partial Pivoting (Doolittle form)
lu_pivot <- function(A, tol = 1e-12) {
  n <- nrow(A)
  if (ncol(A) != n) stop("Matrix must be square")
  
  M <- A
  L <- diag(n)              # Unit lower triangular
  P <- diag(n)              # Permutation matrix
  pivots <- 1:n             # Track row swaps
  
  for (k in 1:(n-1)) {
    # ---- Partial pivoting ----
    max_idx <- k
    max_val <- abs(M[k, k])
    for (i in (k+1):n) {
      if (abs(M[i, k]) > max_val) {
        max_val <- abs(M[i, k])
        max_idx <- i
      }
    }
    
    if (max_val < tol) stop(paste("Matrix is singular (pivot", k, "<", tol, ")"))
    
    if (max_idx != k) {
      temp <- M[k, ]; M[k, ] <- M[max_idx, ]; M[max_idx, ] <- temp
      if (k > 1) {
        temp <- L[k, 1:(k-1)]; L[k, 1:(k-1)] <- L[max_idx, 1:(k-1)]; L[max_idx, 1:(k-1)] <- temp
      }
      pivots[c(k, max_idx)] <- pivots[c(max_idx, k)]
    }
    
    # ---- Elimination ----
    for (i in (k+1):n) {
      L[i, k] <- M[i, k] / M[k, k]
      M[i, k:n] <- M[i, k:n] - L[i, k] * M[k, k:n]
    }
  }
  
  U <- M
  P <- P[pivots, ]
  
  list(P = P, L = L, U = U, PA = P %*% A, permutation = pivots)
}

# Pretty print helper
print_lu <- function(lu_res, digits = 6, tol = 1e-10) {
  cat("\n\n=== LU with Partial Pivoting - Decimals ===\n")
  cat("P (Permutation):\n"); pander(round(lu_res$P, digits))
  
  cat("\n\nL (Lower, 1s on diagonal):\n"); pander(round(lu_res$L, digits))
  
  cat("\n\nU (Upper, pivots on diagonal):\n"); pander(round(lu_res$U, digits))
  
  cat("\n\nVerification: ")
  diff_mat <- lu_res$PA - lu_res$L %*% lu_res$U
  diff_vec <- as.numeric(diff_mat)
  if (max(abs(diff_vec)) < tol) {
    cat("Exact match!\n")
  } else {
    cat("ERROR! (difference not zero)\n")
    print(diff_mat)
  }
}
