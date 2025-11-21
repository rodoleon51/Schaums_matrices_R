# ------------------------------------------------------------
# lu_pivot_fr.R – EXACT LU with partial pivoting (fractions)
# FINAL VERIFIED VERSION — NO NULL PRINTING, CLEAN FRACTIONS
# ------------------------------------------------------------

if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
library(MASS)

# ============================================================
#  LU decomposition with partial pivoting (fraction form)
# ============================================================

lu_pivot_fr <- function(A) {
  n <- nrow(A)
  if (ncol(A) != n) stop("Matrix must be square")
  
  frac <- MASS::fractions
  
  # Convert A to fraction matrix
  M <- matrix(frac(A), n, n)       # working matrix for U
  L <- diag(n) * frac(1)           # unit lower triangular
  perm <- seq_len(n)               # permutation vector
  
  # ------------------------------------------------------------
  #                  ELIMINATION + PIVOTING
  # ------------------------------------------------------------
  for (k in seq_len(n-1)) {
    
    # --- Partial pivoting (largest abs value in column k) ---
    col_vals <- abs(M[k:n, k])
    max_row  <- which.max(col_vals) + k - 1
    
    if (max_row != k) {
      # Swap rows in U/M
      tmp        <- M[k, ];  M[k, ]  <- M[max_row, ];  M[max_row, ] <- tmp
      
      # Swap corresponding L row segments 1...(k-1)
      if (k > 1) {
        tmp                           <- L[k, 1:(k-1)]
        L[k, 1:(k-1)]                 <- L[max_row, 1:(k-1)]
        L[max_row, 1:(k-1)]           <- tmp
      }
      
      # Update permutation book-keeping
      perm[c(k, max_row)] <- perm[c(max_row, k)]
    }
    
    # --- Elimination ---
    for (i in (k+1):n) {
      L[i, k] <- M[i, k] / M[k, k]             # multiplier
      M[i, k:n] <- M[i, k:n] - L[i, k] * M[k, k:n]
    }
  }
  
  # ------------------------------------------------------------
  #              Final conversions to clean fractions
  # ------------------------------------------------------------
  U <- MASS::fractions(M)
  L <- MASS::fractions(L)
  P <- diag(n)[perm, ]
  
  list(
    P   = P,
    L   = L,
    U   = U,
    PA  = P %*% matrix(frac(A), n, n),
    perm = perm
  )
}

# ============================================================
#  Pretty-printer for LU-fr with optional pander() formatting
# ============================================================

print_lu_fr <- function(res, width = 14, tol = 1e-10, use_pander = FALSE) {
  
  cat("=== LU with Partial Pivoting – Exact Fractions ===\n\n")
  
  # ----------------------------------------------------------
  # Required converter for pander(): turns fraction objects
  # into plain character matrices like "1/3", "-2/5", etc.
  # ----------------------------------------------------------
  frac_to_char_matrix <- function(M) {
    matrix(as.character(M), nrow = nrow(M), ncol = ncol(M),
           dimnames = dimnames(M))
  }
  
  # ----------------------------------------------------------
  # PANDER MODE
  # ----------------------------------------------------------
  if (use_pander) {
    
    if (!requireNamespace("pander", quietly = TRUE)) {
      stop("pander package is required for use_pander = TRUE")
    }
    
    cat("Permutation matrix P:\n")
    pander::pander(res$P)
    
    cat("\nL (1s on diagonal):\n")
    pander::pander(frac_to_char_matrix(res$L))
    
    cat("\nU (pivots on diagonal):\n")
    pander::pander(frac_to_char_matrix(res$U))
    
    # Verification (always numeric → no fraction issue)
    diff <- res$PA - res$L %*% res$U
    diff_vec <- as.numeric(diff)
    
    cat("\nVerification: ")
    if (max(abs(diff_vec)) < tol) {
      cat("Exact match!\n")
    } else {
      cat("ERROR! (difference not zero)\n")
      pander::pander(diff)    # prints numeric matrix
    }
    
    return(invisible(res))
  }
  
  # ----------------------------------------------------------
  # PLAIN TEXT MODE (your original fixed-width printer)
  # ----------------------------------------------------------
  cat("Permutation matrix P:\n")
  print(res$P)
  
  safe_format <- function(x) {
    if (inherits(x, "fractions")) as.character(x) else as.character(x)
  }
  
  print_mat <- function(mat, name) {
    cat("\n", name, ":\n", sep = "")
    for (i in seq_len(nrow(mat))) {
      for (j in seq_len(ncol(mat))) {
        txt <- safe_format(mat[i, j])
        cat(sprintf("%*s", width, txt))
      }
      cat("\n")
    }
  }
  
  print_mat(res$L, "L (1s on diagonal)")
  print_mat(res$U, "U (pivots on diagonal)")
  
  # Verification
  diff <- res$PA - res$L %*% res$U
  diff_vec <- as.numeric(diff)
  
  cat("\nVerification: ")
  if (max(abs(diff_vec)) < tol) {
    cat("Exact match!\n")
  } else {
    cat("ERROR! (difference not zero)\n")
    print(diff)
  }
  
  invisible(res)
}


