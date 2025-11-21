############################################################
# MATRIX TOOLS LIBRARY
# Crout & Doolittle (pivoted & non-pivoted)
# Determinants, Inverses, Solvers, REF/RREF, Diagnostics
############################################################

############################################################
# DIAGNOSTIC & TRIANGULAR HELPERS
############################################################

diag_msg <- function(...) {
  cat("[diag]", ..., "\n")
}

is_upper_tri_matrix <- function(A) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) return(FALSE)
  all(A[lower.tri(A)] == 0)
}

is_lower_tri_matrix <- function(A) {
  A <- as.matrix(A)
  if (nrow(A) != ncol(A)) return(FALSE)
  all(A[upper.tri(A)] == 0)
}


############################################################
# BASIC UTILITIES: FORWARD/BACK SUBSTITUTION, MATVEC
############################################################

forward_substitution <- function(L, b) {
  n <- length(b)
  y <- numeric(n)
  diag_msg("forward_substitution(): n =", n)
  for (i in 1:n) {
    s <- if (i > 1) sum(L[i, 1:(i - 1)] * y[1:(i - 1)]) else 0
    y[i] <- (b[i] - s) / L[i, i]
  }
  y
}

backward_substitution <- function(U, y) {
  n <- length(y)
  x <- numeric(n)
  diag_msg("backward_substitution(): n =", n)
  for (i in n:1) {
    s <- if (i < n) sum(U[i, (i + 1):n] * x[(i + 1):n]) else 0
    x[i] <- (y[i] - s) / U[i, i]
  }
  x
}

matvec <- function(A, x) {
  diag_msg("matvec(): computing A %*% x")
  A %*% x
}


############################################################
# ==========================================================
#               C R O U T   M E T H O D S
# ==========================================================
############################################################

############################################################
# 1. Crout LU (no pivoting)  A = L U, U has 1s on diag
############################################################

crout_decomposition <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  L <- matrix(0, n, n)
  U <- diag(n)
  
  diag_msg("crout_decomposition(): starting, n =", n)
  
  for (k in 1:n) {
    # Compute L[i,k]
    for (i in k:n) {
      s <- if (k > 1) sum(L[i, 1:(k - 1)] * U[1:(k - 1), k]) else 0
      L[i, k] <- A[i, k] - s
    }
    
    if (L[k, k] == 0) {
      diag_msg("crout_decomposition(): zero pivot at k =", k)
      stop("Zero pivot in non-pivoted Crout.")
    }
    
    # Compute U[k,j] for j > k
    if (k < n) {
      for (j in (k + 1):n) {
        s <- if (k > 1) sum(L[k, 1:(k - 1)] * U[1:(k - 1), j]) else 0
        U[k, j] <- (A[k, j] - s) / L[k, k]
      }
    }
  }
  
  list(L = L, U = U)
}


############################################################
# 2. Crout solve Ax = b (no pivoting)
############################################################

solve_lu <- function(A, b, L = NULL, U = NULL) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  n <- nrow(A)
  
  diag_msg("solve_lu(): inspecting A of size", nrow(A), "x", ncol(A))
  
  if (is_upper_tri_matrix(A) || is_lower_tri_matrix(A)) {
    kind <- if (is_upper_tri_matrix(A)) "upper" else "lower"
    diag_msg("solve_lu(): A is", kind, "triangular — using base::solve(A, b)")
    return(as.vector(solve(A, b)))
  }
  
  if (is.null(L) || is.null(U)) {
    lu <- crout_decomposition(A)
    L <- lu$L; U <- lu$U
  }
  
  y <- forward_substitution(L, b)
  x <- backward_substitution(U, y)
  x
}


############################################################
# 3. Crout determinant
############################################################

det_lu <- function(A) {
  A <- as.matrix(A)
  if (is_upper_tri_matrix(A) || is_lower_tri_matrix(A)) {
    diag_msg("det_lu(): A is triangular — using product of diagonal")
    return(prod(diag(A)))
  }
  
  lu <- crout_decomposition(A)
  diag_msg("det_lu(): determinant via prod(diag(L))")
  prod(diag(lu$L))
}


############################################################
# 4. Crout inverse
############################################################

matrix_inverse <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  
  if (is_upper_tri_matrix(A) || is_lower_tri_matrix(A)) {
    diag_msg("matrix_inverse(): A is triangular — using base::solve(A)")
    return(solve(A))
  }
  
  diag_msg("matrix_inverse(): using Crout LU")
  lu <- crout_decomposition(A)
  L <- lu$L; U <- lu$U
  I <- diag(n)
  invA <- matrix(0, n, n)
  
  for (i in 1:n) {
    y <- forward_substitution(L, I[, i])
    invA[, i] <- backward_substitution(U, y)
  }
  
  invA
}


############################################################
# 5. Crout with Partial Pivoting (PA = L U)
############################################################

crout_decomposition_pivoting <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  
  P <- diag(n)
  L <- matrix(0, n, n)
  U <- diag(n)
  
  diag_msg("crout_decomposition_pivoting(): starting, n =", n)
  
  for (k in 1:n) {
    pivot <- which.max(abs(A[k:n, k])) + (k - 1)
    if (A[pivot, k] == 0) {
      diag_msg("crout_decomposition_pivoting(): singular matrix at column", k)
      stop("Matrix is singular in pivoted Crout.")
    }
    
    if (pivot != k) {
      diag_msg("crout_decomposition_pivoting(): swapping rows", k, "<->", pivot)
      A[c(k, pivot), ] <- A[c(pivot, k), ]
      P[c(k, pivot), ] <- P[c(pivot, k), ]
      if (k > 1)
        L[c(k, pivot), 1:(k - 1)] <- L[c(pivot, k), 1:(k - 1)]
    }
    
    for (i in k:n) {
      s <- if (k > 1) sum(L[i, 1:(k - 1)] * U[1:(k - 1), k]) else 0
      L[i, k] <- A[i, k] - s
    }
    
    if (L[k, k] == 0) {
      diag_msg("crout_decomposition_pivoting(): zero pivot at k =", k)
      stop("Zero pivot in pivoted Crout.")
    }
    
    if (k < n) {
      for (j in (k + 1):n) {
        s <- if (k > 1) sum(L[k, 1:(k - 1)] * U[1:(k - 1), j]) else 0
        U[k, j] <- (A[k, j] - s) / L[k, k]
      }
    }
  }
  
  list(P = P, L = L, U = U)
}


############################################################
# ==========================================================
#            D O O L I T T L E   M E T H O D S
# ==========================================================
############################################################

############################################################
# 1. Doolittle LU (no pivoting)
############################################################

doolittle_decomposition <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  
  L <- matrix(0, n, n)
  U <- diag(n)
  
  diag_msg("doolittle_decomposition(): starting, n =", n)
  
  for (k in 1:n) {
    # Compute L[i,k]
    for (i in k:n) {
      s <- if (k > 1) sum(L[i, 1:(k - 1)] * U[1:(k - 1), k]) else 0
      L[i, k] <- A[i, k] - s
    }
    
    if (L[k, k] == 0) {
      diag_msg("doolittle_decomposition(): zero pivot at k =", k,
               "- consider doolittle_decomposition_pivoting().")
      stop("Zero pivot in non-pivoted Doolittle.")
    }
    
    # Compute U[k,j]
    if (k < n) {
      for (j in (k + 1):n) {
        s <- if (k > 1) sum(L[k, 1:(k - 1)] * U[1:(k - 1), j]) else 0
        U[k, j] <- (A[k, j] - s) / L[k, k]
      }
    }
  }
  
  list(L = L, U = U)
}


############################################################
# 2. Doolittle solver (no pivoting)
############################################################

solve_doolittle <- function(A, b) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  
  diag_msg("solve_doolittle(): inspecting A", nrow(A), "x", ncol(A))
  
  if (is_upper_tri_matrix(A) || is_lower_tri_matrix(A)) {
    kind <- if (is_upper_tri_matrix(A)) "upper" else "lower"
    diag_msg("solve_doolittle(): A is", kind, "triangular — using base::solve(A, b)")
    return(as.vector(solve(A, b)))
  }
  
  lu <- doolittle_decomposition(A)
  L <- lu$L; U <- lu$U
  y <- forward_substitution(L, b)
  x <- backward_substitution(U, y)
  x
}


############################################################
# 3. Doolittle determinant (no pivoting)
############################################################

det_doolittle <- function(A) {
  A <- as.matrix(A)
  if (is_upper_tri_matrix(A) || is_lower_tri_matrix(A)) {
    diag_msg("det_doolittle(): A is triangular — using product of diagonal")
    return(prod(diag(A)))
  }
  
  lu <- doolittle_decomposition(A)
  diag_msg("det_doolittle(): determinant via prod(diag(L))")
  prod(diag(lu$L))
}


############################################################
# 4. Doolittle inverse (no pivoting)
############################################################

inverse_doolittle <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  
  if (is_upper_tri_matrix(A) || is_lower_tri_matrix(A)) {
    diag_msg("inverse_doolittle(): A is triangular — using base::solve(A)")
    return(solve(A))
  }
  
  diag_msg("inverse_doolittle(): using non-pivoted Doolittle")
  I <- diag(n)
  invA <- matrix(0, n, n)
  
  for (i in 1:n) {
    invA[, i] <- solve_doolittle(A, I[, i])
  }
  
  invA
}


############################################################
# 5. Doolittle with Partial Pivoting (PA = L U)
############################################################

doolittle_decomposition_pivoting <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  
  P <- diag(n)
  L <- matrix(0, n, n)
  U <- diag(n)
  
  diag_msg("doolittle_decomposition_pivoting(): starting, n =", n)
  
  for (k in 1:n) {
    pivot <- which.max(abs(A[k:n, k])) + (k - 1)
    if (A[pivot, k] == 0) {
      diag_msg("doolittle_decomposition_pivoting(): singular matrix at column", k)
      stop("Matrix is singular in pivoted Doolittle.")
    }
    
    if (pivot != k) {
      diag_msg("doolittle_decomposition_pivoting(): swapping rows", k, "<->", pivot)
      A[c(k, pivot), ] <- A[c(pivot, k), ]
      P[c(k, pivot), ] <- P[c(pivot, k), ]
      if (k > 1)
        L[c(k, pivot), 1:(k - 1)] <- L[c(pivot, k), 1:(k - 1)]
    }
    
    for (i in k:n) {
      s <- if (k > 1) sum(L[i, 1:(k - 1)] * U[1:(k - 1), k]) else 0
      L[i, k] <- A[i, k] - s
    }
    
    if (L[k, k] == 0) {
      diag_msg("doolittle_decomposition_pivoting(): zero pivot at k =", k)
      stop("Zero pivot even after pivoting (Doolittle).")
    }
    
    if (k < n) {
      for (j in (k + 1):n) {
        s <- if (k > 1) sum(L[k, 1:(k - 1)] * U[1:(k - 1), j]) else 0
        U[k, j] <- (A[k, j] - s) / L[k, k]
      }
    }
  }
  
  list(P = P, L = L, U = U)
}


############################################################
# 6. Doolittle solver with pivoting
############################################################

solve_doolittle_pivot <- function(A, b) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  
  diag_msg("solve_doolittle_pivot(): using pivoted Doolittle")
  
  if (is_upper_tri_matrix(A) || is_lower_tri_matrix(A)) {
    kind <- if (is_upper_tri_matrix(A)) "upper" else "lower"
    diag_msg("solve_doolittle_pivot(): A is", kind,
             "triangular — using base::solve(A, b)")
    return(as.vector(solve(A, b)))
  }
  
  lu <- doolittle_decomposition_pivoting(A)
  P <- lu$P; L <- lu$L; U <- lu$U
  b2 <- P %*% b
  y <- forward_substitution(L, b2)
  x <- backward_substitution(U, y)
  x
}


############################################################
# 7. Doolittle determinant with pivoting
############################################################

det_doolittle_pivot <- function(A) {
  A <- as.matrix(A)
  
  if (is_upper_tri_matrix(A) || is_lower_tri_matrix(A)) {
    diag_msg("det_doolittle_pivot(): A is triangular — using product of diagonal")
    return(prod(diag(A)))
  }
  
  lu <- doolittle_decomposition_pivoting(A)
  P <- lu$P; L <- lu$L
  sign <- ifelse(det(P) < 0, -1, 1)
  diag_msg("det_doolittle_pivot(): determinant via sign(P) * prod(diag(L))")
  sign * prod(diag(L))
}


############################################################
# 8. Doolittle inverse with pivoting
############################################################

inverse_doolittle_pivot <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  
  if (is_upper_tri_matrix(A) || is_lower_tri_matrix(A)) {
    diag_msg("inverse_doolittle_pivot(): A is triangular — using base::solve(A)")
    return(solve(A))
  }
  
  diag_msg("inverse_doolittle_pivot(): using pivoted Doolittle")
  I <- diag(n)
  invA <- matrix(0, n, n)
  
  for (i in 1:n) {
    invA[, i] <- solve_doolittle_pivot(A, I[, i])
  }
  
  invA
}


############################################################
# ==========================================================
#                 REF / RREF & RANK
# ==========================================================
############################################################

row_echelon_form <- function(Alist, name, form = c("ref", "rref", "both")) {
  form <- match.arg(form)
  A <- as.matrix(Alist[[name]])
  rows <- nrow(A)
  cols <- ncol(A)
  r <- 1
  pivots <- integer(0)
  
  diag_msg("row_echelon_form(): form =", form, ", size =", rows, "x", cols)
  
  # Build REF
  for (c in 1:cols) {
    if (r > rows) break
    
    pivot_row <- which(A[r:rows, c] != 0)[1]
    if (!is.na(pivot_row)) {
      pivot_row <- pivot_row + r - 1
      if (pivot_row != r) {
        diag_msg("row_echelon_form(): swap rows", r, "<->", pivot_row)
        A[c(r, pivot_row), ] <- A[c(pivot_row, r), ]
      }
      
      pivots <- c(pivots, c)
      
      # Normalize pivot row
      A[r, ] <- A[r, ] / A[r, c]
      
      # Eliminate below
      if (r < rows) {
        for (i in (r + 1):rows) {
          if (A[i, c] != 0) {
            A[i, ] <- A[i, ] - A[i, c] * A[r, ]
          }
        }
      }
      
      r <- r + 1
    }
  }
  
  REF <- A
  
  if (form == "ref") {
    return(REF)
  }
  
  # Build RREF from REF
  RREF <- REF
  if (length(pivots) > 0) {
    for (p in seq_along(pivots)) {
      c <- pivots[p]
      pivot_row <- which(RREF[, c] == 1)[1]
      if (!is.na(pivot_row)) {
        if (pivot_row > 1) {
          for (i in 1:(pivot_row - 1)) {
            if (RREF[i, c] != 0) {
              RREF[i, ] <- RREF[i, ] - RREF[i, c] * RREF[pivot_row, ]
            }
          }
        }
      }
    }
  }
  
  if (form == "rref") {
    return(RREF)
  }
  
  # BOTH
  list(REF, RREF, pivots)
}

rank_lu <- function(A) {
  diag_msg("rank_lu(): computing rank via REF")
  REF <- row_echelon_form(list(A = A), "A", form = "ref")
  # Rank = number of non-zero rows
  rank <- sum(apply(REF, 1, function(row) any(abs(row) > 1e-12)))
  rank
}


############################################################
# ==========================================================
#                 FULL ANALYSIS WRAPPER
# ==========================================================
############################################################

run_full_analysis <- function(A, b = NULL, name = "Matrix") {
  library(pander)
  panderOptions("digits", 6)
  
  section <- function(title) {
    cat("\n\n")
    pander(paste0("## ", title))
    cat("\n")
  }
  
  A <- as.matrix(A)
  
  section(paste("Analysis for", name))
  pander("### Input Matrix A")
  pander(A)
  
  if (!is.null(b)) {
    pander("### Input Vector b")
    pander(b)
  }
  
  # Triangular diagnostics
  isU <- is_upper_tri_matrix(A)
  isL <- is_lower_tri_matrix(A)
  if (isU || isL) {
    kind <- if (isU) "upper triangular" else "lower triangular"
    diag_msg("run_full_analysis(): A is detected as", kind)
    pander(paste("Diagnostic: A is", kind, "— some operations use direct solve()."))
  } else {
    diag_msg("run_full_analysis(): A is not triangular")
  }
  
  ## ---------------- Crout (no pivot) ----------------
  section("Crout LU Decomposition (no pivot)")
  lu_c <- crout_decomposition(A)
  pander("L (Crout)"); pander(lu_c$L)
  pander("U (Crout)"); pander(lu_c$U)
  
  if (!is.null(b)) {
    section("Solve Ax = b (Crout, no pivot)")
    x_c <- solve_lu(A, b)
    pander("x"); pander(x_c)
    pander("A %*% x"); pander(A %*% x_c)
  }
  
  section("det(A) via Crout")
  pander(det_lu(A))
  
  section("Inverse(A) via Crout")
  inv_c <- matrix_inverse(A)
  pander(inv_c)
  pander("A %*% A^-1"); pander(A %*% inv_c)
  
  ## -------------- Doolittle (no pivot) --------------
  section("Doolittle LU Decomposition (no pivot)")
  dool_ok <- TRUE
  lu_d <- NULL
  tryCatch({
    lu_d <- doolittle_decomposition(A)
    pander("L (Doolittle)"); pander(lu_d$L)
    pander("U (Doolittle)"); pander(lu_d$U)
  }, error = function(e) {
    dool_ok <<- FALSE
    pander("Non-pivoted Doolittle failed:")
    pander(as.character(e))
  })
  
  if (!is.null(b) && dool_ok) {
    section("Solve Ax = b (Doolittle, no pivot)")
    x_d <- solve_doolittle(A, b)
    pander("x"); pander(x_d)
    pander("A %*% x"); pander(A %*% x_d)
  }
  
  if (dool_ok) {
    section("det(A) via Doolittle (no pivot)")
    pander(det_doolittle(A))
    
    section("Inverse(A) via Doolittle (no pivot)")
    inv_d <- inverse_doolittle(A)
    pander(inv_d)
    pander("A %*% A^-1"); pander(A %*% inv_d)
  }
  
  ## -------------- Crout with Pivoting ---------------
  section("Crout LU with Partial Pivoting (PA = L U)")
  lu_cp <- crout_decomposition_pivoting(A)
  pander("P (Crout pivot)"); pander(lu_cp$P)
  pander("L (Crout pivot)"); pander(lu_cp$L)
  pander("U (Crout pivot)"); pander(lu_cp$U)
  pander("P %*% A"); pander(lu_cp$P %*% A)
  pander("L %*% U"); pander(lu_cp$L %*% lu_cp$U)
  
  if (!is.null(b)) {
    section("Solve Ax = b (Crout Pivoted)")
    b2 <- lu_cp$P %*% b
    y_cp <- forward_substitution(lu_cp$L, b2)
    x_cp <- backward_substitution(lu_cp$U, y_cp)
    pander("x"); pander(x_cp)
    pander("A %*% x"); pander(A %*% x_cp)
  }
  
  ## --------- Doolittle with Partial Pivoting --------
  section("Doolittle LU with Partial Pivoting (PA = L U)")
  lu_dp <- doolittle_decomposition_pivoting(A)
  pander("P (Doolittle pivot)"); pander(lu_dp$P)
  pander("L (Doolittle pivot)"); pander(lu_dp$L)
  pander("U (Doolittle pivot)"); pander(lu_dp$U)
  pander("P %*% A"); pander(lu_dp$P %*% A)
  pander("L %*% U"); pander(lu_dp$L %*% lu_dp$U)
  
  if (!is.null(b)) {
    section("Solve Ax = b (Doolittle Pivoted)")
    x_dp <- solve_doolittle_pivot(A, b)
    pander("x"); pander(x_dp)
    pander("A %*% x"); pander(A %*% x_dp)
  }
  
  section("det(A) via Doolittle Pivoted")
  pander(det_doolittle_pivot(A))
  
  section("Inverse(A) via Doolittle Pivoted")
  inv_dp <- inverse_doolittle_pivot(A)
  pander(inv_dp)
  pander("A %*% A^-1"); pander(A %*% inv_dp)
  
  ## ----------------- REF / RREF / Rank ---------------
  section("REF and RREF")
  ref_res <- row_echelon_form(list(A = A), "A", form = "both")
  pander("REF"); pander(ref_res[[1]])
  pander("RREF"); pander(ref_res[[2]])
  pander("Pivot columns"); pander(ref_res[[3]])
  
  section("Rank of A")
  pander(rank_lu(A))
  
  section("End of Analysis")
}
