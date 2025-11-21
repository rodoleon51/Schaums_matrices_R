# ================================================================
# matrix_tools_b.R
# Advanced Matrix Methods:
#   0. Utilities:
#        - round_small()
#        - print_matrix_clean()
#   1. Classical Gram-Schmidt (CGS)
#   2. Modified Gram-Schmidt (MGS)
#   3. Householder QR (Default)
#   4. Unified QR Decomposition Interface
#   5. QR Algorithm for Eigenvalues (unshifted)
#   6. Shifted QR Algorithm (Wilkinson shift)
# ================================================================


# ------------------------------------------------
# 0. Utility 1: Round very small values to zero
# ------------------------------------------------
round_small <- function(X, tol = 1e-12) {
  X[abs(X) < tol] <- 0
  X
}


# ------------------------------------------------
# 0. Utility 2: Pretty-print a matrix with rounding
# ------------------------------------------------
print_matrix_clean <- function(X, tol = 1e-12, digits = 6) {
  Xr <- round_small(X, tol)
  Xr <- round(Xr, digits = digits)
  print(Xr)
  invisible(Xr)
}



# ------------------------------------------------
# 1. Classical Gram–Schmidt
# ------------------------------------------------
gram_schmidt_cgs <- function(A, diag = FALSE) {
  A <- as.matrix(A)
  m <- nrow(A)
  n <- ncol(A)
  
  Q <- matrix(0, m, n)
  R <- matrix(0, n, n)
  
  if (diag) cat("[diag] Classical Gram-Schmidt\n")
  
  for (k in 1:n) {
    v <- A[, k]
    
    for (j in 1:(k - 1)) {
      R[j, k] <- t(Q[, j]) %*% A[, k]
      v <- v - R[j, k] * Q[, j]
    }
    
    R[k, k] <- sqrt(sum(v^2))
    if (R[k, k] < 1e-12) stop("Columns linearly dependent.")
    
    Q[, k] <- v / R[k, k]
  }
  
  list(Q = round_small(Q), R = round_small(R))
}



# ------------------------------------------------
# 2. Modified Gram–Schmidt (numerically more stable)
# ------------------------------------------------
gram_schmidt_mgs <- function(A, diag = FALSE) {
  A <- as.matrix(A)
  m <- nrow(A)
  n <- ncol(A)
  
  Q <- matrix(0, m, n)
  R <- matrix(0, n, n)
  
  if (diag) cat("[diag] Modified Gram-Schmidt\n")
  
  for (k in 1:n) {
    v <- A[, k]
    
    for (j in 1:(k - 1)) {
      R[j, k] <- t(Q[, j]) %*% v
      v <- v - R[j, k] * Q[, j]
    }
    
    R[k, k] <- sqrt(sum(v^2))
    if (R[k, k] < 1e-12) stop("Columns linearly dependent.")
    
    Q[, k] <- v / R[k, k]
  }
  
  list(Q = round_small(Q), R = round_small(R))
}



# ------------------------------------------------
# 3. Householder QR (Numerically Stable — Industry Standard)
# ------------------------------------------------
qr_householder <- function(A, diag = FALSE) {
  A <- as.matrix(A)
  m <- nrow(A)
  n <- ncol(A)
  
  R <- A
  Q <- diag(m)
  
  if (diag) cat("[diag] Householder QR\n")
  
  for (k in 1:min(m, n)) {
    
    x <- R[k:m, k]
    e1 <- c(1, rep(0, length(x) - 1))
    
    alpha <- -sign(x[1]) * sqrt(sum(x^2))
    v <- x - alpha * e1
    v <- v / sqrt(sum(v^2))
    
    Hk <- diag(length(v)) - 2 * (v %*% t(v))
    
    H <- diag(m)
    H[k:m, k:m] <- Hk
    
    R <- H %*% R
    Q <- Q %*% H
  }
  
  list(Q = round_small(Q), R = round_small(R))
}



# ------------------------------------------------
# 4. Unified QR Interface (Householder is default)
# ------------------------------------------------
qr_decomposition <- function(A,
                             method = c("householder", "mgs", "cgs"),
                             diag = FALSE) {
  method <- match.arg(method)
  
  if (method == "householder") {
    return(qr_householder(A, diag = diag))
  }
  if (method == "mgs") {
    return(gram_schmidt_mgs(A, diag = diag))
  }
  if (method == "cgs") {
    return(gram_schmidt_cgs(A, diag = diag))
  }
}



# ------------------------------------------------
# 5. QR Algorithm for Eigenvalues (unshifted)
# ------------------------------------------------
qr_algorithm <- function(A, max_iter = 1000, tol = 1e-10, diag = FALSE) {
  A <- as.matrix(A)
  n <- nrow(A)
  
  if (n != ncol(A)) stop("QR algorithm requires a square matrix.")
  
  Ak <- A
  
  if (diag) cat("[diag] Starting QR eigenvalue iteration (unshifted)...\n")
  
  for (i in 1:max_iter) {
    
    QR <- qr_decomposition(Ak, method = "householder")
    Q <- QR$Q
    R <- QR$R
    
    Ak_next <- R %*% Q
    
    off_norm <- sqrt(sum(Ak_next[upper.tri(Ak_next)]^2) +
                       sum(Ak_next[lower.tri(Ak_next)]^2))
    
    if (diag && i %% 10 == 0)
      cat(sprintf("[diag] Iter %d: off-diagonal norm = %.6e\n", i, off_norm))
    
    if (off_norm < tol) {
      if (diag) cat("[diag] Converged (unshifted QR).\n")
      return(list(
        values = round_small(diag(Ak_next)),
        matrix = round_small(Ak_next),
        iter   = i
      ))
    }
    
    Ak <- Ak_next
  }
  
  warning("QR algorithm (unshifted) did not converge.")
  
  list(
    values = round_small(diag(Ak)),
    matrix = round_small(Ak),
    iter   = max_iter
  )
}

# ------------------------------------------------
# Put this somewhere above qr_algorithm_shifted() in matrix_tools_b.R
# ------------------------------------------------

wilkinson_shift <- function(A) {
  n <- nrow(A)
  if (n < 2) return(A[n, n])

  # Bottom-right 2×2 block:
  # [a  b]
  # [c  d]
  a <- A[n-1, n-1]
  b <- A[n-1, n]
  c <- A[n,   n-1]
  d <- A[n,   n]

  trace <- a + d
  det2  <- a*d - b*c

  # Discriminant of λ^2 - trace λ + det2 = 0
  disc <- trace^2 - 4 * det2

  # Protect against small negative due to roundoff
  if (disc < 0) disc <- 0

  lambda1 <- (trace + sqrt(disc)) / 2
  lambda2 <- (trace - sqrt(disc)) / 2

  # Choose the eigenvalue closest to bottom-right entry d (Wilkinson rule)
  if (abs(lambda1 - d) < abs(lambda2 - d)) {
    lambda1
  } else {
    lambda2
  }
}

# ------------------------------------------------
# 6. Shifted QR Algorithm (Wilkinson shift)
# ------------------------------------------------
qr_algorithm_shifted <- function(A, max_iter = 1000, tol = 1e-10, diag = FALSE) {
  A <- as.matrix(A)
  n <- nrow(A)
  if (n != ncol(A)) stop("Shifted QR requires square matrix.")
  
  Ak <- A
  
  if (diag) cat("[diag] Starting SHIFTED QR eigenvalue iteration...\n")
  
  for (i in 1:max_iter) {
    
    if (n == 1) {
      if (diag) cat("[diag] n = 1, trivial case.\n")
      return(list(
        values = round_small(diag(Ak)),
        matrix = round_small(Ak),
        iter   = i
      ))
    }
    
    # Robust Wilkinson shift
    mu <- wilkinson_shift(Ak)
    
    # QR step on A - mu I
    QR <- qr_decomposition(Ak - mu * diag(n), method = "householder")
    Q  <- QR$Q
    R  <- QR$R
    
    Ak_next <- R %*% Q + mu * diag(n)
    
    off_norm <- sqrt(
      sum(Ak_next[upper.tri(Ak_next)]^2) +
        sum(Ak_next[lower.tri(Ak_next)]^2)
    )
    
    # Safety check: if something *still* goes wrong, fall back
    if (is.nan(off_norm) || is.infinite(off_norm)) {
      warning("Shifted QR produced NaN/Inf off-diagonal norm; falling back to unshifted QR.")
      return(qr_algorithm(A, max_iter = max_iter, tol = tol, diag = diag))
    }
    
    if (diag && i %% 10 == 0)
      cat(sprintf("[diag] Iter %d: off-norm = %.6e\n", i, off_norm))
    
    if (off_norm < tol) {
      if (diag) cat("[diag] Converged (shifted QR).\n")
      return(list(
        values = round_small(diag(Ak_next)),
        matrix = round_small(Ak_next),
        iter   = i
      ))
    }
    
    Ak <- Ak_next
  }
  
  warning("Shifted QR did not converge within max_iter.")
  list(
    values = round_small(diag(Ak)),
    matrix = round_small(Ak),
    iter   = max_iter
  )
}


# ------------------------------------------------
# Utility 3: Pretty-print using pander() + round_small()
# ------------------------------------------------
pander_clean <- function(X, tol = 1e-12, digits = 6) {
  Xc <- round_small(X, tol)
  Xc <- round(Xc, digits)
  pander::pander(Xc)
}
