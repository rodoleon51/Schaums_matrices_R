############################################################
# test_matrix_tools.R
# Hybrid test suite for matrix_tools.R
# - Human-readable sections with pander
# - Automated checks with stopifnot / all.equal
############################################################
# Load required packages
library(pander)
panderOptions("digits", 6)

# Load toolbox
source("functions/matrix_tools.R")

# Helper functions
section <- function(title) {
  cat("\n\n")
  pander(paste0("## ", title))
  cat("\n")
}

ae <- function(x, y, tol = 1e-8) {
  isTRUE(all.equal(x, y, tolerance = tol))
}

pander("# Test Suite for matrix_tools.R")


############################################################
# TEST MATRICES
############################################################

section("Test Matrices")

# Define vectors BEFORE printing (to avoid pander interference)
b4   <- c(3, -1, 7, 13)
b3   <- c(1,  2, 3)
b_piv<- c(2,  3, -6)
b_tri<- c(1,  2, 3)

A4 <- matrix(c(
  1,  2, -2,  3,
  -1,  1,  0,  2,
  3, -3,  4,  1,
  2,  1,  1, -2
), 4, 4, byrow = TRUE)

A3 <- matrix(c(
  4,  1, -2,
  1,  2,  3,
  -1,  3,  1
), 3, 3, byrow = TRUE)

A_piv <- matrix(c(
  0,  1,  1,
  5,  1, -1,
  2, -3, -3
), 3, 3, byrow = TRUE)

U3 <- matrix(c(
  2,  1,  3,
  0,  1,  2,
  0,  0,  3
), 3, 3, byrow = TRUE)

L3 <- matrix(c(
  4,  0,  0,
  1,  2,  0,
  -1,  3,  1
), 3, 3, byrow = TRUE)

pander("### A4")
pander(A4)
pander("### b4")
pander(b4)

pander("### A3")
pander(A3)
pander("### b3")
pander(b3)

pander("### A_piv (needs pivoting)")
pander(A_piv)
pander("### b_piv")
pander(b_piv)

pander("### U3 (upper triangular)")
pander(U3)
pander("### L3 (lower triangular)")
pander(L3)


############################################################
# 1. Crout LU (no pivot)
############################################################

section("Crout LU Decomposition (no pivot) on A4")

lu_c <- crout_decomposition(A4)
L_c  <- lu_c$L
U_c  <- lu_c$U

pander("L (Crout)")
pander(L_c)
pander("U (Crout)")
pander(U_c)

# Check reconstruction
section("Check: Crout L * U ≈ A4")
LUc <- L_c %*% U_c
pander(LUc)
stopifnot(ae(as.numeric(LUc), as.numeric(A4)))

# Solve Ax=b
section("Solve Ax = b (Crout, no pivot)")
x_c <- solve_lu(A4, b4)
pander("x_c")
pander(x_c)

Ax_c <- A4 %*% x_c
pander("A4 %*% x_c")
pander(Ax_c)

# **FIXED**: compare numeric to numeric
stopifnot(ae(as.numeric(Ax_c), as.numeric(b4)))

# Determinant
section("det(A4) via Crout vs base::det")
det_c <- det_lu(A4)
pander(det_c)
stopifnot(ae(det_c, det(A4)))
if (det_c < 1e-16) {
  cat(round(det_c,0))
  cat("\nThe matrix is singular = non invertible.")
}

# Inverse
section("Inverse(A4) via Crout vs base::solve")
inv_c <- matrix_inverse(A4)
pander(inv_c)
stopifnot(ae(as.numeric(inv_c), as.numeric(solve(A4))))
stopifnot(ae(as.numeric(A4 %*% inv_c), as.numeric(diag(4))))


############################################################
# 2. Doolittle LU (no pivot) on A3
############################################################

section("Doolittle LU Decomposition (no pivot) on A3")

lu_d <- doolittle_decomposition(A3)
L_d  <- lu_d$L
U_d  <- lu_d$U

pander("L (Doolittle)")
pander(L_d)
pander("U (Doolittle)")
pander(U_d)

# Check reconstruction
section("Check: Doolittle L * U ≈ A3")
LUd <- L_d %*% U_d
pander(LUd)
stopifnot(ae(as.numeric(LUd), as.numeric(A3)))

# Solve Ax=b
section("Solve Ax = b (Doolittle, no pivot)")
x_d <- solve_doolittle(A3, b3)
pander("x_d")
pander(x_d)

Ax_d <- A3 %*% x_d
stopifnot(ae(as.numeric(Ax_d), as.numeric(b3)))

# Determinant
section("det(A3) via Doolittle vs base::det")
det_d <- det_doolittle(A3)
pander(det_d)
stopifnot(ae(det_d, det(A3)))
if (det_d < 1e-16) {
  cat(round(det_d,0))
  cat("\nThe matrix is singular = non invertible.")
  }

# Inverse
section("Inverse(A3) via Doolittle vs base::solve")
inv_d <- inverse_doolittle(A3)
pander(inv_d)
stopifnot(ae(as.numeric(inv_d), as.numeric(solve(A3))))
stopifnot(ae(as.numeric(A3 %*% inv_d), as.numeric(diag(3))))


############################################################
# 3. Crout with Partial Pivoting on A_piv
############################################################

section("Crout LU with Partial Pivoting on A_piv")

lu_cp <- crout_decomposition_pivoting(A_piv)
P_cp  <- lu_cp$P
L_cp  <- lu_cp$L
U_cp  <- lu_cp$U

pander("P (Crout pivot)")
pander(P_cp)
pander("L (Crout pivot)")
pander(L_cp)
pander("U (Crout pivot)")
pander(U_cp)

section("Check: P %*% A_piv ≈ L * U (Crout pivoted)")
PA_cp <- P_cp %*% A_piv
LU_cp <- L_cp %*% U_cp
pander("P_cp %*% A_piv")
pander(PA_cp)
pander("L_cp %*% U_cp")
pander(LU_cp)

stopifnot(ae(as.numeric(PA_cp), as.numeric(LU_cp)))

# Solve Ax = b with pivoted Crout
section("Solve A_piv x = b_piv (Crout pivoted)")
b2_cp <- P_cp %*% b_piv
y_cp  <- forward_substitution(L_cp, b2_cp)
x_cp  <- backward_substitution(U_cp, y_cp)
pander("x_cp")
pander(x_cp)
stopifnot(ae(as.numeric(A_piv %*% x_cp), as.numeric(b_piv)))


############################################################
# 4. Doolittle with Partial Pivoting on A_piv
############################################################

section("Doolittle LU with Partial Pivoting on A_piv")

lu_dp <- doolittle_decomposition_pivoting(A_piv)
P_dp  <- lu_dp$P
L_dp  <- lu_dp$L
U_dp  <- lu_dp$U

pander("P (Doolittle pivot)")
pander(P_dp)
pander("L (Doolittle pivot)")
pander(L_dp)
pander("U (Doolittle pivot)")
pander(U_dp)

section("Check: P %*% A_piv ≈ L * U (Doolittle pivot)")
PA_dp <- P_dp %*% A_piv
LU_dp <- L_dp %*% U_dp
pander("P_dp %*% A_piv")
pander(PA_dp)
pander("L_dp %*% U_dp")
pander(LU_dp)
stopifnot(ae(as.numeric(PA_dp), as.numeric(LU_dp)))

# Solve Ax = b with pivoted Doolittle
section("Solve A_piv x = b_piv (Doolittle pivoted)")
x_dp <- solve_doolittle_pivot(A_piv, b_piv)
pander("x_dp")
pander(x_dp)
stopifnot(ae(as.numeric(A_piv %*% x_dp), as.numeric(b_piv)))

# Determinant (pivoted Doolittle) vs base::det
section("det(A_piv) via Doolittle pivoted vs base::det")
det_dp <- det_doolittle_pivot(A_piv)
pander(det_dp)
stopifnot(ae(det_dp, det(A_piv)))

# Inverse via pivoted Doolittle vs solve()
section("Inverse(A_piv) via Doolittle pivoted vs base::solve")
inv_dp <- inverse_doolittle_pivot(A_piv)
pander(inv_dp)
stopifnot(ae(as.numeric(inv_dp), as.numeric(solve(A_piv))))
stopifnot(ae(as.numeric(A_piv %*% inv_dp), as.numeric(diag(3))))


############################################################
# 5. Triangular Matrix Handling
############################################################

section("Triangular Matrix Handling")

pander("U3 (upper triangular)")
pander(U3)
pander("L3 (lower triangular)")
pander(L3)

## Upper triangular
section("Upper triangular: solve vs base::solve")
xU_crout  <- solve_lu(U3, b_tri)
xU_dool   <- solve_doolittle(U3, b_tri)
xU_dool_p <- solve_doolittle_pivot(U3, b_tri)
xU_base   <- as.vector(solve(U3, b_tri))

stopifnot(ae(as.numeric(xU_crout),  xU_base))
stopifnot(ae(as.numeric(xU_dool),   xU_base))
stopifnot(ae(as.numeric(xU_dool_p), xU_base))

## Lower triangular
section("Lower triangular: solve vs base::solve")
xL_crout  <- solve_lu(L3, b_tri)
xL_dool   <- solve_doolittle(L3, b_tri)
xL_dool_p <- solve_doolittle_pivot(L3, b_tri)
xL_base   <- as.vector(solve(L3, b_tri))

stopifnot(ae(as.numeric(xL_crout),  xL_base))
stopifnot(ae(as.numeric(xL_dool),   xL_base))
stopifnot(ae(as.numeric(xL_dool_p), xL_base))


############################################################
# 6. REF, RREF, and Rank
############################################################

section("REF, RREF, and Rank")

A_rank <- matrix(c(
  1, 2, 3,
  2, 4, 6,
  1, 1, 1
), 3, 3, byrow=TRUE)

pander("A_rank")
pander(A_rank)

ref_res <- row_echelon_form(list(A=A_rank), "A", form="both")
REF  <- ref_res[[1]]
RREF <- ref_res[[2]]
piv  <- ref_res[[3]]

pander("REF(A_rank)")
pander(REF)
pander("RREF(A_rank)")
pander(RREF)
pander("Pivot columns")
pander(piv)

rank_our <- rank_lu(A_rank)
rank_qr  <- qr(A_rank)$rank

stopifnot(rank_our == rank_qr)


############################################################
# 7. run_full_analysis smoke test
############################################################

section("run_full_analysis() Smoke Tests")

run_full_analysis(A3,  b3,   name="A3, b3")
run_full_analysis(U3,  b_tri, name="U3, triangular")


############################################################
# END
############################################################

section("All tests completed successfully.")
pander("If no errors were raised, the toolbox is consistent with base R and internal checks.")
