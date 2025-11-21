# Define the matrix
E <- matrix(c(1, 2, 3,
              4, 9, 7),
            nrow = 2, byrow = TRUE)



row_echelon_no_swap <- function(mat, tol = 1e-10) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  for (i in 1:min(nr, nc)) {
    # Only swap rows if the current pivot is zero or smaller than tolerance
    if (abs(mat[i, i]) <= tol) {
      pivot_row <- i + which.max(abs(mat[i:nr, i])) - 1
      if (abs(mat[pivot_row, i]) > tol && pivot_row != i) {
        mat[c(i, pivot_row), ] <- mat[c(pivot_row, i), ]
      }
    }
    
    # Normalize the pivot row
    if (abs(mat[i, i]) > tol) {
      mat[i, ] <- mat[i, ] / mat[i, i]
    }
    
    # Eliminate entries below the pivot
    if (i < nr) {
      for (j in (i + 1):nr) {
        if (abs(mat[j, i]) > tol) {
          mat[j, ] <- mat[j, ] - mat[j, i] * mat[i, ]
        }
      }
    }
  }
  
  return(mat)
}

# Test the matrix E:
E <- matrix(c(1, 2, 3,
              4, 9, 7),
            nrow = 2, byrow = TRUE)

row_echelon_no_swap(E)
