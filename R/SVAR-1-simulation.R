

create_svar_data <- function(C, B, Y0, EE) {

  # A or B model ??
  # long-run or short run restrictions ??
  # Other properties ?? nonconstant ... so many things

  # --> keep it simple: short-run restrictions is what I need at the moment

  # provide structural coefficients: A_0, C, B (check notation in KilianLütkepohl..)


  # turn into reduced form
  A <- solve(A_0) %*% A_1
  EE <- solve(A_0) %*% E

  # normalise errors?
  # EE <- EE / sqrt(diag(A_0 %*% t(A_0)))

  create_varp_data(A, Y_0, EE)[, (1:n) + nburnin]
}
