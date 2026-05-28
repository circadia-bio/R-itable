# Shared fixtures for all itable tests.
# Loaded automatically by testthat before any test file runs.
#
# Pedigree layout (make_test_pedigree):
#
#   Founders (sex):  1(M) 2(F)   3(M) 4(F)   5(M) 6(F)
#   F2:              7(M)=1x2    8(F)=3x4    9(M)=5x6   10(F)=5x6
#
#   NOTE: 9 and 10 are full siblings -- do NOT mate them.
#   F3 (study subjects 11-15):
#     11(M) = 7(M) x  8(F)   <- couple A x couple B offspring
#     12(F) = 7(M) x  8(F)   <- full sib of 11
#     13(M) = 9(M) x  8(F)   <- 9 is unrelated to 8's family via different founders
#     14(F) = 9(M) x  8(F)   <- full sib of 13
#     15(M) = 7(M) x 10(F)   <- unrelated to 11/12 via different mothers
#
#   Sex assignments:  8=F, 10=F are mothers; 7=M, 9=M are fathers. Consistent.
#   No consanguineous matings => all F3 diagonal = 1.

make_test_pedigree <- function() {
  data.frame(
    id  = 1:15,
    pat = c(0L,0L,0L,0L,0L,0L, 1L,3L,5L,5L, 7L,7L,9L,9L,7L),
    mom = c(0L,0L,0L,0L,0L,0L, 2L,4L,6L,6L, 8L,8L,8L,8L,10L),
    sex = c(1L,2L,1L,2L,1L,2L, 1L,2L,1L,2L, 1L,2L,1L,2L,1L)
  )
}

make_test_data <- function(seed = 42) {
  set.seed(seed)
  study_ids <- 11:15
  data.frame(
    IID     = study_ids,
    age     = round(runif(length(study_ids), 20, 70)),
    sex_num = c(1L, 2L, 1L, 2L, 1L),
    bmi     = round(rnorm(length(study_ids), 25, 4), 1),
    hdl     = round(rnorm(length(study_ids), 55, 12), 1)
  )
}

# Larger dataset for herit_vc / herit_batch tests.
# n_families nuclear families, each with 2 founders + 2 offspring.
# Phenotype simulated with shared family effect (h2 ~ 0.5).
make_family_data <- function(n_families = 50, seed = 123) {
  set.seed(seed)

  ped_list <- lapply(seq_len(n_families), function(f) {
    base <- (f - 1L) * 4L
    data.frame(
      id  = base + 1:4,
      pat = c(0L, 0L, base + 1L, base + 1L),
      mom = c(0L, 0L, base + 2L, base + 2L),
      sex = c(1L, 2L, sample(1:2, 2, replace = TRUE))
    )
  })
  ped_df <- do.call(rbind, ped_list)

  family_effect <- rnorm(n_families, 0, 1)
  pheno <- numeric(nrow(ped_df))
  for (f in seq_len(n_families)) {
    rows        <- ((f - 1L) * 4L + 1L):((f - 1L) * 4L + 4L)
    pheno[rows] <- family_effect[f] + rnorm(4L, 0, 1)
  }

  off_rows <- ped_df$pat != 0L
  data_df  <- data.frame(
    IID     = ped_df$id[off_rows],
    age     = round(runif(sum(off_rows), 20, 70)),
    sex_num = ped_df$sex[off_rows],
    trait1  = pheno[off_rows],
    trait2  = pheno[off_rows] + rnorm(sum(off_rows), 0, 0.5)
  )

  list(ped = ped_df, data = data_df, study_ids = data_df$IID)
}
