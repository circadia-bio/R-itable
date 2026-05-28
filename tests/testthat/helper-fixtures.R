# Shared fixtures for all itable tests.
# Loaded automatically by testthat before any test file runs.
#
# make_test_pedigree() -- small pedigree for build_grm unit tests.
#
#   Four unrelated founder couples (F1):
#     Couple A: 1(M) x 2(F)
#     Couple B: 3(M) x 4(F)
#     Couple C: 5(M) x 6(F)
#     Couple D: 7(M) x 8(F)
#
#   F2 (one offspring per couple, all unrelated to each other):
#     9(M)  = 1 x 2   (from couple A)
#    10(F)  = 3 x 4   (from couple B)
#    11(M)  = 5 x 6   (from couple C)
#    12(F)  = 7 x 8   (from couple D)
#
#   F3 study subjects (IDs 13-16):
#    13(M) = 9  x 10  <- parents from A and B -- unrelated to 15/16
#    14(F) = 9  x 10  <- full sibling of 13
#    15(M) = 11 x 12  <- parents from C and D -- unrelated to 13/14
#    16(F) = 11 x 12  <- full sibling of 15
#
#   13 and 15 share NO founders => A["13","15"] = 0 exactly.
#   13 and 14 are full siblings  => A["13","14"] = 0.5 exactly.
#   All F3 diagonal = 1 (no inbreeding).

make_test_pedigree <- function() {
  data.frame(
    id  = 1:16,
    pat = c(0L,0L,0L,0L,0L,0L,0L,0L, 1L,3L,5L,7L, 9L,9L,11L,11L),
    mom = c(0L,0L,0L,0L,0L,0L,0L,0L, 2L,4L,6L,8L, 10L,10L,12L,12L),
    sex = c(1L,2L,1L,2L,1L,2L,1L,2L, 1L,2L,1L,2L,  1L,2L,1L,2L)
  )
}

make_test_data <- function(seed = 42) {
  set.seed(seed)
  study_ids <- 13:16
  data.frame(
    IID     = study_ids,
    age     = round(runif(length(study_ids), 20, 70)),
    sex_num = c(1L, 2L, 1L, 2L),
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
