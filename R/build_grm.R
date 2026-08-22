#' Build an additive genetic relationship matrix from a pedigree
#'
#' Constructs the additive genetic relationship matrix **A** (= 2 x kinship
#' matrix) for a set of study subjects, using a pedigree that may include
#' additional founder parents not in the study sample.
#'
#' @param ped_df A data frame with at least the columns `id`, `pat`, `mom`,
#'   and `sex`. Missing parents should be `NA` or `0`. `sex` should be
#'   numeric: `1` = male, `2` = female (any other value is recoded to `1`
#'   with a warning).
#' @param study_ids Character or integer vector of IDs for the study subjects
#'   whose sub-matrix should be extracted. Defaults to all IDs in `ped_df`.
#' @param id_col,pat_col,mom_col,sex_col Column names for the four required
#'   pedigree fields. Defaults are `"id"`, `"pat"`, `"mom"`, `"sex"`.
#' @param mz_col Optional column name in `ped_df` identifying monozygotic
#'   twin/multiple groups (e.g. a shared family label present only for MZ
#'   sets, blank/`NA` for everyone else, as in a `MZTWIN` column). All
#'   pairwise combinations of individuals sharing the same non-blank,
#'   non-`NA` value of `mz_col` are passed to [kinship2::pedigree()] as
#'   monozygotic (`code = 1`), overriding the default full-sibling kinship
#'   of 0.25 with the correct MZ kinship of 0.5 (A = 1.0). Groups of size
#'   \eqn{\geq} 3 (e.g. MZ triplets) are supported: every pair within the
#'   group is marked. `NULL` (default) disables this and reproduces prior
#'   behaviour exactly.
#'
#' @return A symmetric numeric matrix of dimension
#'   `length(study_ids) x length(study_ids)`, with diagonal entries equal to
#'   1 for non-inbred individuals and row/column names matching `study_ids`.
#'
#' @details
#' The function calls [kinship2::kinship()] on the full pedigree (including
#' founders), then multiplies by 2 to obtain the additive relationship matrix,
#' and finally subsets to `study_ids`. Including founders in the pedigree
#' ensures that kinship coefficients between study subjects connected only
#' through founders are estimated correctly.
#'
#' For twin cohorts, `pat`/`mom` are typically synthetic per-pair founder
#' labels rather than genotyped individuals; such founders should still be
#' present as their own rows in `ped_df` (with `pat`/`mom` set to `NA`) so
#' that [kinship2::pedigree()] can resolve the family structure. Without an
#' `mz_col`, twins are treated as ordinary full siblings (A = 0.5); supplying
#' `mz_col` corrects MZ pairs to A = 1.0 while leaving DZ pairs (which
#' already have the correct full-sibling kinship) unchanged.
#'
#' @examples
#' # Minimal two-generation pedigree: two couples, four offspring
#' ped <- data.frame(
#'   id  = 1:8,
#'   pat = c(0, 0, 0, 0, 1, 1, 3, 3),
#'   mom = c(0, 0, 0, 0, 2, 2, 4, 4),
#'   sex = c(1, 2, 1, 2, 1, 2, 1, 2)
#' )
#' A <- build_grm(ped, study_ids = 5:8)
#' round(A, 3)
#'
#' # With an MZ-twin column: subjects 7 and 8 are MZ, not ordinary full sibs
#' ped$mztwin <- c(NA, NA, NA, NA, NA, NA, "MZ1", "MZ1")
#' A_mz <- build_grm(ped, study_ids = 5:8, mz_col = "mztwin")
#' round(A_mz, 3)  # A[7,8] = 1.0 instead of 0.5
#'
#' @importFrom kinship2 pedigree kinship
#' @importFrom utils head combn
#' @export
build_grm <- function(ped_df,
                      study_ids = NULL,
                      id_col    = "id",
                      pat_col   = "pat",
                      mom_col   = "mom",
                      sex_col   = "sex",
                      mz_col    = NULL) {

  # -- Input checks -----------------------------------------------------------
  required <- c(id_col, pat_col, mom_col, sex_col)
  missing  <- setdiff(required, names(ped_df))
  if (length(missing)) {
    rlang::abort(c(
      "Required pedigree columns not found:",
      paste0("  Missing: ", paste(missing, collapse = ", ")),
      "i" = "Set `id_col`, `pat_col`, `mom_col`, `sex_col` to match your data."
    ))
  }
  if (!is.null(mz_col) && !mz_col %in% names(ped_df)) {
    rlang::abort(paste0("`mz_col` column '", mz_col, "' not found in `ped_df`."))
  }

  ids <- ped_df[[id_col]]
  pat <- ped_df[[pat_col]]
  mom <- ped_df[[mom_col]]
  sex <- ped_df[[sex_col]]

  # Coerce blank-string missing parents to NA. kinship2::pedigree() accepts
  # NA for founders' parents regardless of whether `id` is numeric or
  # character (verified equivalent to the numeric-only `0` convention used
  # by earlier versions of this function), so NA is used uniformly here.
  # This matters for cohorts (e.g. twin studies) where `pat`/`mom` are
  # character founder labels rather than numeric IDs.
  pat[!is.na(pat) & pat == ""] <- NA
  mom[!is.na(mom) & mom == ""] <- NA

  # Coerce sex: must be 1 or 2
  sex_num <- suppressWarnings(as.integer(sex))
  bad_sex <- is.na(sex_num) | !(sex_num %in% 1:2)
  if (any(bad_sex)) {
    rlang::warn(c(
      paste0(sum(bad_sex),
             " individual(s) have unrecognised sex values - recoded to 1 (male)."),
      "i" = "Expected 1 (male) or 2 (female)."
    ))
    sex_num[bad_sex] <- 1L
  }

  # -- MZ-twin relation matrix --------------------------------------------------
  relation <- NULL
  if (!is.null(mz_col)) {
    mz_grp <- ped_df[[mz_col]]
    grp    <- split(ids, mz_grp)
    grp    <- grp[!is.na(names(grp)) & trimws(names(grp)) != ""]
    grp    <- Filter(function(g) length(g) >= 2, grp)
    if (length(grp)) {
      pair_mat <- do.call(rbind, lapply(grp, function(g) t(utils::combn(g, 2))))
      relation <- data.frame(id1 = pair_mat[, 1], id2 = pair_mat[, 2], code = 1)
    } else {
      rlang::warn(paste0(
        "`mz_col` = '", mz_col, "' had no groups of size >= 2; ",
        "no MZ relatedness was applied."
      ))
    }
  }

  # -- Build pedigree & kinship -----------------------------------------------
  ped  <- if (is.null(relation)) {
    kinship2::pedigree(id    = ids,
                       dadid = pat,
                       momid = mom,
                       sex   = sex_num)
  } else {
    kinship2::pedigree(id       = ids,
                       dadid    = pat,
                       momid    = mom,
                       sex      = sex_num,
                       relation = relation)
  }
  phi2 <- 2 * kinship2::kinship(ped)

  # -- Subset to study subjects -----------------------------------------------
  if (is.null(study_ids)) study_ids <- ids
  study_ids <- as.character(study_ids)

  not_found <- setdiff(study_ids, rownames(phi2))
  if (length(not_found)) {
    rlang::abort(c(
      paste0(length(not_found), " study_id(s) not found in the pedigree:"),
      paste0("  First few: ", paste(head(not_found, 5), collapse = ", "))
    ))
  }

  phi2[study_ids, study_ids]
}
