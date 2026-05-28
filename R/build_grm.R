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
#' @importFrom kinship2 pedigree kinship
#' @importFrom utils head
#' @export
build_grm <- function(ped_df,
                      study_ids = NULL,
                      id_col    = "id",
                      pat_col   = "pat",
                      mom_col   = "mom",
                      sex_col   = "sex") {

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

  ids <- ped_df[[id_col]]
  pat <- ped_df[[pat_col]]
  mom <- ped_df[[mom_col]]
  sex <- ped_df[[sex_col]]

  # Coerce parents: NA / missing -> 0
  pat[is.na(pat)] <- 0L
  mom[is.na(mom)] <- 0L

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

  # -- Build pedigree & kinship -----------------------------------------------
  ped  <- kinship2::pedigree(id    = ids,
                             dadid = pat,
                             momid = mom,
                             sex   = sex_num)
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
