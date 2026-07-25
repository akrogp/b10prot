# Proteomics Score Processing Functions

LIMIT_PROB <- 1e-300
LIMIT_COLOG <- -log10(LIMIT_PROB)

#' Compute Logarithmic Transformation for Probabilities
#'
#' This function computes a custom logarithmic transformation of probabilities.
#' If the probability is less than `LIMIT_PROB`, it returns `LIMIT_COLOG`.
#' Otherwise, it returns the negative base-10 logarithm of the probability.
#'
#' @param prob A numeric vector of probabilities.
#'
#' @return A numeric vector where values less than `LIMIT_PROB` are replaced
#' by `LIMIT_COLOG`, and other values are transformed using `-log10()`.
#'
#' @examples
#' # Example with a single probability value
#' colog(0.001)
#'
#' # Example with a vector of probabilities
#' colog(c(0.001, 0.05, 0.1))
#'
#' @export
colog <- function(prob) {
  ifelse(prob < LIMIT_PROB, LIMIT_COLOG, -log10(prob))
}

arrange_score <- function(data, score, lower_better = TRUE) {
  if( lower_better )
    data %>% arrange({{score}})
  else
    data %>% arrange(desc({{score}}))
}

diff_score <- function(score1, score2, lower_better = TRUE) {
  if(lower_better)
    score1 - score2
  else
    score2 - score1
}

#' Target-Decoy Approach for FDR Estimation
#'
#' This function applies the Target-Decoy Approach (TDA) to estimate the
#' False Discovery Rate (FDR) based on a scoring metric for the potential
#' identifications. It calculates the p-value, local confidence score (LP),
#' and q-value (FDR) for each identification.
#'
#' @param data A data frame containing the identification data. It must
#' include a logical column `isDecoy` indicating whether each row corresponds
#' to a decoy identification.
#' @param score The column name of the score used to rank identifications.
#' This should be an unquoted column name.
#' @param lower_better A logical value indicating whether lower scores are
#' better (default is `TRUE`).
#' @param unbiased When set to `TRUE`, one is added to the decoy and target
#' counts before FDR calculation. This adjustment can reduce bias in small datasets.
#'
#' @return A data frame with the original data and additional columns:
#' \describe{
#'   \item{decoys}{The cumulative number of decoys up to each identification.}
#'   \item{targets}{The cumulative number of targets up to each identification.}
#'   \item{pval}{The p-value estimated using the target-decoy approach.}
#'   \item{LP}{The local confidence score, computed using the \code{colog()} function.}
#'   \item{FDR}{The false discovery rate (FDR) for each score threshold.}
#'   \item{qval}{The cumulative minimum FDR (q-value).}
#' }
#'
#' @examples
#' # Example usage with a sample dataset
#' sample_data <- data.frame(
#'   score = c(0.01, 0.02, 0.03, 0.04, 0.04, 0.05, 0.06, 0.07),
#'   isDecoy = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE)
#' )
#' target_decoy_approach(sample_data, score, lower_better = TRUE)
#'
#' @export
target_decoy_approach <- function(data, score, lower_better = TRUE, unbiased = FALSE ) {
  check_required_cols(data, c("isDecoy", as_name(enquo(score))))
  DECOYS <- sum(data$isDecoy)
  data %>%
    arrange_score({{score}}, lower_better) %>%
    mutate(decoys = cumsum(isDecoy), targets = cumsum(!isDecoy)) %>%
    mutate(decoys = if(unbiased) decoys + 1 else decoys) %>%
    mutate(target = if(unbiased) targets + 1 else targets) %>%
    group_by({{score}}) %>%
    mutate(decoys = max(decoys), targets = max(targets)) %>%
    ungroup() %>%
    mutate(pval = (decoys + ifelse(isDecoy, -0.5, 0.5)) / DECOYS) %>%
    mutate(LP = colog(pval)) %>%
    mutate(FDR = decoys/targets) %>%
    arrange_score({{score}}, !lower_better) %>%
    mutate(qval = cummin(FDR)) %>%
    arrange_score({{score}}, lower_better)
}

#' Compute Global False Discovery Rate (FDR)
#'
#' This function calculates the global False Discovery Rate (FDR) as a percentage
#' based on the total number of target and decoy identifications in the data.
#'
#' @param data A data frame containing the identification data. It must include
#' a logical column `isDecoy` indicating whether each row corresponds to a decoy identification.
#'
#' @return A data frame with the total counts of target and decoy identifications,
#' and the global FDR expressed as a percentage.
#'
#' @examples
#' # Example usage with a sample dataset
#' sample_data <- data.frame(
#'   isDecoy = c(TRUE, FALSE, FALSE, TRUE, FALSE)
#' )
#' global_fdr(sample_data)
#'
#' @export
global_fdr <- function(data) {
  data %>%
    check_required_cols(c("isDecoy")) %>%
    group_by(isDecoy) %>%
    summarise(count = n()) %>%
    pivot_wider(names_from = isDecoy, values_from = count) %>%
    rename(Target=1, Decoy=2) %>%
    mutate(`Global FDR (%)` = Decoy/Target*100)
}

#' Refined False Discovery Rate (FDR) Calculation
#'
#' This function computes refined False Discovery Rate (FDR) estimates using a
#' competitive approach between target and decoy identifications. It provides
#' three types of FDR calculations: FDRn, FDRp, and FDRr, which adjust for
#' different competitive scenarios between targets and decoys.
#'
#' @param data A data frame containing the identification data, including
#' columns for the reference level, a score, and whether each identification
#' is a decoy (`isDecoy`).
#' @param levelRef The column name containing the reference level for each
#' identification (e.g., protein or gene reference). This should be an unquoted
#' column name.
#' @param score The column name of the score used to rank the identifications.
#' This should be an unquoted column name.
#' @param lower_better A logical value indicating whether lower scores are
#' better (default is `TRUE`).
#' @param affix A string indicating the suffix of prefix used to identify decoy
#' entries in the reference level column. Default is `"_REVERSED"`.
#'
#' @return A data frame with the original data and additional columns for the
#' refined FDR estimates:
#' \describe{
#'   \item{FDRn}{Normal FDR estimation as cumulative minimum (q-value).}
#'   \item{FDRp}{Picked FDR estimation as cumulative minimum (q-value).}
#'   \item{FDRr}{Refined FDR estimationas cumulative minimum (q-value).}
#'   \item{to}{Target-only identifications count.}
#'   \item{do}{Decoy-only identifications count.}
#'   \item{td}{Count of identifications with the same target and decoy scores.}
#'   \item{tb}{Target-best identifications count.}
#'   \item{db}{Decoy-best identifications count.}
#' }
#'
#' @examples
#' # Example usage with a sample dataset
#' sample_data <- data.frame(
#'   proteinRef = c("P1", "P1_REVERSED", "P2", "P3", "P3_REVERSED"),
#'   score = c(0.1, 0.2, 0.3, 0.5, 0.4),
#'   isDecoy = c(FALSE, TRUE, FALSE, FALSE, TRUE)
#' )
#' refined_fdr(sample_data, levelRef = proteinRef, score = score, lower_better = TRUE)
#'
#' @seealso
#' - [Protein Probability Model for High-Throughput Protein Identification by Mass Spectrometry-Based Proteomics](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00819) for more information on the refined FDR estimation.
#'
#' @export
refined_fdr <- function(data, levelRef, score, lower_better = TRUE, affix = "_REVERSED") {
  competition <-
    data %>%
    check_required_cols(c("isDecoy", as_name(enquo(levelRef)), as_name(enquo(score)))) %>%
    mutate(competitionRef = str_remove({{levelRef}}, affix)) %>%
    select(competitionRef, isDecoy, {{score}}) %>%
    mutate(isDecoy = ifelse(isDecoy, "decoy", "target")) %>%
    pivot_wider(names_from = isDecoy, values_from = {{score}}) %>%
    mutate(decoy = (if("decoy" %in% names(.)) decoy else NA)) %>%
    mutate(region = case_when(
      is.na(decoy) ~ "to",
      is.na(target) ~ "do",
      target == decoy ~ "td",
      diff_score(target, decoy, lower_better) < 0 ~ "tb",
      TRUE ~ "db"
    ))

  data %>%
    rowwise() %>%
    mutate(
      to = sum(
        diff_score(competition$target, {{score}}, lower_better) <= 0 &
          (
            competition$region == "to" |
              (
                competition$region == "tb" &
                  diff_score(competition$decoy, {{score}}, lower_better) > 0
              )
          ),
        na.rm = TRUE),
      do = sum(
        diff_score(competition$decoy, {{score}}, lower_better) <= 0 &
          (
            competition$region == "do" |
              (
                competition$region == "db" &
                  diff_score(competition$target, {{score}}, lower_better) > 0
              )
          ),
        na.rm = TRUE),
      td = sum(
        diff_score(competition$target, {{score}}, lower_better) <= 0 &
          competition$region == "td",
        na.rm = TRUE),
      tb = sum(
        diff_score(competition$target, {{score}}, lower_better) <= 0 &
          diff_score(competition$decoy, {{score}}, lower_better) <= 0 &
          competition$region == "tb",
        na.rm = TRUE),
      db = sum(
        diff_score(competition$target, {{score}}, lower_better) <= 0 &
          diff_score(competition$decoy, {{score}}, lower_better) <= 0 &
          competition$region == "db",
        na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      FDRn = (do + db + tb + td) / (to + db + tb + td),
      FDRp = (do + db + td) / (to + tb + td),
      FDRr = (do + 2*db + td) / (db + tb + to + td)
    ) %>%
    arrange_score({{score}}, !lower_better) %>%
    mutate(across(c(FDRn, FDRp, FDRr), cummin)) %>%
    arrange_score({{score}}, lower_better)
}

#' Picked Protein Group FDR
#'
#' Implements the picked protein group false discovery rate (FDR) strategy as
#' described by The et al. (2022).
#'
#' @param data A tibble containing protein group information. Must include at
#' least the columns groupRef, proteinMaster, proteinRefs, isDecoy and the scoring
#' column provided in score.
#' @param score The column name of the score used to rank the identifications.
#'   This should be an unquoted column name.
#' @param lower_better Logical; if TRUE, lower scores indicate better identifications (default `TRUE`).
#' @param affix String indicating the suffix/prefix used to identify decoy entries (default `"_REVERSED"`).
#'
#' @return A tibble containing the subset of protein groups that passed the
#' picked protein group filtering and the subsequent target-decoy FDR analysis.
#' The output will typically include the original columns from data and any
#' additional FDR-related metrics added by target_decoy_approach().
#'
#' @details
#' Protein groups are first sorted by a scoring column, and for each leading protein
#' any subsequent group containing it within its members is removed.
#'
#' Internally, the function:
#' \enumerate{
#' \item Normalizes protein identifiers by removing the affix (decoy tag).
#' \item Orders protein groups according to the chosen score and direction.
#' \item Compares each group against subsequent groups to do the competition.
#' \item Retains only the winner groups.
#' \item Applies a target-decoy approach to estimate FDR.
#' }
#'
#' @references
#' The M, Samaras P, Kuster B, Wilhelm M.
#' Reanalysis of ProteomicsDB Using an Accurate, Sensitive, and Scalable False
#' Discovery Rate Estimation Approach for Protein Groups.
#' Mol Cell Proteomics. 2022 Dec;21(12):100437. \doi{10.1016/j.mcpro.2022.100437}
#'
#' @examples
#' library(dplyr)
#' library(stringr)
#'
#' df <- tibble(
#'   groupRef = c("a", "b", "c", "d", "e", "f_REVERSED", "f"),
#'   proteinMaster = c("a", "b", "c", "d", "e", "f_REVERSED", "f"),
#'   proteinRefs = c(
#'     "a;x",
#'     "a_REVERSED;b_REVERSED",
#'     "b_REVERSED;c_REVERSED;x_REVERSED",
#'     "d",
#'     "e;d",
#'     "f_REVERSED",
#'     "f"
#'   ),
#'   score = 7:1
#' ) %>%
#'   mutate(isDecoy = str_detect(proteinRefs, "_REVERSED"))
#'
#' result <- b10prot:::picked_gfdr(df, score, lower_better = FALSE)
#' glimpse(result)
#'
#' @seealso [target_decoy_approach()]
#'
#' @keywords internal
picked_gfdr <- function(data, score, lower_better = TRUE, affix = "_REVERSED") {
  data %>%
    check_required_cols(c("groupRef", "proteinMaster", "proteinRefs", "isDecoy", as_name(enquo(score)))) %>%
    select(groupRef, proteinMaster, proteinRefs, {{score}}) %>%
    mutate(proteinMaster = str_remove(proteinMaster, affix)) %>%
    mutate(proteinRefs = str_remove(proteinRefs, affix)) %>%
    arrange_score({{score}}, lower_better) %>%
    mutate(idx = row_number()) %>%
    cross_join(., .) %>%
    filter(idx.y >= idx.x) %>%
    mutate(proteinRefs.y = str_split(proteinRefs.y, ";")) %>%
    rowwise() %>%
    mutate(keep.y = !(proteinMaster.x %in% proteinRefs.y)) %>%
    ungroup() %>%
    mutate(keep.y = ifelse(groupRef.x == groupRef.y, TRUE, keep.y)) %>%
    group_by(groupRef.y) %>%
    summarise(keep = all(keep.y)) %>%
    rename(groupRef = groupRef.y) %>%
    filter(keep) %>%
    select(-keep) %>%
    inner_join(data, by = join_by(groupRef)) %>%
    target_decoy_approach({{score}}, lower_better)
}

#' Refined Group-Level False Discovery Rate (FDRr) Calculation
#'
#' This function computes refined False Discovery Rate (FDR) estimates at the protein
#' group level using a competitive approach. It extends `refined_fdr` by handling
#' multiple proteins per identification in the target-decoy competition.
#'
#' @param data A data frame containing identification data. It must include:
#'   \itemize{
#'     \item `proteinMaster`: the main protein identifier for the group.
#'     \item `proteinRefs`: all proteins associated with the identification, separated by ";".
#'     \item `isDecoy`: logical indicating whether the identification is a decoy.
#'     \item A score column used for ranking identifications.
#'   }
#' @param score The column name of the score used to rank the identifications.
#'   This should be an unquoted column name.
#' @param lower_better Logical; if TRUE, lower scores indicate better identifications (default `TRUE`).
#' @param affix String indicating the suffix/prefix used to identify decoy entries (default `"_REVERSED"`).
#'
#' @return A data frame containing the original data along with additional columns:
#' \describe{
#'   \item{FDRn}{Normal FDR estimation as cumulative minimum (q-value).}
#'   \item{FDRp}{Picked FDR estimation as cumulative minimum (q-value).}
#'   \item{FDRr}{Refined FDR estimation as cumulative minimum (q-value).}
#'   \item{to}{Target-only identifications count at the protein group level.}
#'   \item{do}{Decoy-only identifications count at the protein group level.}
#'   \item{td}{Count of identifications with equal target and decoy scores.}
#'   \item{tb}{Target-best identifications count at the protein group level.}
#'   \item{db}{Decoy-best identifications count at the protein group level.}
#' }
#'
#' @details
#' The function first expands multiple protein references per identification into
#' individual rows, removes decoy affixes, and then determines which protein references
#' are in competition with the proteinMaster. It calculates the FDR estimates by
#' summing target and decoy identifications across regions (to, do, td, tb, db),
#' following the competitive approach described in the cited paper.
#'
#' @seealso
#' - [Refined FDR for Single Proteins](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00819)
#' - [Competitive Protein Group FDR](https://pubmed.ncbi.nlm.nih.gov/36328188)
#'
#' @export
refined_gfdr <- function(data, score, lower_better = TRUE, affix = "_REVERSED") {
  competition <-
    data %>%
    check_required_cols(c("isDecoy", "proteinMaster", "proteinRefs", as_name(enquo(score)))) %>%
    separate_rows(proteinRefs, sep = ";") %>%
    mutate(competitionRef = str_remove(proteinMaster, affix)) %>%
    mutate(proteinRef = str_remove(proteinRefs, affix)) %>%
    group_by(proteinRef) %>%
    mutate(pair = proteinRef %in% competitionRef) %>%
    ungroup() %>%
    filter(pair) %>%
    select(proteinRef, isDecoy, {{score}}) %>%
    mutate(isDecoy = ifelse(isDecoy, "decoy", "target")) %>%
    pivot_wider(names_from = isDecoy, values_from = {{score}}) %>%
    mutate(decoy = (if("decoy" %in% names(.)) decoy else NA)) %>%
    mutate(region = case_when(
      is.na(decoy) ~ "to",
      is.na(target) ~ "do",
      target == decoy ~ "td",
      diff_score(target, decoy, lower_better) < 0 ~ "tb",
      TRUE ~ "db"
    ))

  data %>%
    rowwise() %>%
    mutate(
      to = sum(
        diff_score(competition$target, {{score}}, lower_better) <= 0 &
          (
            competition$region == "to" |
              (
                competition$region == "tb" &
                  diff_score(competition$decoy, {{score}}, lower_better) > 0
              )
          ),
        na.rm = TRUE),
      do = sum(
        diff_score(competition$decoy, {{score}}, lower_better) <= 0 &
          (
            competition$region == "do" |
              (
                competition$region == "db" &
                  diff_score(competition$target, {{score}}, lower_better) > 0
              )
          ),
        na.rm = TRUE),
      td = sum(
        diff_score(competition$target, {{score}}, lower_better) <= 0 &
          competition$region == "td",
        na.rm = TRUE),
      tb = sum(
        diff_score(competition$target, {{score}}, lower_better) <= 0 &
          diff_score(competition$decoy, {{score}}, lower_better) <= 0 &
          competition$region == "tb",
        na.rm = TRUE),
      db = sum(
        diff_score(competition$target, {{score}}, lower_better) <= 0 &
          diff_score(competition$decoy, {{score}}, lower_better) <= 0 &
          competition$region == "db",
        na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      FDRn = (do + db + tb + td) / (to + db + tb + td),
      FDRp = (do + db + td) / (to + tb + td),
      FDRr = (do + 2*db + td) / (db + tb + to + td)
    ) %>%
    arrange_score({{score}}, !lower_better) %>%
    mutate(across(c(FDRn, FDRp, FDRr), cummin)) %>%
    arrange_score({{score}}, lower_better)
}

#' LP Gamma (LPG) Metrics Calculation
#'
#' This function calculates various LP Gamma (LPG) metrics for a specified
#' reference level, such as a protein, a gene or a protein group. The metrics are
#' derived from the coLogarithm of Probability (LP) of their constituting peptides
#' and include probabilities derived from maximum LP, sum of LP values, and
#' filtered LP sums using a specified peptide-level FDR threshold.
#'
#' @param data A data frame containing identification data, including
#' peptide-related columns for coLogarithm of Probability (`LP`), q-values
#' (`qval`), and a logical column `isDecoy`.
#' @param levelRef The column name of the reference level to group by, such as
#' a protein or a gene identifier. This should be an unquoted column name.
#' @param threshold A numeric value representing the FDR threshold for
#' peptide-level q-values (default is `0.01`).
#' @param extra_cols A character vector of optional columns to include in the
#'   summarisation step (if present in the input data). By default, it attempts
#'   to include `proteinType`, `proteinCount`, `proteinRefs`, and `proteinMaster`.
#'
#' @return A data frame containing the calculated LPG metrics with the following columns:
#' \describe{
#'   \item{isDecoy}{Indicates whether the group contains any decoy identification.}
#'   \item{n}{The total number of peptide identifications for the group.}
#'   \item{m}{The number of peptide identifications with a q-value below the threshold.}
#'   \item{LPM}{The maximum coLogarithm of Probability (`LP`) for the group.}
#'   \item{LPS}{The sum of coLogarithm of Probability (`LP`) for the group.}
#'   \item{LPF}{The sum of coLogarithm of Probability for identifications with a q-value below the threshold.}
#'   \item{LPGM}{The LP Gamma value based on the maximum `LP`.}
#'   \item{LPGS}{The LP Gamma value based on the sum of `LP` values.}
#'   \item{LPGF}{The LP Gamma value based on the filtered sum of `LP` values for confident identifications.}
#' }
#'
#' @examples
#' # Example usage with a sample dataset
#' sample_data <- data.frame(
#'   levelRef = c("P1", "P1", "P2", "P2", "P3"),
#'   LP = c(1.5, 2.0, 0.5, 1.0, 1.2),
#'   qval = c(0.01, 0.02, 0.005, 0.03, 0.01),
#'   isDecoy = c(FALSE, FALSE, TRUE, FALSE, FALSE)
#' )
#' lpg(sample_data, levelRef, threshold = 0.01)
#'
#' @seealso
#' - [Protein Probability Model for High-Throughput Protein Identification by Mass Spectrometry-Based Proteomics](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00819) for more information on the LPG scores.
#'
#' @export
lpg <- function(data, levelRef, threshold = 0.01, extra_cols = c()) {
  extra_cols <- c(extra_cols, "proteinType", "proteinCount", "proteinRefs", "proteinMaster")
  data %>%
    check_required_cols(c("isDecoy", "LP", "qval", as_name(enquo(levelRef)))) %>%
    group_by({{levelRef}}) %>%
    summarise(
      isDecoy = any(isDecoy),
      n = n(),
      m = sum(ifelse(qval <= threshold, 1, 0)),
      LPM = max(LP),
      LPS = sum(LP),
      LPF = sum(ifelse(qval <= threshold, LP, 0)),
      LPGM = colog(1 - (1 - 10^(-LPM))^n),
      LPGS = colog(1 - pgamma(LPS*log(10),n)),
      LPGF = ifelse(m == 0, LPGM, colog((1 - pgamma(LPF*log(10),m)) * choose(n, m))),
      across(any_of(extra_cols), first),
      .groups = "drop"
    )
}

#' Plot Rank of Decoy Scores
#'
#' This function creates a rank plot of decoy scores based on various LP (coLogarithm
#' of Probability) metrics, including LPM, LPS, LPF, and LPG scores (LPGM, LPGS, LPGF).
#' The plot displays the ranked scores of decoys with a reference line for comparison.
#'
#' @param data A data frame containing identification data, including columns
#' for decoy status (`isDecoy`) and any of the different LP metrics (LPM, LPS,
#' LPF, LPGM, LPGS, LPGF).
#'
#' @return A ggplot object showing the rank plot of decoy scores for the different
#' metrics. Each score type is displayed in a separate facet with the rank plotted
#' on the x-axis and the score on the y-axis. The red diagonal line represents a
#' reference for ideal ranking.
#'
#' @export
plot_rank <- function(data) {
  decoys <-
    data %>%
    filter(isDecoy)
  N <- nrow(decoys)
  decoys %>%
    pivot_longer(matches("^LP."), names_to = "score", values_to = "value") %>%
    group_by(score) %>%
    mutate(rank = rank(-value)) %>%
    ungroup() %>%
    mutate(score = factor(score, levels = c("LPM", "LPS", "LPF", "LPGM", "LPGS", "LPGF"))) %>%
    ggplot(aes(x=colog(rank/N), y=value)) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_point(color = "blue") +
    facet_wrap(vars(score), ncol = 3, scales = "free_y") +
    ggtitle("Decoy scores distribution")
}
