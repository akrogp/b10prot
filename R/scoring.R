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
  data %>%
    arrange(case_when(
      lower_better ~ {{score}},
      TRUE ~ desc({{score}})
    ))
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
target_decoy_approach <- function(data, score, lower_better = TRUE ) {
  DECOYS <- sum(data$isDecoy)
  data %>%
    arrange_score({{score}}, lower_better) %>%
    mutate(decoys = cumsum(isDecoy), targets = cumsum(!isDecoy)) %>%
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
#'
#' @return A data frame containing the calculated LPG metrics:
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
lpg <- function(data, levelRef, threshold = 0.01) {
  data %>%
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
      across(everything(), first),
      .groups = "drop"
    ) %>%
    select(-LP, -qval)
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
