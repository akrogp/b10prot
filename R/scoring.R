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
    )
}

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
