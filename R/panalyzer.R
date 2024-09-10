# Protein inference with PAnalyzer

add_class <- function(obj, new_class) {
  class(obj) <- c(new_class, class(obj))
  obj
}

#' Run PAnalyzer on Peptide-to-Protein Data
#'
#' Executes the PAnalyzer tool on a provided peptide-to-protein dataset. This
#' function writes the input data to a temporary file, runs the PAnalyzer tool
#' using Java, and reads the output into R.
#'
#' @param pep2prot A data frame containing peptide-to-protein mappings. The
#'   data frame must include the columns `peptideRef` and `proteinRef`.
#'
#' @return A data frame with the original data and additional columns:
#'   - `peptideType`: Type of the peptide determined by PAnalyzer. Possible values
#'     include `"unique"`, `"discriminating"`, and `"non-discriminating"`.
#'   - `proteinType`: Type of the protein determined by PAnalyzer. Possible values
#'     include `"conclusive"`, `"indistinguishable"`, `"ambiguous"`, and
#'     `"non-conclusive"`. All the proteins in a group will have the same
#'     `proteinType` value.
#'   - `groupRef`: Group reference from PAnalyzer.
#'
#' @seealso
#' - [PAnalyzer: A software tool for protein inference in shotgun proteomics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-288) for more information on the tool.
panalyzer <- function(pep2prot) {
  jar_path <- system.file("PAnalyzer.jar", package = "b10prot")
  tmp_input <- tempfile("pa_input_", fileext = ".tsv")
  tmp_output <- tempfile("pa_output_", fileext = ".tsv")

  write_delim(pep2prot, tmp_input, delim = "\t")
  sprintf("java -jar %s %s peptideRef proteinRef %s", jar_path, tmp_input, tmp_output) %>%
    system()
  output <- read_delim(tmp_output, delim = "\t") %>%
    add_class("panalyzer")

  file.remove(tmp_input)
  file.remove(tmp_output)

  return(output)
}

#' Summarize PAnalyzer Results
#'
#' Provides a summary of the PAnalyzer results by grouping the data based on
#' whether a protein is a target or decoy, and by its protein type. The summary
#' includes the number of distinct proteins and protein groups for both target
#' and decoy entries.
#'
#' @param panalyzer A data frame containing the results from the `panalyzer` function.
#'   This data frame must contain the following columns:
#'   - `isDecoy`: Logical indicator of whether the protein is a decoy (`TRUE`) or a target (`FALSE`).
#'   - `proteinType`: The type of protein as classified by PAnalyzer (e.g., `"conclusive"`, `"indistinguishable"`, `"ambiguous"`, `"non-conclusive"`).
#'   - `proteinRef`: Protein reference identifier.
#'   - `groupRef`: Protein group reference.
#'
#' @return A data frame summarizing the PAnalyzer results with the following columns:
#'   - `Type`: The type of protein (`"conclusive"`, `"indistinguishable"`, `"ambiguous"`, `"non-conclusive"`).
#'   - `TargetProteins`: The number of target proteins of the given type.
#'   - `DecoyProteins`: The number of decoy proteins of the given type.
#'   - `TargetGroups`: The number of target protein groups of the given type.
#'   - `DecoyGroups`: The number of decoy protein groups of the given type.
#'
#' @export
summary.panalyzer <- function(panalyzer) {
  panalyzer %>%
    group_by(isDecoy, proteinType) %>%
    summarise(proteins = n_distinct(proteinRef), groups = n_distinct(groupRef), .groups = "drop") %>%
    pivot_wider(names_from = isDecoy, values_from = c(proteins, groups)) %>%
    rename(Type = 1, TargetProteins = 2, DecoyProteins = 3, TargetGroups = 4, DecoyGroups = 5)
}

#' Plot Protein Groups Composition
#'
#' Visualizes the composition of selected protein groups based on their associated peptides.
#' The function filters the specified protein groups, arranges them, and displays a plot showing
#' the peptides associated with each protein in the group. Peptides are colored according to
#' their classification as unique, discriminating, non-discriminating, or non-significant.
#'
#' @param panalyzer A data frame containing the results from the `panalyzer` function. The data frame
#'   must include columns such as `groupRef`, `peptideRef`, `peptideType`, `proteinRef` and `proteinType`.
#' @param groupRefs A vector of group reference identifiers (`groupRef`) representing the protein
#'   groups to be visualized.
#'
#' @return A ggplot object visualizing the selected protein groups and their associated peptides.
#'
#' @examples
#' data(example_panalyzer, package = "b10prot")
#' plot_groups(example_panalyzer, groupRefs = 1:5)
#'
#' @seealso
#' - [PAnalyzer: A software tool for protein inference in shotgun proteomics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-288) for an example plot.
#'
#' @export
plot_groups <- function(panalyzer, groupRefs) {
  colors <- c(
    "unique" = "lightblue",
    "discriminating" = "lightgreen",
    "non-discriminating" = "lightgrey",
    "non-significant" = "white"
  )

  tmp_pep <-
    panalyzer %>%
    filter(groupRef %in% groupRefs) %>%
    distinct(peptideRef) %>%
    as_vector()

  tmp_grp <-
    panalyzer %>%
    filter(peptideRef %in% tmp_pep) %>%
    distinct(groupRef) %>%
    as_vector()

  panalyzer %>%
    filter(groupRef %in% tmp_grp) %>%
    arrange(groupRef) %>%
    mutate(peptideType = factor(peptideType, levels = c("unique", "discriminating", "non-discriminating", "non-significant"))) %>%
    mutate(proteinRef = factor(proteinRef, levels = sort(unique(panalyzer$proteinRef), TRUE))) %>%
    mutate(groupRef = paste0("Protein group #", groupRef, " (", proteinType, ")")) %>%
    group_by(proteinRef) %>%
    mutate(start = min(peptideRef), end = max(peptideRef)) %>%
    ungroup() %>%
    ggplot(aes(x = peptideRef, y = proteinRef, fill = peptideType)) +
    geom_segment(aes(x = start, xend = end, y = proteinRef, yend = proteinRef)) +
    geom_tile(width = 0.9, height = 0.9, color = "black") +
    geom_text(aes(label = ifelse(nchar(peptideRef) < 3, peptideRef, ""))) +
    scale_x_discrete(name = "Peptides", labels = function(x) str_trunc(x, 15)) +
    scale_y_discrete(name = "Proteins") +
    labs(title = "Protein groups composition") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    facet_col(vars(groupRef), scales = "free_y", space = "free")
}
