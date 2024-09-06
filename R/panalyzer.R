add_class <- function(obj, new_class) {
  class(obj) <- c(new_class, class(obj))
  obj
}

panalyzer <- function(pep2prot) {
  write_delim(pep2prot, "tmp_pa_input.tsv", delim = "\t")
  system("java -jar PAnalyzer.jar tmp_pa_input.tsv peptideRef proteinRef tmp_pa_output.tsv")
  read_delim("tmp_pa_output.tsv", delim = "\t") %>%
    add_class("panalyzer")
}

summary.panalyzer <- function(panalyzer) {
  panalyzer %>%
    group_by(isDecoy, proteinType) %>%
    summarise(proteins = n_distinct(proteinRef), groups = n_distinct(groupRef), .groups = "drop") %>%
    pivot_wider(names_from = isDecoy, values_from = c(proteins, groups)) %>%
    rename(Type = 1, TargetProteins = 2, DecoyProteins = 3, TargetGroups = 4, DecoyGroups = 5)
}

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
