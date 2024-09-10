# Utility functions for the different steps of a proteomics identification workflow

load_psms <- function(mzid_file, verbose = TRUE) {
    if( verbose )
        sprintf("Loading %s...", mzid_file) %>%
            print()
    data_mzid <-
        openIDfile(mzid_file)
    if( verbose )
        mzidInfo(data_mzid) %>%
            print()
    data_psms <-
        data_mzid %>%
        psms() %>%
        as_tibble()
        #mutate(SpectraSource = basename(mzidInfo(data_mzid)$SpectraSource), .before = 1)
    data_scores <-
        data_mzid %>%
        score() %>%
        as_tibble()
    data_psms %>%
        bind_cols(data_scores %>% select(-spectrumID))
}

#' Load PSMs from mzIdentML Files
#'
#' This function loads Peptide-to-Spectrum Matches (PSMs) from mzIdentML files in a specified directory,
#' combines them, and adds a standard column for the PSM score and protein reference.
#'
#' @param path A character string specifying the directory where the mzIdentML files are located.
#'   Default is the current working directory (`"."`).
#' @param pattern A character string representing the file pattern to search for. Default is `".mzid"`,
#'   which targets mzIdentML files.
#' @param psm_score An optional character string specifying the column name to be used as the PSM score.
#'   If `NULL`, the last column in the loaded data will be used as the PSM score.
#'
#' @return A `data.frame` containing the combined PSM data, with added columns for the PSM score
#'   (`psmScore`) and the protein reference (`proteinRef`).
#'
#' @export
wf_load_psms <- function(path = ".", pattern = ".mzid", psm_score = NULL, verbose = FALSE) {
    psms <-
        # Find files with PSMs
        list.files(path, pattern, full.names = TRUE) %>%
        # Loads PSMs from each single file
        map(~ load_psms(.x, verbose = verbose)) %>%
        # Concatenate PSMs from different files (eg. fractions)
        bind_rows() %>%
        # Create a standard column for the psm score
        mutate(psmScore =
            if(is.null(psm_score)) .data[[last(colnames(.))]]
            else .data[[psm_score]]
        ) %>%
        # Create a standard column with the protein reference
        mutate(proteinRef = DatabaseAccess)
}

#' Aggregate PSMs into Peptides
#'
#' This function aggregates Peptide-to-Spectrum Matches (PSMs) into peptides, selecting the best
#' PSM score for each peptide and maintaining decoy information. Peptides that are present in both
#' the target and decoy subsets are removed from the final results.
#'
#' @param psms A `data.frame` containing the PSM data, which must include columns `peptideRef`, `psmScore`, and `isDecoy`.
#' @param lower_better A logical value indicating whether a lower PSM score is better. Default is `TRUE`.
#'
#' @return A `data.frame` containing one row per unique peptide with columns:
#' \item{`peptideRef`}{The peptide identifier.}
#' \item{`pepScore`}{The best score for the peptide.}
#' \item{`isDecoy`}{A logical indicating whether the peptide belongs to the decoy set.}
#'
#' @seealso [wf_load_psms] for loading PSMs.
#'
#' @export
wf_psm2pep <- function(psms, lower_better = TRUE) {
    psms %>%
        # peptideRef is the primary key and will be unique in the list
        group_by(peptideRef) %>%
        summarise(
            # We use the best PSM per peptide
            pepScore = if(lower_better) min(psmScore) else max(psmScore),
            # Check if the peptide is present in both the decoy and target subsets
            subsets = n_distinct(isDecoy),
            # We keep decoy information (mandatory)
            isDecoy = any(isDecoy),
            .groups = "drop"
        ) %>%
        # Remove peptides present in both target and decoy subsets
        filter(subsets == 1) %>%
        select(-subsets)
}

#' Map Peptides to a Specified Reference Level
#'
#' This function maps peptides to a specified reference level (e.g., protein, gene, or protein group)
#' based on Peptide-to-Spectrum Matches (PSMs). It also calculates how many reference entities
#' (e.g., proteins, genes, or groups) are matched by each peptide.
#'
#' @param psms A `data.frame` containing PSM data, including `peptideRef` (peptide identifier) and the reference level specified by `levelRef`.
#' @param levelRef The column name in `psms` that contains the reference level to which peptides should be mapped.
#'
#' @return A `data.frame` containing:
#' \item{`peptideRef`}{The peptide identifier.}
#' \item{`levelRef`}{The reference level identifier.}
#' \item{`shared`}{The number of reference entities matched by each peptide.}
#'
#' @export
wf_pep2level <- function(psms, levelRef) {
    psms %>%
        # We obtain a list of peptide-to-protein relations from PSMs
        group_by(peptideRef, {{levelRef}}) %>%
        summarise(.groups = "drop") %>%
        # Calculate the number of proteins matched by each peptide
        group_by(peptideRef) %>%
        mutate(shared = n()) %>%
        ungroup()
}

#' Perform Protein Grouping Based on Peptide-to-Protein Relations
#'
#' This function builds protein groups using the peptides that pass a
#' specified peptide-level FDR threshold. Peptides that do not pass the threshold
#' are grouped separately and assigned a negative group identifier. It utilizes
#' the PAnalyzer algorithm to infer protein groups and assigns peptide and protein types.
#'
#' @param pep2prot A data frame containing peptide-to-protein relations, including
#'   the columns `peptideRef`, `proteinRef`, `qval`, and `isDecoy`.
#' @param threshold A numeric value specifying the peptide-level FDR threshold
#'   for grouping (default is 0.01).
#'
#' @return A data frame with the inferred protein groups, with additional columns:
#'   - `peptideType`: Type of peptide (`"unique"`, `"discriminating"`, `"non-discriminating"`).
#'   - `proteinType`: Type of protein (`"conclusive"`, `"indistinguishable"`, `"ambiguous"`, `"non-conclusive"`).
#'   - `groupRef`: Group identifier for the proteins.
#'   - `shared`: The number of groups matched by each peptide.
#'
#' @seealso \code{\link{panalyzer}} for protein grouping, \code{\link{wf_pep2level}} for obtaining peptide-to-protein relations.
#' @export
wf_grouping <- function(pep2prot, threshold = 0.01) {
    pa_ok <-
        pep2prot %>%
        # Builg groups using only peptides passing the FDR threshold
        filter(qval <= 0.01) %>%
        # Filter unnecessary information
        select(peptideRef, proteinRef, isDecoy) %>%
        panalyzer() %>%
        right_join(pep2prot) %>%
        # Fill missing protein-related fields in peptide-to-protein relations not used for building the groups
        group_by(proteinRef) %>%
        mutate(across(c(proteinType, groupRef), ~ first(na.omit(.x)))) %>%
        ungroup() %>%
        mutate(peptideType = ifelse(is.na(peptideType),"non-significant",peptideType)) %>%
        # Calculate the number of groups matched by each peptide
        group_by(peptideRef) %>%
        mutate(shared = ifelse(any(is.na(groupRef)), NA, n_distinct(groupRef))) %>%
        ungroup() %>%
        filter(!is.na(groupRef))

    pa_ko <-
        pep2prot %>%
        filter(!peptideRef %in% pa_ok$peptideRef) %>%
        select(peptideRef, proteinRef, isDecoy) %>%
        panalyzer() %>%
        inner_join(pep2prot) %>%
        group_by(peptideRef) %>%
        mutate(shared = n_distinct(groupRef)) %>%
        ungroup() %>%
        mutate(groupRef = -groupRef)

    bind_rows(pa_ok, pa_ko) %>%
        add_class("panalyzer")
}

#' Assign Peptides to Protein Groups
#'
#' This function processes peptide-to-protein group relationships by summarizing
#' the number of discriminating and total peptides per protein, and then assigning
#' peptides to their respective protein groups. Proteins within a group are ordered
#' by the number of discriminating peptides, followed by the total number of peptides.
#'
#' @param pep2prot2group A data frame containing peptide-to-protein-to-group relations,
#'   including `peptideRef`, `proteinRef`, `groupRef`, and `peptideType`.
#'
#' @return A data frame where peptides are assigned to protein groups, with additional columns:
#'   - `proteinCount`: The number of proteins in each group.
#'   - `proteinRefs`: A concatenation of all protein references in the group.
#'   - `proteinMaster`: The first protein of the group.
#'
#' @seealso \code{\link{wf_grouping}} for generating peptide-to-protein-to-group relations,
#'   and \code{\link{panalyzer}} for protein grouping.
#' @export
wf_pep2group <- function(pep2prot2group) {
    pep2prot2group %>%
        group_by(proteinRef) %>%
        mutate(
            discPeptides = sum(peptideType == "discriminating"),
            totalPeptides = n()
        ) %>%
        ungroup() %>%
        group_by(groupRef) %>%
        arrange(desc(discPeptides), desc(totalPeptides), proteinRef) %>%
        mutate(proteinCount = n_distinct(proteinRef)) %>%
        mutate(proteinRefs = paste0(unique(proteinRef), collapse = ";")) %>%
        mutate(proteinMaster = first(proteinRef)) %>%
        ungroup() %>%
        group_by(peptideRef, groupRef) %>%
        summarise(across(everything(), first), .groups = "drop")
}
