load_psms <- function(mzid_file) {
    sprintf("Loading %s...", mzid_file) %>%
        print()
    data_mzid <-
        openIDfile(mzid_file)
    mzidInfo(data_mzid) %>%
        print()
    data_psms <-
        data_mzid %>%
        psms() %>%
        as_tibble() %>%
        mutate(SpectraSource = basename(mzidInfo(data_mzid)$SpectraSource), .before = 1)
    data_scores <-
        data_mzid %>%
        score() %>%
        as_tibble()
    data_psms %>%
        bind_cols(data_scores %>% select(-spectrumID))
}

wf_load_psms <- function(path = ".", pattern = ".mzid", psm_score = NULL) {
    psms <-
        # Find files with PSMs
        list.files(path, pattern, full.names = TRUE) %>%
        # Loads PSMs from each single file
        map(load_psms) %>%
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
            isDecoy = any(isDecoy)
        ) %>%
        # Remove peptides present in both target and decoy subsets
        filter(subsets == 1) %>%
        select(-subsets)
}

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
