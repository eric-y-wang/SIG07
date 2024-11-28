library(tidyverse)
library(Matrix)
library(glmGamPoi)
# load parallelization libraries
library(furrr)
plan(multicore)
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB

# load data
counts <- readMM("/data/rudensky/EYW/SIG07/scanpy_outs/SIG07_doublets_CR_RNA_counts.mtx")
genes <- read_csv("/data/rudensky/EYW/SIG07/scanpy_outs/SIG07_doublets_CR_RNA_obs_genes.csv")
cells <- read_csv("/data/rudensky/EYW/SIG07/scanpy_outs/SIG07_doublets_CR_RNA_obs_cells.csv")
obs <- read_csv("/data/rudensky/EYW/SIG07/scanpy_outs/SIG07_doublets_CR_RNA_obs.csv")
degSig <- read_csv("/data/rudensky/EYW/git_projects/SIG07/analysis_outs/deg_sig_wilcoxon.csv")[,-1]

# add gene and cell names to counts
colnames(counts) <- genes$`...1`
rownames(counts) <- cells$cell_barcode

# transpose counts and pre-filter counts
counts <- t(counts)
counts <- counts[rownames(counts) %in% degSig$namse,]

get_gene_list <- function(ligands, deg_data, p_val_cutoff = 0.1) {
    # Function Purpose:
    # This function identifies a unique list of significant genes associated with specified ligand combinations
    # from differential expression data based on an adjusted p-value cutoff.
    # This function identifies DEGs that are present in any of the ligands or interaction (union).

    # Arguments:
    # - ligands: A character vector of two ligand names to define the groups of interest.
    # - deg_data: A data frame containing differential expression results, expected to have
    #             columns 'group', 'pvals_adj', and 'names'.
    # - p_val_cutoff: A numeric value (default = 0.1) specifying the adjusted p-value cutoff for significance.

    # Return:
    # - A character vector of unique gene names meeting the criteria.

    # Create a vector of group names based on the input ligands in a specific order.
    groupPull <- c(
        paste0(ligands[1], "_linker"),         # Group 1: First ligand + "_linker"
        paste0(ligands[2], "_linker"),         # Group 2: Second ligand + "_linker"
        paste0(ligands[1], "_", ligands[2]),   # Group 3: First ligand + "_" + Second ligand
        paste0(ligands[2], "_", ligands[1]),   # Group 4: Second ligand + "_" + First ligand
        "linker_linker"                        # Group 5: Default "linker_linker" group
    )

    # Filter for relevant groups and pull out significant genes.
    relevantDeg = deg_data %>% 
                    filter(group %in% groupPull) %>% 
                    filter(pvals_adj < p_val_cutoff)  

    # Extract the unique names of genes from the filtered data.
    genes = relevantDeg$names %>% unique()

    # Return the list of significant genes.
    return(genes)
}

glm_interaction_fit <- function(ligands, genes, counts, obs, pseudocount = 0.01) {
    # Function Purpose:
    # This function fits a generalized linear model (GLM) with an interaction term for two ligands.
    # It incorporates covariates like mitochondrial content, cell cycle scores, and scaling factors.
    # The model is fit using the `glmGamPoi` package.

    # Arguments:
    # - ligands: A character vector of two ligand names (e.g., c("ligandA", "ligandB")).
    # - genes: A character vector of gene names to subset from the count matrix.
    # - counts: A matrix of gene expression counts (rows: genes, columns: cells).
    # - obs: A data frame with metadata for cells (e.g., ligand calls, cell cycle scores).
    # - pseudocount: A small positive value added to the counts to avoid log-transform issues (default = 0.01).

    # Return:
    # - A fitted GLM object from `glmGamPoi`.

    # Get names of groups to pull.
    # The order of groupPull is important for downstream steps.
    groupPull <- c(
        paste0(ligands[1], "_linker"),         # Group 1: ligand 1 + "_linker"
        paste0(ligands[2], "_linker"),         # Group 2: ligand 2 + "_linker"
        paste0(ligands[1], "_", ligands[2]),   # Group 3: ligand 1 + "_" + ligand 2
        paste0(ligands[2], "_", ligands[1]),   # Group 4: ligand 2 + "_" + ligand 1
        "linker_linker"                        # Group 5: default "linker_linker"
    )

    # Subset the obs data frame to include only relevant groups.
    obsSub <- obs %>% filter(ligand_call_oBC_CR %in% groupPull)

    # Subset the counts matrix to include only the selected genes and cells from obsSub.
    # The transposition (`as.matrix`) is necessary for compatibility with `glmGamPoi`.
    countsSub <- counts[genes, obsSub$cell_barcode] %>% as.matrix()

    # Add a small pseudocount to the count data to prevent issues with log-transforming zero counts.
    countsSub <- countsSub + pseudocount

    # Define binary indicators for the ligands:
    # - ligand1 is TRUE for cells exposed to ligand 1 (including double ligand groups).
    # - ligand2 is TRUE for cells exposed to ligand 2 (including double ligand groups).
    ligand1 <- obsSub$ligand_call_oBC_CR %in% c(groupPull[c(1, 3, 4)])
    ligand2 <- obsSub$ligand_call_oBC_CR %in% c(groupPull[c(2, 3, 4)])

    # Extract additional covariates for the model.
    group <- obsSub$group_call_CR %>% as.factor()        # Categorical group assignments.
    percent.mito <- obsSub$pct_counts_mt %>% as.numeric()  # Percentage of mitochondrial counts.
    s.score <- obsSub$S_score %>% as.numeric()          # Cell cycle S-phase score.
    g2m.score <- obsSub$G2M_score %>% as.numeric()      # Cell cycle G2M-phase score.
    scaling.factor <- obsSub$total_counts %>% as.numeric()  # Scaling factor for library size.

    # Combine covariates into a data frame for the model.
    model.df <- cbind(
        ligand1,
        ligand2,
        group,
        percent.mito,
        s.score,
        g2m.score,
        scaling.factor
    ) %>% as.data.frame()

    # Define the formula for the GLM.
    # The interaction term `ligand1*ligand2` allows for detecting combinatorial effects.
    model.formula <- as.formula(paste0(
        '~ ligand1*ligand2 +',  # Interaction between ligands 1 and 2.
        'group + ',             # Group-specific effects.
        'percent.mito + ',      # Mitochondrial content.
        's.score + ',           # S-phase cell cycle score.
        'g2m.score'             # G2M-phase cell cycle score.
    ))

    # Fit the generalized linear model with glmGamPoi.
    fit <- glm_gp(
        countsSub,                # Subsetted count matrix.
        design = model.formula,   # Model formula specifying covariates and interactions.
        col_data = model.df,      # Data frame containing model covariates.
        offset = log(model.df$scaling.factor),  # Offset term (log of scaling factor).
        size_factors = FALSE,     # Use explicit scaling factors instead of size factor normalization.
        on_disk = FALSE           # Fit the model in memory.
    )

    # Return the fitted model.
    return(fit)
}

run_glm_test <- function(ligands, counts, obs, deg_data) {
    # Function Purpose:
    # This function performs differential expression (DE) testing for genes that show an interaction 
    # between two specified ligands using a generalized linear model (GLM).

    # Arguments:
    # - ligands: A character vector of two ligand names (e.g., c("ligandA", "ligandB")).
    # - counts: A matrix of gene expression counts (rows: genes, columns: cells).
    # - obs: A data frame with metadata for cells (e.g., ligand calls, cell cycle scores).
    # - deg_data: A data frame with differential expression results, used to filter genes to test.

    # Return:
    # - A tibble with DE results, including interaction terms for the ligands.

    # Step 1: Get genes to test
    IRdisplay::display("Step 1: Extracting genes to test based on ligands.")
    genes <- get_gene_list(ligands, deg_data)
    IRdisplay::display(paste("Number of genes to test:", length(genes)))

    # Step 2: Fit the interaction model
    IRdisplay::display("Step 2: Fitting the interaction model using glm_interaction_fit.")
    fit <- glm_interaction_fit(
        ligands = ligands,
        genes = genes,
        counts = counts,
        obs = obs
    )
    IRdisplay::display("Model fitting complete.")

    # Step 3: Perform DE testing
    IRdisplay::display("Step 3: Performing differential expression testing for the interaction term.")
    res <- test_de(fit, contrast = "ligand1:ligand2") %>%  # Test interaction term
        as_tibble() %>%                                   # Convert to tibble
        mutate(interaction = paste0(ligands[1], "_", ligands[2]))  # Annotate with interaction name
    IRdisplay::display("DE testing complete.")

    # Step 4: Return results
    IRdisplay::display("Returning DE results.")
    return(res)
}

# get unique ligands 
calls <- obs$ligand_call_oBC_CR %>%
    unique() %>%
    str_split("_") %>% 
    unlist() %>%
    unique()
calls <- calls[calls != 'linker']
calls

# Generate non-repeating pairwise combinations
pairwise_combinations <- combn(calls, 2, simplify = FALSE)

# Run the function in parallel and combine results
results <- future_map_dfr(
  pairwise_combinations, 
  ~ run_glm_test(ligands = .x, counts = counts, obs = obs, deg_data = degSig),
  .progress = TRUE
)

write_csv(results,"/data/rudensky/EYW/git_projects/SIG07/analysis_outs/glmGamPoi_interaction_deg.csv")