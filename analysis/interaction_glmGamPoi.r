library(tidyverse)
library(Matrix)
library(glmGamPoi)
library(furrr)

# Set up parallel processing
plan(multicore)
options(future.globals.maxSize = 8 * 1024^3)  # Increase memory limit to 8 GB

# Load data
cat("Loading data...\n")
counts <- readMM("/data/rudensky/EYW/SIG07/scanpy_outs/SIG07_doublets_CR_RNA_counts.mtx")
genes <- read_csv("/data/rudensky/EYW/SIG07/scanpy_outs/SIG07_doublets_CR_RNA_obs_genes.csv")
cells <- read_csv("/data/rudensky/EYW/SIG07/scanpy_outs/SIG07_doublets_CR_RNA_obs_cells.csv")
obs <- read_csv("/data/rudensky/EYW/SIG07/scanpy_outs/SIG07_doublets_CR_RNA_obs.csv")
degSig <- read_csv("/data/rudensky/EYW/git_projects/SIG07/analysis_outs/deg_sig_wilcoxon.csv")[, -1]
cat("Data loaded successfully.\n")

# Add gene and cell names to counts
cat("Processing counts matrix...\n")
colnames(counts) <- genes$`...1`
rownames(counts) <- cells$cell_barcode

# Transpose counts and pre-filter counts
counts <- t(counts)
counts <- counts[rownames(counts) %in% degSig$names, ]
cat("Counts matrix processed.\n")

# Define functions
get_gene_list <- function(ligands, deg_data, p_val_cutoff = 0.1) {
    cat("Getting gene list for ligands: ", paste(ligands, collapse = ", "), "\n")
    groupPull <- c(
        paste0(ligands[1], "_linker"),
        paste0(ligands[2], "_linker"),
        paste0(ligands[1], "_", ligands[2]),
        paste0(ligands[2], "_", ligands[1]),
        "linker_linker"
    )
    relevantDeg <- deg_data %>%
        filter(group %in% groupPull) %>%
        filter(pvals_adj < p_val_cutoff)
    genes <- relevantDeg$names %>% unique()
    cat("Number of genes found: ", length(genes), "\n")
    return(genes)
}

glm_interaction_fit <- function(ligands, genes, counts, obs, pseudocount = 0.01) {
    cat("Fitting GLM for ligands: ", paste(ligands, collapse = ", "), "\n")
    groupPull <- c(
        paste0(ligands[1], "_linker"),
        paste0(ligands[2], "_linker"),
        paste0(ligands[1], "_", ligands[2]),
        paste0(ligands[2], "_", ligands[1]),
        "linker_linker"
    )
    obsSub <- obs %>% filter(ligand_call_oBC_CR %in% groupPull)

    # Check if there are valid cells and genes
    obsSub <- obsSub[obsSub$cell_barcode %in% colnames(counts), ]
    genes <- genes[genes %in% rownames(counts)]

    if (nrow(obsSub) == 0 || length(genes) == 0) {
        cat("No valid cells or genes for ligands: ", paste(ligands, collapse = ", "), "\n")
        return(NULL)
    }

    countsSub <- counts[genes, obsSub$cell_barcode, drop = FALSE] %>% as.matrix()
    countsSub <- countsSub + pseudocount

    ligand1 <- as.numeric(obsSub$ligand_call_oBC_CR %in% groupPull[c(1, 3, 4)])
    ligand2 <- as.numeric(obsSub$ligand_call_oBC_CR %in% groupPull[c(2, 3, 4)])
    group <- obsSub$group_call_CR %>% as.factor()
    percent.mito <- obsSub$pct_counts_mt %>% as.numeric()
    s.score <- obsSub$S_score %>% as.numeric()
    g2m.score <- obsSub$G2M_score %>% as.numeric()
    scaling.factor <- obsSub$total_counts %>% as.numeric()

    model.df <- cbind(
        ligand1,
        ligand2,
        group,
        percent.mito,
        s.score,
        g2m.score,
        scaling.factor
    ) %>% as.data.frame()

    model.formula <- as.formula(
        '~ ligand1*ligand2 + group + percent.mito + s.score + g2m.score'
    )

    fit <- glm_gp(
        countsSub,
        design = model.formula,
        col_data = model.df,
        offset = log(model.df$scaling.factor),
        size_factors = FALSE,
        on_disk = FALSE
    )
    cat("GLM fitting complete for ligands: ", paste(ligands, collapse = ", "), "\n")
    return(fit)
}

run_glm_test <- function(ligands, counts, obs, deg_data) {
    genes <- get_gene_list(ligands, deg_data)

    # Check if any genes were found
    if (length(genes) == 0) {
        # If no genes are found, log a message and return NULL
        message("No genes to test for ligands: ", paste(ligands, collapse = ", "))
        return(NULL)
    }

    # Step 2: Fit the GLM with interaction terms
    # Use the `glm_interaction_fit` function to fit a generalized linear model.
    # The model includes interaction terms for the two ligands and adjusts for covariates.
    fit <- tryCatch(
        glm_interaction_fit(
            ligands = ligands,
            genes = genes,
            counts = counts,
            obs = obs
        ),
        error = function(e) {
            # Handle errors in model fitting: log the error and return NULL
            message("Error fitting model for ligands: ", paste(ligands, collapse = ", "), ": ", e$message)
            return(NULL)
        }
    )

    # Check if the model fitting was successful
    if (is.null(fit)) {
        # If model fitting failed, return NULL
        return(NULL)
    }

    # Step 3.1: Extract DE test results interaction
    # Use the `test_de` function to perform differential expression testing for the interaction term.
    res1 <- tryCatch(
        test_de(fit, contrast = `ligand1:ligand2`) %>%
            as_tibble() %>%                           
            mutate(interaction = paste0(ligands[1], "_", ligands[2])),
        error = function(e) {
            # Handle errors in DE testing: log the error and return NULL
            message("Error in test_de for ligand interaction: ", paste(ligands, collapse = ", "), ": ", e$message)
            return(NULL)
        }
    )

    # Step 3.2: Extract DE test results ligand 1
    # Use the `test_de` function to perform differential expression testing for the interaction term.
    res2 <- tryCatch(
        test_de(fit, contrast = `ligand1`) %>%
            as_tibble() %>%                           
            mutate(interaction = paste0(ligands[1],"_",ligands[2]),
		   single_ligand = paste0(ligands[1])),
        error = function(e) {
            # Handle errors in DE testing: log the error and return NULL
            message("Error in test_de for ligand 1: ", paste(ligands, collapse = ", "), ": ", e$message)
            return(NULL)
        }
    )

    # Step 3.3: Extract DE test results ligand 1
    # Use the `test_de` function to perform differential expression testing for the interaction term.
    res3 <- tryCatch(
        test_de(fit, contrast = `ligand2`) %>%
            as_tibble() %>%                           
            mutate(interaction = paste0(ligands[1],"_",ligands[2]),
		   single_ligand = paste0(ligands[2])),
        error = function(e) {
            # Handle errors in DE testing: log the error and return NULL
            message("Error in test_de for ligands 2: ", paste(ligands, collapse = ", "), ": ", e$message)
            return(NULL)
        }
    )

    # Step 4: Extract model coefficients
    # Extract the coefficients from the fitted model and convert them to a tibble.
    coef_tibble <- tryCatch(
        as_tibble(fit$Beta) %>%
            mutate(genes = rownames(fit$Beta),
                   interaction = paste0(ligands[1], "_", ligands[2])),  # Add ligand pair as metadata
        error = function(e) {
            # Handle errors in coefficient extraction: log the error and return NULL
            message("Error extracting coefficients for ligands: ", paste(ligands, collapse = ", "), ": ", e$message)
            return(NULL)
        }
    )

    # Combine DE results and coefficients into a single list
    # `de_results` contains the differential expression results
    # `coefficients` contains the wide-format coefficients
    combined_results <- list(
        de_results = res1,
	de_results_singles = bind_rows(res2,res3),
        coefficients = coef_tibble
    )

    # Return the combined results
    return(combined_results)
}

# Get unique ligands
cat("Generating pairwise ligand combinations...\n")
calls <- obs$ligand_call_oBC_CR %>%
    unique() %>%
    str_split("_") %>%
    unlist() %>%
    unique()
calls <- calls[calls != "linker"]
pairwise_combinations <- combn(calls, 2, simplify = FALSE)
cat("Pairwise combinations generated.\n")

# Run the function in parallel and combine results
cat("Running GLM tests in parallel...\n")
results <- future_map(
    pairwise_combinations,
    ~ run_glm_test(ligands = .x, counts = counts, obs = obs, deg_data = degSig),
    .progress = TRUE
)

# Combine and save results
cat("Combining and saving results...\n")
interaction_de_results <- bind_rows(lapply(results, function(x) x$de_results))
single_de_results <- bind_rows(lapply(results, function(x) x$de_results_singles))
all_coefficients <- bind_rows(lapply(results, function(x) x$coefficients))
write_csv(interaction_de_results, "/data/rudensky/EYW/git_projects/SIG07/analysis_outs/glmGamPoi_interaction_deg.csv")
write_csv(single_de_results, "/data/rudensky/EYW/git_projects/SIG07/analysis_outs/glmGamPoi_singles_deg.csv")
write_csv(all_coefficients, "/data/rudensky/EYW/git_projects/SIG07/analysis_outs/glmGamPoi_coefficients.csv")
cat("Results saved successfully.\n")
