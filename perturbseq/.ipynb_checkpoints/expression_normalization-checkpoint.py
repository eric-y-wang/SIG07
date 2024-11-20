# Perturbseq library for loading and manipulating single-cell experiments
# Copyright (C) 2019  Thomas Norman

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import pandas as pd
import numpy as np
from sklearn import preprocessing as pre
from pandas.api.types import is_numeric_dtype
from six.moves import zip as izip
from time import time
import gc
from tqdm import tqdm_notebook
import scipy as sp

# MODIFIED FUNCTIONS

def normalize_matrix_to_control(matrix, control_matrix, scale_by_total=True, median_umi_count=None):
    """ Normalize expression distribution relative to a population of control (unperturbed) cells.
    The normalization proceeds by first normalizing the UMI counts (if umi_count != None) within each cell to the median
    UMI count within the population. The control population is normalized to the same amount. 
    The expression within the population is then Z-normalized with respect to the control 
    distribution: i.e., for each gene the control mean is subtracted, and then these values are
    divided by the control standard deviation.
        
    Args:
        matrix: gene expression matrix to normalize (output from cellranger)
        control_matrix: gene expression matrix of control population
        scale_by_total: Rescale UMI counts within each cell to population median (default: True)
        median_umi_count: If provided, set the median to normalize to. This is useful for example
            when normalizing multiple independent lanes/gemgroups.
    
    Returns:
        DataFrame of normalized expression data
    """
    # Convert sparse matrices to dense if necessary
    if sp.sparse.issparse(matrix):
        print('     Densifying matrix...')
        matrix = matrix.todense()
    if sp.sparse.issparse(control_matrix):
        print('     Densifying control matrix...')
        control_matrix = control_matrix.todense()

    # Convert to numpy arrays for processing if needed
    if isinstance(matrix, pd.DataFrame):
        matrix = matrix.values
    if isinstance(control_matrix, pd.DataFrame):
        control_matrix = control_matrix.values

    # Normalize to median UMI count if specified
    if scale_by_total:
        print("     Determining scale factors...")
        reads_per_bc = matrix.sum(axis=1)
        
        if median_umi_count is None:
            median_reads_per_bc = np.median(reads_per_bc)
        else:
            median_reads_per_bc = median_umi_count
        
        scaling_factors = median_reads_per_bc / reads_per_bc
        
        print("     Normalizing matrix to median")
        m = matrix.astype(np.float64)
        # Normalize expression within each cell by median total count
        m *= scaling_factors[:, np.newaxis]
        if np.mean(median_reads_per_bc) < 5000:
            print("Scaling with a small number of reads. Are you sure this is what you want?")
            
        control_reads_per_bc = control_matrix.sum(axis=1)
        
        print("     Normalizing control matrix to median")
        control_scaling_factors = median_reads_per_bc / control_reads_per_bc
        c_m = control_matrix.astype(np.float64)
        c_m *= control_scaling_factors[:, np.newaxis]
    else:
        m = matrix.astype(np.float64)
        c_m = control_matrix.astype(np.float64)

    control_mean = c_m.mean(axis=0)
    control_std = c_m.std(axis=0)
    
    print("     Scaling matrix to control")
    # Center and rescale the expression of each gene to average 0 and std 1
    m_out = (m - control_mean) / control_std
    
    print("     Done.")
    return pd.DataFrame(m_out, columns=np.arange(m.shape[1]), index=np.arange(m.shape[0]))

def normalize_to_gemgroup_control_residuals_adata(adata, control_cells_query, **kwargs):
    """
    Normalizes a multi-lane 10x experiment stored in an AnnData object.
    Cells within each gemgroup are normalized to the control cells within the same gemgroup.

    Args:
        adata: AnnData object to normalize
        control_cells_query: String query to identify control cell population
            to normalize with respect to.
        **kwargs: Additional arguments passed to groupby on adata.obs, useful for refined slicing.

    Returns:
        DataFrame of normalized expression data with appropriate indices.
    """
    # Check if 'gem_group' exists and determine unique gem groups
    if 'gem_group' not in adata.obs.columns:
        raise ValueError("'gem_group' column not found in adata.obs.")
    
    unique_gem_groups = adata.obs['gem_group'].unique()
    
    if len(unique_gem_groups) == 1:
        # Only one gem group, normalize directly
        print("Single gem group detected. Normalizing without grouping.")
        gemgroup_pop = adata.obs
        gemgroup_control_pop = adata[adata.obs.query(control_cells_query).index].obs

        gemgroup_data = adata[gemgroup_pop.index].X.copy()
        control_data = adata[gemgroup_control_pop.index].X.copy()

        normalized_matrix = normalize_matrix_to_control(
            gemgroup_data,
            control_data,
            scale_by_total=False
        )
        
        # Assign indices and columns
        normalized_matrix.index = gemgroup_pop.index
        normalized_matrix.columns = adata.var_names

        return normalized_matrix
    else:
        # Multiple gem groups, normalize for each gem group
        print("Multiple gem groups detected. Processing each group separately.")
        
        gemgroup_iterator = zip(
            adata.obs.groupby('gem_group', observed=False),
            adata[adata.obs.query(control_cells_query).index].obs.groupby('gem_group', observed=False)
        )
        
        gem_group_matrices = dict()
        for (i, gemgroup_pop), (_, gemgroup_control_pop) in gemgroup_iterator:
            print('Processing gem group {0}'.format(i))
            t = time()

            # Select the data for the current gemgroup
            gemgroup_data = adata[gemgroup_pop.index].X.copy()
            control_data = adata[gemgroup_control_pop.index].X.copy()

            gem_group_matrices[i] = normalize_matrix_to_control(
                gemgroup_data,
                control_data,
                scale_by_total=False
            )
            
            # Assign the correct index
            gem_group_matrices[i].index = gemgroup_pop.index
            gem_group_matrices[i].columns = adata.var_names

            print(time() - t)
        
        print('Merging submatrices...')
        
        # Merge all matrices into a DataFrame
        normalized_matrix = pd.concat(gem_group_matrices.values())
        
        # Ensure the rows match the original AnnData index
        normalized_matrix = normalized_matrix.loc[adata.obs_names]
        
        return normalized_matrix

def normalize_to_control_residuals_adata(adata, control_cells_query, groupby_column, **kwargs):
    """
    Normalizes a multi-lane 10x experiment stored in an AnnData object.
    Cells within each group (defined by a column in `.obs`) are normalized to the control cells within the same group.

    Args:
        adata: AnnData object to normalize.
        control_cells_query: String query to identify control cell population to normalize with respect to.
        groupby_column: Column name in `adata.obs` to use for grouping cells (e.g., "gem_group").
        **kwargs: Additional arguments passed to groupby on `adata.obs`, useful for refined slicing.

    Returns:
        DataFrame of normalized expression data with appropriate indices.
    """
    # Check if the groupby column exists in .obs
    if groupby_column not in adata.obs.columns:
        raise ValueError(f"'{groupby_column}' column not found in adata.obs.")
    
    unique_groups = adata.obs[groupby_column].unique()
    
    if len(unique_groups) == 1:
        # Only one group, normalize directly
        print(f"Single {groupby_column} detected. Normalizing without grouping.")
        group_pop = adata.obs
        control_indices = adata.obs.query(control_cells_query).index

        if control_indices.empty:
            raise ValueError(f"No control cells found for the query: {control_cells_query}")

        group_data = adata[group_pop.index].X.copy()
        control_data = adata[control_indices].X.copy()

        normalized_matrix = normalize_matrix_to_control(
            group_data,
            control_data,
            scale_by_total=False
        )
        
        # Assign indices and columns
        normalized_matrix.index = group_pop.index
        normalized_matrix.columns = adata.var_names

        return normalized_matrix
    else:
        # Multiple groups, normalize for each group
        print(f"Multiple {groupby_column} detected. Processing each group separately.")
        
        group_iterator = zip(
            adata.obs.groupby(groupby_column, observed=False),
            adata[adata.obs.query(control_cells_query).index].obs.groupby(groupby_column, observed=False)
        )
        
        group_matrices = dict()
        for (group_name, group_pop), (_, group_control_pop) in group_iterator:
            print(f'Processing {groupby_column} {group_name}')
            
            # Check if control population is empty for the group
            if group_control_pop.empty:
                print(f"Warning: No control cells found for {groupby_column} {group_name}. Skipping normalization.")
                continue

            # Select the data for the current group
            group_data = adata[group_pop.index].X.copy()
            control_data = adata[group_control_pop.index].X.copy()

            if control_data.shape[0] == 0:
                print(f"Warning: Control data matrix is empty for {groupby_column} {group_name}. Skipping normalization.")
                continue

            group_matrices[group_name] = normalize_matrix_to_control(
                group_data,
                control_data,
                scale_by_total=False
            )
            
            # Assign the correct index
            group_matrices[group_name].index = group_pop.index
            group_matrices[group_name].columns = adata.var_names
        
        print('Merging submatrices...')
        
        # Merge all matrices into a DataFrame
        normalized_matrix = pd.concat(group_matrices.values(), axis=0)
        
        # Ensure the rows match the original AnnData index
        normalized_matrix = normalized_matrix.loc[adata.obs_names]
        
        return normalized_matrix

    
# UNMODIFIED FUNCTIONS

def strip_low_expression(pop, threshold=0):
    """Remove genes with low or zero expression to reduce memory usage. Modifies the
    target CellPopulation in place.
    
    Args:
        pop: CellPopulation instance
        threshold: all genes with expression <= threshold will be removed
    """
    retain = pop.genes.query('mean > @threshold').index
    remove = pop.genes.query('mean <= @threshold').index
    if len(remove) == 0:
        print('No genes have expression below threshold.')
        return
    pop.matrix = pop.matrix[retain]
    pop.genes.loc[remove, 'in_matrix'] = False
    # set all numeric properties to nan for genes that have been removed 
    for col in np.setdiff1d(pop.genes.columns, ['gene_name', 'in_matrix']):
        if is_numeric_dtype(pop.genes[col]):
            pop.genes.loc[remove, col] = np.nan
    # force garbage collection
    gc.collect()
                
def equalize_UMI_counts(matrix, median_umi_count=None):
    """Normalize all cells in an expression matrix to a specified UMI count
    """
    reads_per_bc = matrix.sum(axis=1)
    if median_umi_count is None:
        median_reads_per_bc = np.median(reads_per_bc)
    else:
        median_reads_per_bc = median_umi_count
    scaling_factors = median_reads_per_bc / reads_per_bc
    m = matrix.astype(np.float64)
    # Normalize expression within each cell by median total count
    m = m.mul(scaling_factors, axis=0)
    if np.mean(median_reads_per_bc) < 5000:
        print("Scaling with a small number of reads. Are you sure this is what you want?")
    return m

def log_normalize_expression(pop, scale_by_total=True, pseudocount=1):
    """ Normalize expression distribution by log transformation.
    The normalization proceeds by first (optionally) normalizing the UMI counts within each cell 
    to the median UMI count within the population. The expression within the population is then 
    log-normalized: i.e., transformed according to Y = log2(X + 1)
       
    Args:
        pop: population to normalize
        scale_by_total: Rescale UMI counts within each cell to population median (default: True)
        pseudocount: offset for 0 values (default: 1)
        
    Returns:
        DataFrame of log-normalized expression data
    """
    matrix = pop.matrix
    
    if (scale_by_total):
        m = equalize_UMI_counts(matrix)
    else:
        m = matrix.astype(np.float64) 

    m = np.log2(m + pseudocount)
    
    return pd.DataFrame(m, columns=m.columns, index=m.index)
    
def z_normalize_expression(pop, scale_by_total=True):
    """ Normalize expression distribution by Z-scoring.
    The normalization proceeds by first normalizing the UMI counts within each cell to the 
    median UMI count within the population. The expression within the population is then 
    Z-normalized: i.e., for each gene the mean is subtracted, and then these values are divided 
    by that gene's standard deviation.
       
    Args:
        pop: population to normalize
        scale_by_total: Rescale UMI counts within each cell to population median (default: True)
    
    Returns:
        DataFrame of Z-normalized expression data
    """
    matrix = pop.matrix
    
    if (scale_by_total):
        m = equalize_UMI_counts(matrix)
    else:
        m = matrix.astype(np.float64) 

    # Now center and rescale the expression of each gene to average 0 and std 1
    m_out = pre.scale(m.as_matrix(), axis=0)
    
    return pd.DataFrame(m_out, columns=m.columns, index=m.index)

def normalize_to_control(pop, control_cells, scale_by_total=True, median_umi_count=None, **kwargs):
    """ Normalize expression distribution relative to a population of control (unperturbed) cells.
    The normalization proceeds by first normalizing the UMI counts within each cell to the median
    UMI count within the population. The control population is normalized to the same amount. 
    The expression within the population is then Z-normalized with respect to the control 
    distribution: i.e., for each gene the control mean is subtracted, and then these values are
    divided by the control standard deviation.
        
    Args:
        pop: CellPopulation to normalize
        control_cells: metadata condition passed to pop.where to identify control cell population
            to normalize with respect to
        scale_by_total: Rescale UMI counts within each cell to population median (default: True)
        median_umi_count: If provided, set the median to normalize to. This is useful for example
        **kwargs: all other arguments are passed to pop.groupby, and hence to pop.where. These can
            be used to do more refined slicing of the population if desired.
            when normalizing multiple independent lanes/gemgroups.
    
    Returns:
        DataFrame of normalized expression data
        
    Example:
        >>>pop.normalized_matrix = normalize_expression_to_control(pop, 'perturbed == "control"')
    """
    matrix = pop.where(**kwargs)
    control_matrix = pop.where(cells=control_cells, **kwargs)
    return normalize_matrix_to_control(matrix, control_matrix, scale_by_total=scale_by_total, median_umi_count=median_umi_count)

def normalize_to_gemgroup_control(pop, control_cells, median_umi_count=None, **kwargs):
    """Normalizes a multi-lane 10x experiment. Cell within each gemgroup are normalized to the 
    control cells within the same gemgroup.
        
    Args:
        pop: CellPopulation to normalize
        control_cells: metadata condition passed to pop.where to identify control cell population
            to normalize with respect to
        median_umi_count: Value to normalize UMI counts to across lanes. If None (the default)
            then all cells are normalized to the median UMI count of control cells within the
            whole experiment.
        **kwargs: all other arguments are passed to pop.groupby, and hence to pop.where. These can
            be used to do more refined slicing of the population if desired.
    All other arguments are passed to normalize_matrix_to_control
    
    Returns:
        DataFrame of normalized expression data    
    
    Example:
        normalized_matrix = normalize_to_gemgroup_control(pop,
                                                          control_cells='guide_identity == "control"')
        will produce a normalized expression matrix where cells in each gemgroup are Z-normalized 
        with respect to the expression distribution of cells in that lane bearing the guide_identity
        "control".
    """
    
    gemgroup_iterator = izip(pop.groupby('gem_group', **kwargs),
                        pop.groupby('gem_group', cells=control_cells, **kwargs))
    
    gem_group_matrices = dict()

    if median_umi_count is None:
        median_umi_count = pop.where(cells=control_cells).sum(axis=1).median()
    print('Normalizing all cells to {0} UMI...'.format(median_umi_count))
    
    for (i, gemgroup_pop), (_, gemgroup_control_pop) in gemgroup_iterator:
        print('Processing gem group {0}'.format(i))
        t = time()
        gem_group_matrices[i] = normalize_matrix_to_control(gemgroup_pop,
                                                            gemgroup_control_pop,
                                                            median_umi_count=median_umi_count)
        print(time() - t)
    print('Merging submatrices...')
    return pd.concat(gem_group_matrices.values()).loc[pop.matrix.index]

def normalize_matrix_by_key(pop, key):
    subpop_matrices = dict()

    for name, matrix in tqdm_notebook(pop.groupby(key)):
        print(name)
        subpop_matrices[name] = equalize_UMI_counts(matrix)
    print('Merging submatrices...')
    return pd.concat(subpop_matrices.values()).loc[pop.matrix.index]

def inherit_normalized_matrix(pop, parent_pop):
    """Subset a parent population's normalized expression matrix to the cells within a given population
    """
    pop.normalized_matrix = parent_pop.normalized_matrix.loc[pop.matrix.index, pop.matrix.columns]