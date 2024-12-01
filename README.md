# SIG07

### <ins>Goal<ins>
* Test 10x GEM-X flex chemistry for use with Signal-seq
* Test feasibility of combinatorial studies in a plate-based array format.
* Refine analyses of combinatorial ligand data

### <ins>Experimental Design<ins>
**Cells**: TCRb+CD4+CD44+CD62L- naive CD4 T cells were sorted from 8W old B6 mice straight from JAX and kept in autoclaved cages on DND. Mice were harvested and split into 3 biological replicates, 1F and 1M mice per group. M and F mice were all taken from the same cage. Cells were isolated from skin-draining LNs, mesenteric LNs, and spleen. Cells were activated for 24h with 2.5 ug/mL plate-bound anti-CD3, 3 ug/mL soluble anti-CD28, and 30 U/mL hIL2. After 24h, cells were lifted from plate, washed once, and rested in cRMPI for 12h.

**Treatment**: Cells were treated in 384-well plates at 10,000 cells per well in 500 ul cRPMI. Cells were spinfected with virus (final conc. ~1.4X) 2000g x 30m. Viral combinations were randomly distributed along the plate. Viral supernatent was removed, cells were resuspended in 50 uL cRPMI, and placed into TC incubator for 5h. Cells were pooled, stained with Zombie NIR, and fixed per GEM-X Flex protocol 19h x 4C.

**Sequencing**: Fixed samples were sorted for live cells (meant to sort for virus high cells but set the wrong gates) and kept at 4C x 3 days. Samples were submitted to SAIL for hybridization and GEM processing.

**Caveats in Hindsight**: 
* Due to the way I did my combinatorial additions (column-wise addition followed by row-wise addition using viruses with the same barcode), each barcode pairing is essentially from two indistinguishable wells. For example, A2 and B1 would have the same barcode pair and would be indistinguishable.
* If I had sorted for virus high cells, around 70% of cells were virus high.

### <ins>Analysis<ins>
#### <ins>Determine barcode recovery efficiency with Flex chemistry<ins>
The workflow of Flex is easier than with universal GEM-X because you can use fixed cells. In addition, it is *significantly cheaper* making it an attractive option for large-scale perturbation screens. However, I need to verify that the barcodes are sensitively and specifically detected in this assay.
##### Files:

#### <ins>Examine variability between biological replicates<ins>
While I can't deconvolute exact well identitiy (see caveats), I can at least examine the correlation between biological replicates across plates.
##### Files: 

#### <ins>Differential Transcriptome Analysis<ins>
Examine transcriptional similarity between ligands and look at differentially expressed genes.
##### Files:
* correlation_analysis.ipynb: clustering and inter-ligand comparisons using heatmaps and MDE
* deg_analysis_SIG07.ipynb: running and testing different DEG approaches

#### <ins>Linear modeling to recover additive vs emergent interaction effects.<ins>
Using gene-level and overall expression values, we can quantify the interaction effects across different ligand pairings.
##### Files:
* interaction_glm_notebook.ipynb: testing different glm approaches and coding
* interaction_glmGamPoi.r: R script to run GamPoi GLM in parallel format on cluster
* interaction_glmGamPoi_viz.Rmd: R markdown visualization of glmGamPoi results

#### <ins>Factor analysis<ins>
Perhaps could be used to identify robust downstream signaling networks that correspond to specific signaling intermediates? Can use spectra, topic modeling, etc.

#### <ins>Test usefulness of trajectory analysis in mapping combinatorial stimulations.<ins>
For a subset of interactions, it would be interesting to place them along a trajectory. There are many scenarios one can imagine (e.g. do all combinations lead to distinct branches? Are there terminal branches that are only reached with combinatorial additions?). It would also be useful if trajectory analysis with palantir could be used to identify gene modules that change with combinatorial additions.

#### 5. <ins>Test non-linear ICA model (Ola)<ins>
