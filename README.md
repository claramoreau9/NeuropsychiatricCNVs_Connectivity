## Genetic heterogeneity shapes brain connectivity in psychiatry

#### Authors: Clara A Moreau, Annabelle Harvey, Kumar Kuldeep,...
clara.moreau@umontreal.ca

#### Last authors: Bellec Pierre and Jacquemont Sebastien 

#### Notebooks (still work in progess)

[Interactive brainmaps](https://claramoreau9.github.io/Braimaps_Github.html) To navigate into the 29 CWAS results (beta values reported on brain maps for each FC-Profiles, before and after global signal adjustement) using [nilearn](https://www.frontiersin.org/articles/10.3389/fninf.2014.00014/full)</p>
[Generate Functional Connectivity profile (.md)](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/Publication_Generate_FC_alteration_mC.md) This notebook requests [stats.py](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/stats.py) and [tools.py](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/tools.py) </p>
[Get CWAS Results (.md)](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/Publication_Results_FC_alteration_mC_brainmap.md) This notebook needs [stats.py](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/stats.py) and [tools.py](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/tools.py) </p>

#### Results
[FC-profiles : 29 CWAS](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/tree/master/results_tables) This folder contains all the beta values / p-values / q-values (n=2080 connections) for the 29 FC-profiles.</p>
[Correlation values - gene expression](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/df_CorrPerGene_PvalLabelShuffle_PvalBrainSmash_MIST64_THAL.xlsxs) This excel file contains correlation values for Fig6 - as well as p-values obtained with [BrainSMASH](https://www.sciencedirect.com/science/article/pii/S1053811920305243) or label-shuffling (permutation test by performing similar contrasts in 5000 randomly sampled groups).</p>

#### External repositories linked to this manuscript
[Permutation Tests and cross-validations by A Harvey](https://github.com/harveyaa/cross_cnv_paper_permutations) This repository contains all the scripts and notebooks for generating null models of betamaps (from this manuscript) through permutation and testing for significance.</p>



%%[Gene expression analyses by K. Kumar](https://github.com/kkumar-iitkgp-livia/GeneExp_and_CNV_FCsignatures)</p>

[Mind-GenesParallelCNV by JL Martineau](https://github.com/MartineauJeanLouis/MIND-GENESPARALLELCNV) to compute CNV calling parallel tasks in the most efficient method, in UKBiobank</p>

[Previous notebooks by C. Moreau and S. Urchs](https://github.com/surchs/Neuropsychiatric_CNV_code_supplement) </p>

#### Scripts
[to preprocess resting-state_fMRI_data with NIAK (.m)](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/to_preprocess_restingstate_fMRI_data.m) </p>
[to_generate_connectomes with NIAK using MIST64 (.m)](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/to_generate_connectomes.m) </p>
[to generate heatmaps for CellType and AHBA modules Fig 6 (.R)](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/Rscript_HeatMap_CellType_and_AHBAmodules.R) </p>
[to create chord diagrams Fig 5(.R)](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/to_create_chord_diag.r) </p>

#### MIST ATLAS 
[Multiscale dashboard](https://simexp.github.io/multiscale_dashboard/index.html)
