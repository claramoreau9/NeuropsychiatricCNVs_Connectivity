## The general impact of haploinsufficiency on brain connectivity underlies the pleiotropic effect of neuropsychiatric CNVs

[Preprint on medRXiv](https://www.medrxiv.org/content/10.1101/2020.03.18.20038505v1) </p>

#### Authors: Moreau Clara, Huguet Guillaume, Urchs Sebastian, and colleagues
clara.moreau@umontreal.ca
sebastian.urchs@mail.mcgill.ca 

#### Last authors: Bellec Pierre, and Jacquemont Sebastien 

### One sentence summary: Neuropsychiatric CNVs across the genome reorganize brain connectivity architecture along dominant patterns contributing to complex idiopathic conditions.

### Abstract
<p align="justify"> Polygenicity and pleiotropy are key genomic features of psychiatric disorders, but how they impact large-scale brain connectivity has not been investigated. 
 </p>
<p align="justify"> We analyzed resting-state fMRI data from 32,988 individuals to characterize the effects on functional connectivity of 16 multigenic copy number variants (CNVs), 10 polygenic scores, 6 cognitive and brain morphometry traits as well as 4 psychiatric conditions. </p>
<p align="justify"> Effect-sizes of CNVs on connectivity were correlated with their impact on cognition but this rapidly tapered off as CNVs increased in the number of encompassed genes. Accordingly, polygenic scores had 6-fold lower effect-sizes on connectivity compared to CNVs. Pleiotropy measured by genetic and transcriptomic correlations across 38 pairs of conditions and traits showed stablesignificant concordance with connectomic correlations observed between the same conditions and traits, suggesting a substantial causal genetic component for shared connectivity.
Despite substantial specificities, different forms of molecular risk and psychiatric conditions converged on the overconnectivity of the thalamus and somatomotor network.  </p>

### Notebooks and scripts </p>
*More codes are going to be released in May 2021* </p>

#### Aim 1 Characterize the impact of gene dosage on connectivity for 16 multigenic copy number variants, 10 polygenic scores, 6 cognitive and brain morphometry traits, 4 psychiatric conditions.

Connectome Wide Association Study: Linear model contrasting CNV carriers with controls, adjusting models for mean connectivity </p>
##### Scripts
[to_preprocess_restingstate_fMRI_data.m](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/to_preprocess_restingstate_fMRI_data.m) </p>
[to_generate_connectomes.m](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/to_generate_connectomes.m) </p>
[Publication_Generate_FC_alteration_mC.md](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/Publication_Generate_FC_alteration_mC.md)</p>
[Publication_Results_FC_alteration_mC_brainmap.md](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/Publication_Results_FC_alteration_mC_brainmap.md)</p>


#### Aim 2 Test similarities between whole-brain connectomes of CNVs and idiopathic psychiatric conditions
Pearson correlation between beta maps obtained from each CWAS and individual connectomes of either ASD, SZ, ADHD cases or controls 
##### Scripts
[Publication_similarities_analyses.md](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/Publication_similarities_analyses.md) </p>

#### Aim 3 General effect of haploinsufficiency on FC and latent component shared across CNVs and IPCs
Linear model testing the effect of pLI scores on FC in the CNVs carriers and non-carriers population </p>
Correlation between similarity to pLI-FC-profile and IQ, ADOS, SRS measures of each individual </p>
##### Scripts
[Publication_Results_FC_alteration_mC_brainmap.md](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/Publication_Results_FC_alteration_mC_brainmap.md)</p>
[Publication_similarities_analyses.md](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/Publication_similarities_analyses.md) </p>
[to_create_chord_diag.r](https://github.com/claramoreau9/NeuropsychiatricCNVs_Connectivity/blob/master/to_create_chord_diag.r) </p>
