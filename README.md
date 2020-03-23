## The general impact of haploinsufficiency on brain connectivity underlies the pleiotropic effect of neuropsychiatric CNVs

[Preprint on medRXiv](https://www.medrxiv.org/content/10.1101/2020.03.18.20038505v1) </p>

#### Authors: Moreau Clara, Huguet Guillaume, Urchs Sebastian, and colleagues
clara.moreau@umontreal.ca
sebastian.urchs@mail.mcgill.ca 

#### Last authors: Bellec Pierre, and Jacquemont Sebastien 

### One sentence summary: Neuropsychiatric CNVs across the genome reorganize brain connectivity architecture along dominant patterns contributing to complex idiopathic conditions.

### Abstract
<p align="justify"> Large effect-size mutations such as copy number variants (CNVs) have the potential to provide key insights into the underlying biological mechanisms linking deleterious genetic variants to brain architecture and neuropsychiatric disorders. To date, the effect of CNVs on functional brain connectivity remains mostly unstudied, and findings are derived from analyses conducted one mutation at a time. The lack of systematic cross-CNVs comparisons hinders our understanding of any potential general mechanisms linking CNVs to effects on brain organization. </p>
<p align="justify"> We performed connectome-wide analyses using resting-state functional MRI data from 436 carriers of neuropsychiatric CNVs at the 1q21.1, 15q11.2, 16p11.2, 22q11.2 loci, 4 “neutral effect” CNVs, 66 carriers of scarcer neuropsychiatric CNVs, 756 individuals with idiopathic autism spectrum disorder (ASD), schizophrenia, attention deficit hyperactivity disorder, and 5,377 controls.
Neuropsychiatric CNVs showed global shifts of mean connectivity. The effect size of CNVs on relative connectivity (adjusted for the mean) was correlated with the known level of neuropsychiatric risk conferred by CNVs. Individuals with idiopathic schizophrenia and ASD had similarities in connectivity with neuropsychiatric CNVs. We reported a linear relationship between connectivity and intolerance to haploinsufficiency measured for all genes encompassed by CNVs across 18 loci. This profile involved the thalamus, the basal ganglia, somatomotor and frontoparietal networks and was correlated with lower general intelligence and higher autism severity scores. An exploratory factor analysis confirmed the contribution of these regions to three latent components shared across CNVs and neuropsychiatric disorders. </p>
<p align="justify"> We posit that deleting genes intolerant to haploinsufficiency reorganize connectivity along general dimensions irrespective of where deletions occur in the genome. This haploinsufficiency brain signature opens new avenues to understand polygenicity in psychiatric conditions and the pleiotropic effect of CNVs on cognition and risk for neuropsychiatric disorders. </p>

### Notebooks and scripts </p>
*More codes are going to be released in April 2020* </p>

#### Aim 1 Characterize the impact of gene dosage on connectivity for CNVs at eight genomic loci

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
