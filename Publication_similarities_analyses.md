#### Clara Moreau and Sebastian Urchs, March 2020
##### Contact: clara.moreau@umontreal.ca , sebastian.urchs@mail.mcgill.ca 


```python
import sys
sys.path.append('../../../')
import patsy as pat
import numpy as np
import pandas as pd
import pathlib as pal
import scipy as sp
from sklearn import linear_model as sln
import cnvfc
import seaborn as sbn
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests as stm
from matplotlib import pyplot as plt
```

### Regress individual connectomes


```python
root_p = pal.Path('../../../data/')
pheno_p = root_p / 'pheno/PHENO_all9cohorts_connectomes_meanGH.csv'
connectome_p = root_p / 'preprocessed/connectome/all9cohorts/python/'
connectome_t = 'connectome_{}_cambridge64.npy'
out_p = root_p / 'processed/residual_connectomes_publication/'
if not out_p.is_dir():
    out_p.mkdir()

conn_mask = np.tril(np.ones((64, 64))).astype(bool)
pheno = pd.read_csv(pheno_p, sep=';')

regressors = ['age', 'C(sex)', 'FD_scrubbed', 'C(site_cohort)',  'mean']


available_subject_mask = [True if (connectome_p / connectome_t.format(row.ID)).resolve().exists() else False 
                          for rid, row in pheno.iterrows()]

available_subjects = [str(i) for i in pheno.loc[available_subject_mask].ID.values]
pheno = pheno.loc[available_subject_mask]
paths = [(connectome_p / connectome_t.format(row.ID)).resolve() 
         for rid, row in pheno.iterrows()]
conn_stack = np.array([np.load(p)[conn_mask] for p in paths])
regressors_str = ' + '.join(regressors)
```


```python
# Control the connectomes for nuisance covariates
model = pat.dmatrix(regressors_str, data=pheno, return_type='dataframe')
# Run the glm
mod = sln.LinearRegression(fit_intercept=False, normalize=False)
res = mod.fit(model, conn_stack)
# Get the residuals
resid = conn_stack - res.predict(model)

np.save(out_p / 'icc_residual_connectomes_connectomes_mean.npy', resid)

# Store the list of subjects used to generate the connectomes
with (out_p / 'icc_subjects_included_in_residual_connectomes_mean.tsv').open('w') as f:
    f.write('\n'.join(available_subjects))
```

## Make weights for whole brain and regional similarity analyses

### 1) Whole-brain


```python
def make_flat_weights(data, pattern):
    if not data.shape[1:] == pattern.shape:
        raise Exception(f'data and pattern shape mismatch: {data.shape[1:]}, {pattern.shape}')
    n_nodes = pattern.shape[0]
    n_data = data.shape[0]
    weights = np.array([np.corrcoef(data[data_id, ...].flatten(), pattern.flatten())[0, 1]
                        for data_id in range(n_data)])
    return weights
```


```python
out_p = root_p / 'processed/weights_publication/'
if not out_p.is_dir():
    out_p.mkdir()
connectomes_p = root_p / 'processed/residual_connectomes_publication/icc_residual_connectomes_connectomes_mean.npy'
conn_mask = np.tril(np.ones((64, 64))).astype(bool)
```

### to run for each CNV of interest


```python
profile_p = root_p / 'processed/fc_profiles_publication/cnv_ukbb_del22q_vs_con_mean.tsv'
profile = pd.read_csv(profile_p, sep='\t')
profile_mat = cnvfc.tools.conn2mat(profile.betas.values, conn_mask)

connectomes = np.load(connectomes_p)
n_sub = connectomes.shape[0]
# Cast the vectorized connectomes back to matrices
connectome_mat = np.array([cnvfc.tools.conn2mat(connectomes[i, :], conn_mask)
                           for i in range(connectomes.shape[0])])

w = make_flat_weights(connectome_mat, profile_mat)

np.save(out_p / 'icc_del22q_flatweights_mean.npy', w)
```

### pLI


```python
profile_p = root_p / 'processed/fc_profiles_publication/icc_ukbb_pLIdel_mc.tsv'
profile = pd.read_csv(profile_p, sep='\t')
profile_mat = cnvfc.tools.conn2mat(profile.betas.values, conn_mask)

connectomes = np.load(connectomes_p)
n_sub = connectomes.shape[0]
# Cast the vectorized connectomes back to matrices
connectome_mat = np.array([cnvfc.tools.conn2mat(connectomes[i, :], conn_mask)
                           for i in range(connectomes.shape[0])])

w = make_flat_weights(connectome_mat, profile_mat)

np.save(out_p / 'icc_pLIdel_flatweights_mean.npy', w)
```

### 2) Regional similarity analyses


```python
profile_p = root_p / 'processed/fc_profiles_publication/cnv_ukbb_del22q_vs_con_mean.tsv'

connectomes_p = root_p / 'processed/residual_connectomes_publication/icc_residual_connectomes_connectomes_mean.npy'

conn_mask = np.tril(np.ones((64, 64))).astype(bool)
profile = pd.read_csv(profile_p, sep='\t')
profile_mat = cnvfc.tools.conn2mat(profile.betas.values, conn_mask)

connectomes = np.load(connectomes_p)
n_sub = connectomes.shape[0]
# Cast the vectorized connectomes back to matrices
connectome_mat = np.array([cnvfc.tools.conn2mat(connectomes[i, :], conn_mask)
                           for i in range(connectomes.shape[0])])

w = cnvfc.stats.make_weights(connectome_mat, profile_mat)
np.save(out_p / 'icc_del22q_weights.npy', w)

w.shape
```




    (6643, 64)




```python
profile_p = root_p / 'processed/fc_profiles_publication/icc_ukbb_pLIdel_mc.tsv'

connectomes_p = root_p / 'processed/residual_connectomes_publication/icc_residual_connectomes_connectomes_mean.npy'

conn_mask = np.tril(np.ones((64, 64))).astype(bool)
profile = pd.read_csv(profile_p, sep='\t')
profile_mat = cnvfc.tools.conn2mat(profile.betas.values, conn_mask)

connectomes = np.load(connectomes_p)
n_sub = connectomes.shape[0]
# Cast the vectorized connectomes back to matrices
connectome_mat = np.array([cnvfc.tools.conn2mat(connectomes[i, :], conn_mask)
                           for i in range(connectomes.shape[0])])

w = cnvfc.stats.make_weights(connectome_mat, profile_mat)
np.save(out_p / 'icc_plidel_weights.npy', w)

w.shape
```




    (6643, 64)



## Get Results

#### Define functions


```python
def find_subset(pheno, column, cases=None):
    # TODO check pandas type input
    subset_mask = np.array(~pheno[column].isnull())
    if cases is not None and not not cases:
        all_cases = pheno.loc[subset_mask][column].unique()
        try:
            case_available = np.array([True if case in all_cases else False for case in cases])
        except TypeError as e:
            raise Exception(f'the attribute "cases" needs to be iterable but is: {type(cases)}') from e
        if not all(case_available):
            if not any(case_available):
                raise Exception(f'none of the requested cases of "{column}" are available')
            else:
                warnings.warn(
                    f'\nnot all requested cases of "{column}" are available: {list(zip(cases, case_available))}',
                    RuntimeWarning)
        case_masks = np.array([pheno[column] == case for case in cases])
        subset_mask = np.any(case_masks, 0)
        # Return the masked instances of the requested cases
        cases_dict = {case: case_masks[idx][subset_mask] for idx, case in enumerate(cases)}
        return subset_mask, cases_dict
    else:
        return subset_mask

def categorical_enrichment_flat(pheno, data, group, case=None, control=None):
    sub_mask, cases_mask = find_subset(pheno, group, [case, control])
    sub_data = data[sub_mask, ...]
    n_data = sub_data.shape[0]
    n_case = np.sum(cases_mask[case])
    n_control = np.sum(cases_mask[control])   
    results = {key: list() for key in ['U', 'p_value',
                                       f'median_{group}_{case}',
                                       f'median_{group}_{control}',
                                       'rank_biserial_correlation']}
    u_right, p = sp.stats.mannwhitneyu(sub_data[cases_mask[case]], sub_data[cases_mask[control]],
    alternative='two-sided')
    u_left = n_case * n_control - u_right
    u_min = np.min([u_left, u_right])
    # Compute the median for the case and control groups
    median_case = np.median(sub_data[cases_mask[case]])
    median_control = np.median(sub_data[cases_mask[control] ])
    # Compute rank biserial correlation
    r_b = 1 - (2 * u_min) / (n_case * n_control)
    # Determine if cases > controls or reverse
    case_gt_con = u_right > u_min
    if not case_gt_con:
        r_b = -r_b
    # Store the results
    results['U'].append(u_min)
    results['p_value'].append(p)
    results[f'median_{group}_{case}'].append(median_case)
    results[f'median_{group}_{control}'].append(median_control)
    results['rank_biserial_correlation'].append(r_b)
    results_table = pd.DataFrame(data=results)
    # Correct for multiple comparisons
    (fdr_pass, qval, _, _) = stm(results_table.p_value, alpha=0.05, method='fdr_bh')
    results_table['q_value'] = qval
    return results_table
```

### Whole Brain


```python
root_p = pal.Path('../../../')
subjects_p = root_p / 'data/processed/residual_connectomes_publication/icc_subjects_included_in_residual_connectomes_mean.tsv'
flatweights_p = root_p / 'data/processed/weights_publication/icc_del22q_flatweights_mean.npy'
flatweights = np.load(flatweights_p)

pheno_p = root_p / 'data/pheno/PHENO_all9cohorts_connectomes_meanGH.csv'
pheno = pd.read_csv(pheno_p, sep=';')
out_p = root_p / 'data/tables_publication/'
    
with subjects_p.open('r') as f:
    subjects = f.read().splitlines()
# Limit the dataset to the subjects used in generating the weights
sub_mask = [str(row.ID) in subjects for rid, row in pheno.iterrows()]
pheno = pheno.loc[sub_mask]
```


```python
#ASD
diag_enrich = categorical_enrichment_flat(pheno, flatweights, group='CNV', case='ASD', control='CON_idiop')
print(diag_enrich)
diag_enrich.to_csv(out_p / 'icc_del22qWB_diagnosis_enrichment_in_ASD.tsv', sep='\t', index=False)
```

             U   p_value  median_CNV_ASD  median_CNV_CON_idiop  \
    0  91700.0  0.000915        0.014683              -0.01066   
    
       rank_biserial_correlation   q_value  
    0                   0.141988  0.000915  



```python
#SZ
diag_enrich = categorical_enrichment_flat(pheno, flatweights, group='CNV', case='SZ', control='CON_idiop')
print(diag_enrich)
diag_enrich.to_csv(out_p / 'icc_del22qWB_diagnosis_enrichment_in_SZ.tsv', sep='\t', index=False)
```

             U   p_value  median_CNV_SZ  median_CNV_CON_idiop  \
    0  96384.0  0.000103       0.017531              -0.01066   
    
       rank_biserial_correlation   q_value  
    0                   0.161514  0.000103  



```python
#ADHD
diag_enrich = categorical_enrichment_flat(pheno, flatweights, group='CNV', case='ADHD', control='CON_idiop')
print(diag_enrich)
diag_enrich.to_csv(out_p / 'icc_del22qWB_diagnosis_enrichment_in_ADHD.tsv', sep='\t', index=False)
```

              U   p_value  median_CNV_ADHD  median_CNV_CON_idiop  \
    0  124802.0  0.019199         0.011219              -0.01066   
    
       rank_biserial_correlation   q_value  
    0                   0.090861  0.019199  


### Regional similarities


```python
subjects_p = root_p / 'data/processed/residual_connectomes_publication/icc_subjects_included_in_residual_connectomes_mean.tsv'
label_p = root_p / 'data/parcellation/Parcel_Information/MIST_64.csv'
labels = pd.read_csv(label_p, sep=';')
roi_labels = labels.name.values
alt_label_p = root_p / 'data/MIST_64_short_labels.csv'
with alt_label_p.open('r') as f:
    alt_roi_labels = f.read().splitlines()[1:]
```


```python
# to edit 
weights_p = root_p / 'data/processed/weights_publication/icc_del22q_weights.npy'
weights = np.load(weights_p)
```


```python
diag_enrich = cnvfc.stats.categorical_enrichment(pheno, weights, group='CNV', case='ASD', control='CON_idiop', node_labels=alt_roi_labels)
print(cnvfc.tools.report_enrichment(diag_enrich))
diag_enrich.to_csv(out_p / 'icc_del22qRB_diagnosis_enrichment_in_ASD.tsv', sep='\t', index=False)
```

    Corrected Statistics:
    There are 6 out of 64 regions that show significant enrichment.
    FC in putamen is positively enriched:
        U: 88188 q_value: 9.49E-04 rank_corr: 0.175.
    FC in thalamus is positively enriched:
        U: 84494 q_value: 6.46E-05 rank_corr: 0.209.
    FC in temporal pole is positively enriched:
        U: 88868 q_value: 1.34E-03 rank_corr: 0.168.
    FC in post. cingulate ctx d. is positively enriched:
        U: 87889 q_value: 9.49E-04 rank_corr: 0.178.
    FC in perigenual ant. cingulate ctx is positively enriched:
        U: 90511 q_value: 4.48E-03 rank_corr: 0.153.
    FC in auditory net. is positively enriched:
        U: 91218 q_value: 6.66E-03 rank_corr: 0.146.


### Relation with severity scores


```python
# for continous measure
def continuous_enrichment_flat(pheno, data, covariate):
    sub_mask = find_subset(pheno, covariate)
    sub_pheno = pheno.loc[sub_mask]
    sub_data = data[sub_mask]
    results = {key: list() for key in [ 'pearson_r', 'p_value']}
    r, p = stats.pearsonr(sub_pheno[covariate], sub_data)
    results['pearson_r'].append(r)
    results['p_value'].append(p)
    results_table = pd.DataFrame(data=results)
    # Correct for multiple comparisons
    (fdr_pass, qval, _, _) = stm(results_table.p_value, alpha=0.05, method='fdr_bh')
    results_table['q_value'] = qval
    return results_table
```


```python
#to edit
weights_p = root_p / 'data/processed/weights_publication/icc_plidel_weights.npy'
weights = np.load(weights_p)
```


```python
fiq_enrich = cnvfc.stats.continuous_enrichment(pheno, weights, 'fluid_intelligence_score_f20016_2_0', alt_roi_labels)
#print(cnvfc.tools.report_enrichment(fiq_enrich, corrected=True))
fiq_enrich.to_csv(out_p / 'icc_FSIQ_enrichment_pLIdel_mC.tsv', sep='\t', index=False)
```


```python

```
