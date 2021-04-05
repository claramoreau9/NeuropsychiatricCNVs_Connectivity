#### Clara Moreau and Sebastian Urchs, October 2020
##### Contact: clara.moreau@umontreal.ca 


```python
%matplotlib inline
```


```python
import sys
sys.path.append('../../../')
import cnvfc
import numpy as np
import pandas as pd
import pathlib as pal
import numpy as np
import pandas as pd
import pathlib as pal
import seaborn as sbn
from matplotlib import pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests as stm
```


```python
root_p = pal.Path('../../../data/')

pheno_p = root_p / 'pheno/PHENO_connectomes.csv'

connectome_p = root_p / 'preprocessed/connectome/wholedataset/python/'

connectome_t = 'connectome_{}_cambridge64.npy'
label_p = root_p / 'parcellation/Parcel_Information/MIST_64.csv'
out_p = root_p / 'processed/fc_profiles_publication/'
if not out_p.is_dir():
    out_p.mkdir()
```


```python
conn_mask = np.tril(np.ones((64,64))).astype(bool)
pheno = pd.read_csv(pheno_p, sep=';')
labels = pd.read_csv(label_p, sep=';')
roi_labels = labels.label.values
```

## One CNV at a time


```python
paths = [(connectome_p / connectome_t.format(row.ID)).resolve() for rid, row in pheno.iterrows()]
conn_stack = np.array([np.load(p)[conn_mask] for p in paths])

group = 'CNV'

regressors = ' + '.join(['age',  'C(sex)','FD_scrubbed', 'C(site_cohort)', 'mean']) 
```

#### Generate CNV-FC profiles


```python
glm_del1qcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="DEL1q21_1", control="CON", regressors=regressors, report=True)
table_del1qcon,  table_stand_beta_del1qcon, table_qval_del1qcon = cnvfc.tools.summarize_glm(glm_del1qcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DEL1q21_1 (n=25)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4452
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_dup1qcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="DUP1q21_1", control="CON", regressors=regressors, report=True)
table_dup1qcon,  table_stand_beta_dup1qcon, table_qval_dup1qcon = cnvfc.tools.summarize_glm(glm_dup1qcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DUP1q21_1 (n=16)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4443
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_dup1qTARcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="TAR", control="CON", regressors=regressors, report=True)
table_dup1qTARcon,  table_stand_beta_dup1qTARcon, table_qval_dup1qTARcon = cnvfc.tools.summarize_glm(glm_dup1qTARcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: TAR (n=18)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4445
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_del15qcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="DEL15q11_2", control="CON", regressors=regressors, report=True)
table_del15qcon,  table_stand_beta_del15qcon, table_qval_del15qcon = cnvfc.tools.summarize_glm(glm_del15qcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DEL15q11_2 (n=66)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4493
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_dup15qcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="DUP15q11_2", control="CON", regressors=regressors, report=True)
table_dup15qcon,  table_stand_beta_dup15qcon, table_qval_dup15qcon = cnvfc.tools.summarize_glm(glm_dup15qcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DUP15q11_2 (n=71)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4498
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_del_2q13con = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case='DEL2q13', control='CON', regressors=regressors, report=True)
table_del_2q13con,  table_stand_beta_del_2q13con, table_qval_del_2q13con, = cnvfc.tools.summarize_glm(glm_del_2q13con,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DEL2q13 (n=36)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4463
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_dup_2q13con = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case='DUP2q13', control='CON', regressors=regressors, report=True)
table_dup_2q13con,  table_stand_beta_dup_2q13con, table_qval_dup_2q13con, = cnvfc.tools.summarize_glm(glm_dup_2q13con,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DUP2q13 (n=30)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4457
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_dup15q13con = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="DUP15q13_3", control="CON", regressors=regressors, report=True)
table_dup15q13con,  table_stand_beta_dup15q13con, table_qval_dup15q13con = cnvfc.tools.summarize_glm(glm_dup15q13con,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DUP15q13_3 (n=40)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4467
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_del22qcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="DEL22q11_2", control="CON", regressors=regressors, report=True)
table_del22qcon,  table_stand_beta_del22qcon, table_qval_del22qcon = cnvfc.tools.summarize_glm(glm_del22qcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DEL22q11_2 (n=43)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4470
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_dup22qcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="DUP22q11_2", control="CON", regressors=regressors, report=True)
table_dup22qcon,  table_stand_beta_dup22qcon, table_qval_dup22qcon = cnvfc.tools.summarize_glm(glm_dup22qcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DUP22q11_2 (n=19)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4446
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_del16pcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="DEL16p11_2", control="CON", regressors=regressors, report=True)
table_del16pcon,  table_stand_beta_del16pcon, table_qval_del16pcon = cnvfc.tools.summarize_glm(glm_del16pcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DEL16p11_2 (n=41)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4468
    2080 data points available
    standardized estimators are based on CNV==CON



```python
glm_dup16pcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="DUP16p11_2", control="CON", regressors=regressors, report=True)
table_dup16pcon,  table_stand_beta_dup16pcon, table_qval_dup16pcon = cnvfc.tools.summarize_glm(glm_dup16pcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: DUP16p11_2 (n=31)
    controls: CON (n=4427)
    original sample: n=6643; new sample: n=4458
    2080 data points available
    standardized estimators are based on CNV==CON


#### SAVE FC-PROFILES


```python
#1q21.1 DEL
table_del1qcon.to_csv(out_p / 'cnv_ukbb_del1q_vs_con_mean.tsv', sep='\t')
table_stand_beta_del1qcon.to_csv(out_p / 'cnv_ukbb_del1q_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_del1qcon.to_csv(out_p / 'cnv_ukbb_del1q_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

#1q21.1 DUP
table_dup1qcon.to_csv(out_p / 'cnv_ukbb_dup1q_vs_con_mean.tsv', sep='\t')
table_stand_beta_dup1qcon.to_csv(out_p / 'cnv_ukbb_dup1q_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_dup1qcon.to_csv(out_p / 'cnv_ukbb_dup1q_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

#TAR
table_dup1qTARcon.to_csv(out_p / 'cnv_ukbb_dup1qTAR_vs_con_mean.tsv', sep='\t')
table_stand_beta_dup1qTARcon.to_csv(out_p / 'cnv_ukbb_dup1qTAR_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_dup1qTARcon.to_csv(out_p / 'cnv_ukbb_dup1qTAR_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')


#15q11 del
table_del15qcon.to_csv(out_p / 'cnv_ukbb_del15q_vs_con_mean.tsv', sep='\t')
table_stand_beta_del15qcon.to_csv(out_p / 'cnv_ukbb_del15q_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_del15qcon.to_csv(out_p / 'cnv_ukbb_del15q_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

#15q11 dup
table_dup15qcon.to_csv(out_p / 'cnv_ukbb_dup15q_vs_con_mean.tsv', sep='\t')
table_stand_beta_dup15qcon.to_csv(out_p / 'cnv_ukbb_dup15q_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_dup15qcon.to_csv(out_p / 'cnv_ukbb_dup15q_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

#15q13 dup
table_dup15q13con.to_csv(out_p / 'cnv_ukbb_dup15q13_vs_con_mean.tsv', sep='\t')
table_stand_beta_dup15q13con.to_csv(out_p / 'cnv_ukbb_dup15q13_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_dup15q13con.to_csv(out_p / 'cnv_ukbb_dup15q13_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

#2q13 del
table_del_2q13con.to_csv(out_p / 'cnv_ukbb_del2q13_vs_con_mean.tsv', sep='\t')
table_stand_beta_del_2q13con.to_csv(out_p / 'cnv_ukbb_del2q13_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_del_2q13con.to_csv(out_p / 'cnv_ukbb_del2q13_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

#2q13 dup
table_dup_2q13con.to_csv(out_p / 'cnv_ukbb_dup2q13_vs_con_mean.tsv', sep='\t')
table_stand_beta_dup_2q13con.to_csv(out_p / 'cnv_ukbb_dup2q13_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_dup_2q13con.to_csv(out_p / 'cnv_ukbb_dup2q13_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

#del 22q11.2
table_del22qcon.to_csv(out_p / 'cnv_ukbb_del22q_vs_con_mean.tsv', sep='\t')
table_stand_beta_del22qcon.to_csv(out_p / 'cnv_ukbb_del22q_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_del22qcon.to_csv(out_p / 'cnv_ukbb_del22q_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

# dup 22q11.2
table_dup22qcon.to_csv(out_p / 'cnv_ukbb_dup22q_vs_con_mean.tsv', sep='\t')
table_stand_beta_dup22qcon.to_csv(out_p / 'cnv_ukbb_dup22q_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_dup22qcon.to_csv(out_p / 'cnv_ukbb_dup22q_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

# del 16p11.2
table_del16pcon.to_csv(out_p / 'cnv_ukbb_del16p_vs_con_mean.tsv', sep='\t')
table_stand_beta_del16pcon.to_csv(out_p / 'cnv_ukbb_del16p_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_del16pcon.to_csv(out_p / 'cnv_ukbb_del16p_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

#dup 16p11.2
table_dup16pcon.to_csv(out_p / 'cnv_ukbb_dup16p_vs_con_mean.tsv', sep='\t')
table_stand_beta_dup16pcon.to_csv(out_p / 'cnv_ukbb_dup16p_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_dup16pcon.to_csv(out_p / 'cnv_ukbb_dup16p_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')
```

## Idiopathic FC-PROFILES


```python
glm_ASDcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="ASD", control="CON_idiop", regressors=regressors, report=True)
table_ASDcon,  table_stand_beta_ASDcon, table_qval_ASDcon = cnvfc.tools.summarize_glm(glm_ASDcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: ASD (n=225)
    controls: CON_idiop (n=950)
    original sample: n=6643; new sample: n=1175
    2080 data points available
    standardized estimators are based on CNV==CON_idiop



```python
glm_SZcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="SZ", control="CON_idiop", regressors=regressors, report=True)
table_SZcon,  table_stand_beta_SZcon, table_qval_SZcon = cnvfc.tools.summarize_glm(glm_SZcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: SZ (n=242)
    controls: CON_idiop (n=950)
    original sample: n=6643; new sample: n=1192
    2080 data points available
    standardized estimators are based on CNV==CON_idiop



```python
glm_ADHDcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="ADHD", control="CON_idiop", regressors=regressors, report=True)
table_ADHDcon,  table_stand_beta_ADHDcon, table_qval_ADHDcon = cnvfc.tools.summarize_glm(glm_ADHDcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: ADHD (n=289)
    controls: CON_idiop (n=950)
    original sample: n=6643; new sample: n=1239
    2080 data points available
    standardized estimators are based on CNV==CON_idiop


### SAVE IPC FC-Profiles


```python
# ASD
table_ASDcon.to_csv(out_p / 'cnv_ukbb_asd_vs_con_mean.tsv', sep='\t')
table_stand_beta_ASDcon.to_csv(out_p / 'cnv_ukbb_asd_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_ASDcon.to_csv(out_p / 'cnv_ukbb_asd_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

# SZ
table_SZcon.to_csv(out_p / 'cnv_ukbb_sz_vs_con_mean.tsv', sep='\t')
table_stand_beta_SZcon.to_csv(out_p / 'cnv_ukbb_sz_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_SZcon.to_csv(out_p / 'cnv_ukbb_sz_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

# ADHD
table_ADHDcon.to_csv(out_p / 'cnv_ukbb_adhd_vs_con_mean.tsv', sep='\t')
table_stand_beta_ADHDcon.to_csv(out_p / 'cnv_ukbb_adhd_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_ADHDcon.to_csv(out_p / 'cnv_ukbb_adhd_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')
```

## Beyond case-control, pLI-FC profiles


```python
pheno_p = root_p / 'pheno/pLI/PHENO_connectomes_pli.csv' 

pheno = pd.read_csv(pheno_p, sep=';')
pheno.rename(columns={'ID':'niak_id'}, inplace=True)
available_subject_mask = [True if (connectome_p / connectome_t.format(row.niak_id)).resolve().exists() else False 
                          for rid, row in pheno.iterrows()]
pheno = pheno.loc[available_subject_mask]
paths = [(connectome_p / connectome_t.format(row.niak_id)).resolve() for rid, row in pheno.iterrows()]
conn_stack = np.array([np.load(p)[conn_mask] for p in paths])
```

### General effect of Deletion on FC


```python
contrast = 'pli_del'

regressors = ' + '.join(['age', 'C(sex)', 'FD_scrubbed', 'C(site_cohort)', 'mean'])

glm_pLIdel = cnvfc.stats.glm_wrap_continuous(conn_stack, pheno, contrast, regressors, report=True, fast=False)
table_pLIdel, table_stand_beta_pLIdel, table_qval_pLIdel = cnvfc.tools.summarize_glm(glm_pLIdel, conn_mask, roi_labels)

# Store the results
table_pLIdel.to_csv(out_p / 'icc_ukbb_pLIdel_mc.tsv', sep='\t')
table_stand_beta_pLIdel.to_csv(out_p / 'icc_ukbb_pLIdel_mc_standardized_betas.tsv', sep='\t')
table_qval_pLIdel.to_csv(out_p / 'icc_ukbb_pLIdel_mc_fdr_corrected_pvalues.tsv', sep='\t')

##to generate chord diagram
betas_mat = cnvfc.tools.conn2mat(table_pLIdel.betas.values, conn_mask)
qvals_mat = cnvfc.tools.conn2mat(table_pLIdel.qval.values, conn_mask)

pd.DataFrame(betas_mat).to_csv(out_p / 'cnv_ukbb_table_pLIdel_mc_betas.tsv', sep='\t')
pd.DataFrame(qvals_mat).to_csv(out_p / 'cnv_ukbb_pLIdel_mc_qval.tsv', sep='\t')
```

    Selected sample based on contrast variable pli_del.
    Found 4664 subjects with no missing data for pli_del
    original sample: n=4664; new sample: n=4664
    2080 data points available
    standardized estimators are based on all subjects with no missing data for pli_del


### General effect of Duplication on FC


```python
contrast = 'pli_dup'

regressors = ' + '.join(['age', 'C(sex)', 'FD_scrubbed', 'C(site_cohort)', 'mean'])


glm_pLIdup = cnvfc.stats.glm_wrap_continuous(conn_stack, pheno, contrast, regressors, report=True, fast=False)
table_pLIdup, table_stand_beta_pLIdup, table_qval_pLIdup = cnvfc.tools.summarize_glm(glm_pLIdup, conn_mask, roi_labels)

# Store the results
table_pLIdup.to_csv(out_p / 'icc_ukbb_pLIdup_mc.tsv', sep='\t')
table_stand_beta_pLIdup.to_csv(out_p / 'icc_ukbb_pLIdup_mc_standardized_betas.tsv', sep='\t')
table_qval_pLIdup.to_csv(out_p / 'icc_ukbb_pLIdup_mc_fdr_corrected_pvalues_mc.tsv', sep='\t')

##to generate chord diagram
betas_mat = cnvfc.tools.conn2mat(table_pLIdup.betas.values, conn_mask)
qvals_mat = cnvfc.tools.conn2mat(table_pLIdup.qval.values, conn_mask)
pd.DataFrame(betas_mat).to_csv(out_p / 'cnv_ukbb_table_pLIdup_mc_betas.tsv', sep='\t')
pd.DataFrame(qvals_mat).to_csv(out_p / 'cnv_ukbb_pLIdup_mc_qval.tsv', sep='\t')
```

    Selected sample based on contrast variable pli_dup.
    Found 4692 subjects with no missing data for pli_dup
    original sample: n=4692; new sample: n=4692
    2080 data points available
    standardized estimators are based on all subjects with no missing data for pli_dup



```python

```
