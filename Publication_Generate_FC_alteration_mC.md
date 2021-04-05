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



```python
paths = [(connectome_p / connectome_t.format(row.ID)).resolve() for rid, row in pheno.iterrows()]
conn_stack = np.array([np.load(p)[conn_mask] for p in paths])

group = 'CNV'

regressors = ' + '.join(['age',  'C(sex)','FD_scrubbed', 'C(site_cohort)', 'mean']) 
```

#### Generate CNV-FC profiles

```python
glm_ASDcon = cnvfc.stats.glm_wrap_cc(conn_stack, pheno,  group, case="ASD", control="CON_idiop", regressors=regressors, report=True)
table_ASDcon,  table_stand_beta_ASDcon, table_qval_ASDcon = cnvfc.tools.summarize_glm(glm_ASDcon,conn_mask,roi_labels)
```

    Selected sample based on group variable CNV.
    cases: ASD (n=225)
    controls: CON (n=848)
    2080 data points available
    standardized estimators are based on CNV==CON_idiop



#### SAVE FC-PROFILES

```python
table_ASDcon.to_csv(out_p / 'cnv_ukbb_asd_vs_con_mean.tsv', sep='\t')
table_stand_beta_ASDcon.to_csv(out_p / 'cnv_ukbb_asd_vs_con_standardized_betas_mean.tsv', sep='\t')
table_qval_ASDcon.to_csv(out_p / 'cnv_ukbb_asd_vs_con_fdr_corrected_pvalues_mean.tsv', sep='\t')

```
