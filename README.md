# Reafect
Reafect uses Z-scores obtained from (untargeted) metabolomics data to calculate 'deficient reaction scores'. These scores indicate whether it is more or less likely that a considered reaction is deficient.

Reafect uses KEGG modules and pathways as basis, but reactions can manually be added. 

# Content

**Databases\\**
* kegg_enzyme_with_gene_symbol.csv contains (KEGG) associations of enzymes with genes.
* kegg_enzyme_with_reaction.csv contains (KEGG) associations of enzymes with reactions.
* kegg_enzyme_with_name.csv contains (KEGG) EC enzyme numbers with enzyme names. 
* Z_scores_IEM_patients_Bongaerts.csv: Z-score data of 72 IEM patient samples used in the publication.
* Z_scores_IEM_patients_Miller.csv: Z-score data from Miller et al. (https://doi.org/10.1007/s10545-015-9843-7). Some adjustments were made to prepare the dataset for Reafect. 

**Reafect\\**
Contains all the code used to run Reafect.

*Reafect_reproduce_results_publication.ipynb* reproduces IEM ranking results as presented in our publication


## Manual for running Reafect

First, import the Reafect module.
```python
from Reafect import Reafect
ReafectObj = Reafect(path_to_store_maps = 'Maps/',
              path_to_abc_E = 'Abc_E/',
              path_to_store_results = 'Results/',
              path_to_databases ='Databases/')
```

Reafect uses metabolic modules and pathways from KEGG. These pathways/modules need to be downloaded first.

```python
hsa_modules = ['M00001'] # etc
ReafectObj.get_modules_kegg(hsa_modules)
               
hsa_maps = ['00010'] # etc
ReafectObj.get_pathways_kegg(hsa_maps)
```

An investigator can manually add reactions to all pathways/modules as shown in *Reafect_reproduce_results_publication.ipynb*. After downloading all KEGG modules/pathways, we need to determine the effective Z-scores by searching (reaction) paths between metabolites and reactions in each module/pathway. These effective Z-scores are stored as an symbolic expression (using Sympy) such that this calculation is only done once. 

```python
ReafectObj.prepare_pathways()
```
After determining all effective Z-scores as an symbolic expression, we can substitute the hyper parameters *a*, *b*, *c* for real values. The following method replaces all *a*, *b*, *c* for real values.

```python
ReafectObj.determine_pathways_for_abc_params(a_params= [0.85], b_params = [0.35] ,c_params = [0.75])
```
In order to determine the deficient reaction scores for a sample, the input needs to have the following format (Pandas dataframe):

| network_ID   | Z_score  |
|--------------|------------|
| C00026       | 2.2        |
| C00037       | -0.5       |
| ..           | ..         |
| RM_mn        | 1.25       |

This dataframe can be passed to the *process_patient()* method.

```python
SR_scores = ReafectObj.process_patient( patient_ID = 'S01', 
                                        Z_score_data = Z_scores_per_patient, 
                                        a_param= 0.85, b_param = 0.35, c_param = 0.75 # <-- make sure you ran determine_pathways_for_abc_params() first for these values
                                       )
```
The results are stored in 'Abc_E/'.