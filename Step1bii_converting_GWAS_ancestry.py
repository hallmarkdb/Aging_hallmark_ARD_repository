from itertools import groupby
import pandas as pd
import os
os.chdir("/Aging_hallmark_ARD_repository-main/")

#################################### ARD to Gene Search terms ############################################# 
groups = []
def function_find(string2):
    """Returns the unique upper case words related to ancestry in the GWASCat_Categories.both column 
    (e.g., European, Hispanic, African)"""
    groups2 = []
    for key, group in groupby(string2.split(), lambda x: x[0].isupper()):
        if key:
           groups2.append(' '.join(group))
    groups.append(list(set(groups2)))
    return()


GWASCat_Categories = pd.read_csv("Genetic_data/GWAS_Ancestry_Europe.csv")
GWASCat_Categories_copy = [str(i) for i in GWASCat_Categories.both[0:]]
a = [function_find(i) for i in GWASCat_Categories_copy[0:]]
groups3 = [str(i) for i in groups[0:]]
GWASCat_Categories_update = pd.DataFrame({
    'Initial': GWASCat_Categories.Initial,
    'Replication': GWASCat_Categories.Replication,
    'Groups': groups3,
    'Study Accession': GWASCat_Categories.Accesion,
    'Broad': GWASCat_Categories.Broad
    })
GWASCat_Categories_update = GWASCat_Categories_update.drop_duplicates()
GWASCat_Categories_update.to_csv("Genetic_data/GWASCat_Categories_update.csv")
