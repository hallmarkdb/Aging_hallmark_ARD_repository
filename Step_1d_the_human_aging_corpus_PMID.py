import pandas as pd
import time
from Bio import Entrez
import os

################################### Retrieving human aging corpus PMIDs ###################################

os.chdir("Code_Thesis_HCF/Aging_hallmark_ARD_repository/")

# Create an api_key for your email and then assign the email and api_key to the variables below as strings. 
#my_email = 
#my_api_key = 


def entrez_search(query, count, email=my_email):
    """Retrieves a list of all the PubMed IDs retrieved for a particular search returning a maximum of 100,000"""
    time.sleep(0.1)
    print("Searching...")
    Entrez.email = email
    handle = Entrez.esearch(db='pubmed', sort='relevance', retstart=count, retmax=100000, idtype='acc', term=query, api_key=my_api_key)
    results1 = Entrez.read(handle)
    handle.close()
    results = results1["IdList"]
    return results


def search(query, count=0):
    """Retrieves a list of all the PubMed IDs retrieved for a particular search looping through each 100,000 returned PubMed IDs"""
    count = 0
    total_length = 100000
    results = entrez_search(query, count)
    results_new = results
    while (total_length - len(results)) == 0:
        count += 100000
        total_length = 100000
        results = entrez_search(query, count)
        results_new.extend(results)
    return results_new


# Defining the Human Aging Corpus: have the same output as best match on PubMed
inclusions_main = set(search('"Humans"[MeSH] AND ("1993/12/2"[Date - Publication] : "2017/12/31"[Date - Publication]) AND ("English"[Language] OR "English Abstract"[Publication Type]) AND "journal article"[Publication Type]'))
inclusions_aging = set(search('"aging"[TIAB] OR "ageing"[TIAB] OR "associated with increasing age"[TIAB] OR "age correlated"[TIAB] OR "age dependence"[TIAB] OR "age dependent"[TIAB]  OR "age induced"[TIAB] OR "age related"[TIAB] OR "age associated"[TIAB] OR "age specific"[TIAB] OR "advanced age"[TIAB] OR "aged"[TIAB] OR "health span"[TIAB] OR "health spans"[TIAB] OR "healthspan"[TIAB] OR "healthspans"[TIAB] OR "life span"[TIAB] OR "life spans"[TIAB] OR "lifespans"[TIAB] OR "lifespan"[TIAB] OR "life extension"[TIAB] OR "life extending"[TIAB] OR "length of life"[TIAB] OR "life length"[TIAB] OR "longevity"[TIAB] OR "long lived"[TIAB] OR "long living"[TIAB] OR "senescence"[TIAB] OR "senescent"[TIAB] OR "geriatric"[TIAB] OR "geriatrics"[TIAB] OR "gerontological"[TIAB] OR "senium"[TIAB] OR "senility"[TIAB] OR "senile"[TIAB] OR "gerontology"[TIAB] OR "biogerontology"[TIAB] OR "elderly"[TIAB] OR "old"[TIAB] OR "aged"[MeSH Terms] OR "aging"[MeSH Terms] OR "older"[TIAB] OR "senior citizen"[TIAB]'))
exclusions_article = set(search('"Review"[Publication Type] OR "Introductory Journal Article"[Publication Type] OR "Retraction of Publication"[Publication Type] OR "Letter"[Publication Type] OR "Comment"[Publication Type]'))
exclusions_germline = set(search('"germline mutation"[TIAB] OR "germline mutations"[TIAB] OR "germ line mutation"[TIAB] OR "germ line mutations"[TIAB] OR "hereditary mutation"[TIAB] OR "hereditary mutations"[TIAB] OR "germline alteration"[TIAB] OR "germline alterations"[TIAB] OR "germ line alteration"[TIAB] OR "germ line alterations"[TIAB] OR "hereditary alteration"[TIAB] OR "hereditary alterations"[TIAB] OR "polymorphism"[TIAB] OR "polymorphisms"[TIAB] OR "congenital"[TIAB] OR "familial"[TIAB] OR "family"[TIAB] OR "offspring"[TIAB] OR "inherited"[TIAB] OR "inheritance"[TIAB] OR "congenital"[TIAB] OR "germ-line mutation"[MeSH Terms] OR "polymorphism, single nucleotide"[MeSH Terms] OR "congenital"[Subheading]'))
exlusions_telomerase = set(search('"telomerase activation"[TIAB] OR "TERT activation"[TIAB] OR "HEST2 activation"[TIAB] OR "TSC1 activation"[TIAB] OR "telomerase associated protein 2 activation"[TIAB] OR "telomere reverse transcriptase activation"[TIAB] OR "TP2 activation"[TIAB] OR "activity of telomerase"[TIAB] OR "activity of TERT"[TIAB] OR "activity of telomerase reverse transcriptase"[TIAB] OR "activity of HEST2"[TIAB] OR "activity of TSC1"[TIAB] OR "activity of telomerase associated protein 2"[TIAB] OR "activity of telomere reverse transcriptase"[TIAB] OR "activity of TP2"[TIAB] OR "telomerase activity"[TIAB] OR "TERT activity"[TIAB]  OR "telomerase reverse transcriptase activity"[TIAB] OR "HEST2 activity"[TIAB] OR "TSC1 activity"[TIAB] OR "telomerase associated protein 2 activity"[TIAB] OR "telomere reverse transcriptase activity"[TIAB] OR "TP2 activity"[TIAB] OR "activation of telomerase"[TIAB] OR "activation of TERT"[TIAB] OR "activation of Telomerase Reverse Transcriptase"[TIAB] OR "activation of HEST2"[TIAB] OR "activation of TSC1"[TIAB] OR "activation of telomerase associated protein 2"[TIAB] OR "activation of telomere reverse transcriptase"[TIAB] OR "activation of TP2"[TIAB] OR "telomerase inhibition"[TIAB] OR "TERT inhibition"[TIAB] OR "Telomerase Reverse Transcriptase inhibition"[TIAB] OR "HEST2 inhibition"[TIAB] OR "TSC1 inhibition"[TIAB] OR "telomerase associated protein 2 inhibition"[TIAB] OR "telomere reverse transcriptase inhibition"[TIAB] OR "TP2 inhibition"[TIAB] OR "inhibition of telomerase"[TIAB] OR "inhibition of TERT"[TIAB] OR "inhibition of Telomerase Reverse Transcriptase"[TIAB] OR "inhibition of HEST2"[TIAB] OR "inhibition of TSC1"[TIAB] OR "inhibition of telomerase associated protein 2"[TIAB] OR "inhibition of telomere reverse transcriptase"[TIAB] OR "inhibition of TP2"[TIAB] OR "inhibitor of telomerase"[TIAB] OR "inhibitor of TERT"[TIAB] OR "inhibitor of Telomerase Reverse Transcriptase"[TIAB] OR "inhibitor of HEST2"[TIAB] OR "inhibitor of TSC1"[TIAB] OR "inhibitor of telomerase associated protein 2"[TIAB] OR "inhibitor of telomere reverse transcriptase"[TIAB] OR "inhibitor of TP2"[TIAB] OR "telomerase inhibitor"[TIAB] OR "TERT inhibitor"[TIAB] OR "Telomerase Reverse Transcriptase inhibitor"[TIAB] OR "HEST2 inhibitor"[TIAB] OR "TSC1 inhibitor"[TIAB] OR "telomerase associated protein 2 inhibitor"[TIAB] OR "telomere reverse transcriptase inhibitor"[TIAB] OR "TP2 inhibitor"[TIAB]'))
exclusions_drugs = set(search('"adverse effects"[Subheading] OR "adverse effect"[TIAB] OR "adverse effects"[TIAB]'))
exclusions_sctransplant = set(search('"stem cell transplant"[TIAB] OR "stem cell transplants"[TIAB] OR "stem cell transplantation"[TIAB] OR "stem cell therapies"[TIAB] OR "stem cell therapy"[TIAB] OR "stem cell treatment"[TIAB] OR "stem cell transplantation"[MeSH Terms]'))
exclusions_pluripotent = set(search('"induced pluripotent stem cells"[MeSH Terms] OR "IPS cell"[TIAB] OR "IPS cells"[TIAB] OR "induced pluripotent stem cells"[TIAB] OR "induced pluripotent stem cell"[TIAB] OR "iPSC"[TIAB]'))
exclusions_neoplasticsc = set(search('"Neoplastic Stem Cells"[MeSH Terms] OR "neoplastic stem cells"[TIAB] OR "cell stem tumors"[TIAB] OR "stem cell tumor"[TIAB] OR "tumor stem cell"[TIAB] OR "neoplastic stem cell"[TIAB] OR "tumor initiating cell"[TIAB] OR "tumor stem cell"[TIAB] OR "tumor initiating cells"[TIAB] OR "tumor stem cells"[TIAB] OR "cell stem tumours"[TIAB] OR "stem cell tumour"[TIAB] OR "tumour stem cell"[TIAB] OR "tumour initiating cell"[TIAB] OR "tumour stem cell"[TIAB] OR "tumour initiating cells"[TIAB] OR "tumour stem cells"[TIAB]'))
exclusions_hormoneprep = set(search('"hormone preparation"[TIAB] OR "hormone preparations"[TIAB] OR "hormone product"[TIAB] OR "hormone products"[TIAB] OR "product containing hormone"[TIAB] OR "product containing hormones"[TIAB] OR "drug hormones"[TIAB] OR "drugs hormone"[TIAB] OR "drugs hormones"[TIAB] OR "hormone medication"[TIAB] OR "hormone medications"[TIAB] OR "hormone drug"[TIAB] OR "hormone drugs"[TIAB]'))
exclusions_antiinflamm = set(search('"anti inflammatory agent"[TIAB] OR "antiinflammatory agent"[TIAB] OR "anti inflammatory agents"[TIAB] OR "antiinflammatory agents"[TIAB] OR "antiinflammatories agents"[TIAB] OR "anti inflammatory drug"[TIAB] OR "antiinflammatory drug"[TIAB] OR "anti inflammatory preparation"[TIAB] OR "antiinflammatory preparation"[TIAB] OR "anti inflammatory preparations"[TIAB] OR "antiinflammatory preparations"[TIAB] OR "anti inflammatory drugs"[TIAB]'))

# Deriving a list of included articles and excluded articles
human_aging_corpus_inclusions = inclusions_main.intersection(inclusions_aging)
exclusions_article_AH_corpus = list(set(human_aging_corpus_inclusions).intersection(set(exclusions_article)))
exclusions_germline_AH_corpus = list(set(human_aging_corpus_inclusions).intersection(set(exclusions_germline)))
exclusions_telomerase_AH_corpus = list(set(human_aging_corpus_inclusions).intersection(set(exlusions_telomerase)))
exclusions_drugs_AH_corpus = list(set(human_aging_corpus_inclusions).intersection(set(exclusions_drugs)))
exclusions_sctransplant_AH_corpus = list(set(human_aging_corpus_inclusions).intersection(set(exclusions_sctransplant)))
exclusions_neoplasticsc_AH_corpus = list(set(human_aging_corpus_inclusions).intersection(set(exclusions_neoplasticsc)))
exclusions_pluripotent_AH_corpus = list(set(human_aging_corpus_inclusions).intersection(set(exclusions_pluripotent)))
exclusions_hormoneprep_AH_corpus = list(set(human_aging_corpus_inclusions).intersection(set(exclusions_hormoneprep)))
exclusions_antiinflamm_AH_corpus = list(set(human_aging_corpus_inclusions).intersection(set(exclusions_antiinflamm)))
exclusions_hallmark_corpus = set(exclusions_article).union(set(exclusions_germline), set(exlusions_telomerase), set(exclusions_drugs), set(exclusions_sctransplant), set(exclusions_pluripotent), set(exclusions_neoplasticsc), set(exclusions_hormoneprep), set(exclusions_antiinflamm))
exclusions_hallmark_corpus = list(exclusions_hallmark_corpus)
exclusions_hallmark_corpus_intersect = set(exclusions_hallmark_corpus).intersection(set(human_aging_corpus_inclusions))

# Defining the PMIDs of the "human aging corpus" and saving to .csv file
human_aging_corpus = list(set(human_aging_corpus_inclusions) - set(exclusions_hallmark_corpus))
human_aging_corpus1 = human_aging_corpus[0:1000000]
human_aging_corpus2 = human_aging_corpus[1000000:]
human_aging_corpus1 = pd.DataFrame({'pmid': human_aging_corpus1})
human_aging_corpus2 = pd.DataFrame({'pmid': human_aging_corpus2})
#human_aging_corpus1.to_csv("data/human_aging_corpus/human_aging_corpus1.csv")
#human_aging_corpus2.to_csv("data/human_aging_corpus/human_aging_corpus2.csv")
print("The number of abstracts in the human aging corpus is: ", len(set(human_aging_corpus)))
