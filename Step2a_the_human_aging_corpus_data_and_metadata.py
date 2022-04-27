import pandas as pd
import re
import unidecode
import string
import os
from nltk import sent_tokenize
from nltk.tokenize.treebank import TreebankWordDetokenizer
from nltk.tokenize import RegexpTokenizer
print(string.punctuation)
string.punctuation = '!"#$%&\'*+/:;<=>?@[\\]^_`{}~'
tokenizer = RegexpTokenizer(r'\w+')
import nltk
nltk.download('punkt')
os.chdir("/Aging_hallmark_ARD_repository-main/")

################################### Deriving the titles, abstracts and metadata for the human aging corpus ###################################

def read_csv_Append(number, corpus):
    """Loops through and selects the abstracts with relevant Pubmed IDs for the human aging corpus"""
    if (number >= 100):
        file_string = "pubmed_data/baseline_csv/pubmed19n0" + str(number) + ".csv"
    elif (number >= 10):
        file_string = "pubmed_data/baseline_csv/pubmed19n00" + str(number) + ".csv"
    else:
        file_string = "pubmed_data/baseline_csv/pubmed19n000" + str(number) + ".csv"
    PubmedDoc = pd.read_csv(file_string, dtype={'pmid': 'object'})
    abstractsDf = pd.DataFrame.merge(PubmedDoc, corpus, on="pmid")
    return(abstractsDf)


def abstract_split_clean(pubmed_abstracts_frame, a, b):
    """Removes abstracts without PMID, article title or abstract, appends title to the abstract, returns a new dataframe
    NB. The titles have punctuation at the end (e.g., full stop, question marks)"""
    pubmed_abstracts_frame = pubmed_abstracts_frame.iloc[a:b, 0:5]
    pubmed_abstracts_frame = pubmed_abstracts_frame.dropna(subset=['pmid'])
    pubmed_abstracts_frame = pubmed_abstracts_frame.dropna(subset=['articletitle'])
    pubmed_new = pubmed_abstracts_frame.dropna(subset=['abstract'])
    title_Abstract = pubmed_new.articletitle + " " + pubmed_new.abstract
    new_sentence_dataframe = pd.DataFrame({'pmid_all': pubmed_new['pmid'], 'unedited_sent': title_Abstract,
                                           'mesh_all': pubmed_new['meshHeadings']})
    return(new_sentence_dataframe)


human_aging_corpus1 = pd.read_csv("data/human_aging_corpus/human_aging_corpus1.csv", dtype={'pmid': 'object'})
human_aging_corpus2 = pd.read_csv("data/human_aging_corpus/human_aging_corpus2.csv", dtype={'pmid': 'object'})
human_aging_corpus = human_aging_corpus1.append(human_aging_corpus2)
human_aging_corpus = human_aging_corpus.dropna(subset=["pmid"])
human_aging_corpus = human_aging_corpus.drop_duplicates(subset=["pmid"])
print("There are", human_aging_corpus.pmid.nunique(), "unique PubMed IDs in the human aging corpus")
GWAScat_pmid = pd.read_csv("data/human_aging_corpus/PubMedID_list.csv")
GWAScat_pmid.columns = ['index', 'pmid']
GWAScat_pmid.pmid = GWAScat_pmid.pmid.astype(str)
human_aging_corpus = human_aging_corpus[~human_aging_corpus['pmid'].isin(GWAScat_pmid['pmid'])]
human_aging_corpus = human_aging_corpus.drop_duplicates(subset=['pmid'])
print("There are", human_aging_corpus.pmid.nunique(), "unique PubMed IDs in the human aging corpus")


ah_corpus_abstracts0 = pd.DataFrame()
for q in list(range(1, 100)):
    ah_corpus_abstracts0 = ah_corpus_abstracts0.append(read_csv_Append(q, human_aging_corpus))
    print(ah_corpus_abstracts0.shape)
ah_corpus_abstracts0.to_csv("data/human_aging_corpus/ah_corpus_abstracts_1_to_100.csv")

ah_corpus_abstracts1 = pd.DataFrame()
for q in list(range(100, 200)):
    ah_corpus_abstracts1 = ah_corpus_abstracts1.append(read_csv_Append(q, human_aging_corpus))
    print(ah_corpus_abstracts1.shape)
ah_corpus_abstracts1.to_csv("data/human_aging_corpus/ah_corpus_abstracts_100_to_200.csv")

ah_corpus_abstracts2 = pd.DataFrame()
for q in list(range(200, 300)):
    ah_corpus_abstracts2 = ah_corpus_abstracts2.append(read_csv_Append(q, human_aging_corpus))
    print(ah_corpus_abstracts2.shape)
ah_corpus_abstracts2.to_csv("data/human_aging_corpus/ah_corpus_abstracts_200_to_300.csv")

ah_corpus_abstracts3 = pd.DataFrame()
for q in list(range(300, 400)):
    ah_corpus_abstracts3 = ah_corpus_abstracts3.append(read_csv_Append(q, human_aging_corpus))
    print(ah_corpus_abstracts3.shape)
ah_corpus_abstracts3.to_csv("data/human_aging_corpus/ah_corpus_abstracts_300_to_400.csv")

ah_corpus_abstracts4 = pd.DataFrame()
for q in list(range(400, 500)):
    ah_corpus_abstracts4 = ah_corpus_abstracts4.append(read_csv_Append(q, human_aging_corpus))
    print(ah_corpus_abstracts4.shape)
ah_corpus_abstracts4.to_csv("data/human_aging_corpus/ah_corpus_abstracts_400_to_500.csv")

ah_corpus_abstracts5 = pd.DataFrame()
for q in list(range(500, 600)):
    ah_corpus_abstracts5 = ah_corpus_abstracts5.append(read_csv_Append(q, human_aging_corpus))
    print(ah_corpus_abstracts5.shape)
ah_corpus_abstracts5.to_csv("data/human_aging_corpus/ah_corpus_abstracts_500_to_600.csv")

ah_corpus_abstracts6 = pd.DataFrame()
for q in list(range(600, 700)):
    ah_corpus_abstracts6 = ah_corpus_abstracts6.append(read_csv_Append(q, human_aging_corpus))
    print(ah_corpus_abstracts6.shape)
ah_corpus_abstracts6.to_csv("data/human_aging_corpus/ah_corpus_abstracts_600_to_700.csv")

ah_corpus_abstracts7 = pd.DataFrame()
for q in list(range(700, 800)):
    ah_corpus_abstracts7 = ah_corpus_abstracts7.append(read_csv_Append(q, human_aging_corpus))
    print(ah_corpus_abstracts7.shape)
ah_corpus_abstracts7.to_csv("data/human_aging_corpus/ah_corpus_abstracts_700_to_800.csv")

ah_corpus_abstracts8 = pd.DataFrame()
for q in list(range(800, 900)):
    ah_corpus_abstracts8 = ah_corpus_abstracts8.append(read_csv_Append(q, human_aging_corpus))
    print(ah_corpus_abstracts8.shape)
ah_corpus_abstracts8.to_csv("data/human_aging_corpus/ah_corpus_abstracts_800_to_900.csv")

ah_corpus_abstracts9 = pd.DataFrame()
for q in list(range(900, 973)):
    ah_corpus_abstracts9 = ah_corpus_abstracts9.append(read_csv_Append(q, human_aging_corpus))
    print(ah_corpus_abstracts9.shape)
ah_corpus_abstracts9.to_csv("data/human_aging_corpus/ah_corpus_abstracts_900_to_972.csv")

ah_corpus_abstracts0 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_1_to_100.csv", dtype={'year': 'object', 'nctID': 'object', 'pmid': 'object'})
ah_corpus_abstracts1 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_100_to_200.csv", dtype={'year': 'object', 'nctID': 'object', 'pmid': 'object'})
ah_corpus_abstracts2 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_200_to_300.csv", dtype={'year': 'object', 'nctID': 'object', 'pmid': 'object'})
ah_corpus_abstracts3 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_300_to_400.csv", dtype={'year': 'object', 'nctID': 'object', 'pmid': 'object'})
ah_corpus_abstracts4 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_400_to_500.csv", dtype={'year': 'object', 'nctID': 'object', 'pmid': 'object'})
ah_corpus_abstracts5 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_500_to_600.csv", dtype={'year': 'object', 'nctID': 'object', 'pmid': 'object'})
ah_corpus_abstracts6 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_600_to_700.csv", dtype={'year': 'object', 'nctID': 'object', 'pmid': 'object'})
ah_corpus_abstracts7 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_700_to_800.csv", dtype={'year': 'object', 'nctID': 'object', 'pmid': 'object'})
ah_corpus_abstracts8 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_800_to_900.csv", dtype={'year': 'object', 'nctID': 'object', 'pmid': 'object'})
ah_corpus_abstracts9 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_900_to_972.csv", dtype={'year': 'object', 'nctID': 'object', 'pmid': 'object'})
ah_corpus_abstracts_All = pd.concat([ah_corpus_abstracts0, ah_corpus_abstracts1, ah_corpus_abstracts2, ah_corpus_abstracts3, ah_corpus_abstracts4, ah_corpus_abstracts5, ah_corpus_abstracts6, ah_corpus_abstracts7, ah_corpus_abstracts8, ah_corpus_abstracts9])
del ah_corpus_abstracts0, ah_corpus_abstracts1, ah_corpus_abstracts2, ah_corpus_abstracts3, ah_corpus_abstracts4, ah_corpus_abstracts5
del ah_corpus_abstracts6, ah_corpus_abstracts7, ah_corpus_abstracts8, ah_corpus_abstracts9

ah_corpus_abstracts_All = ah_corpus_abstracts_All.drop_duplicates(subset=["pmid"])
ah_corpus_abstracts_All = ah_corpus_abstracts_All.dropna(subset=["pmid"])

human_aging_corpus.pmid.nunique() - ah_corpus_abstracts_All.pmid.nunique() 
abstracts_missing = set(human_aging_corpus.pmid) - set(ah_corpus_abstracts_All.pmid)

ah_corpus_abstracts_All1 = ah_corpus_abstracts_All.iloc[0:1000000, ]
ah_corpus_abstracts_All2 = ah_corpus_abstracts_All.iloc[1000000:, ]
ah_corpus_abstracts_All1.to_csv("data/human_aging_corpus/ah_corpus_abstracts_All1.csv")
ah_corpus_abstracts_All2.to_csv("data/human_aging_corpus/ah_corpus_abstracts_All2.csv")

ah_corpus_abstracts_All1 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_All1.csv", dtype={'year': 'object', 'pmid': 'object'})
ah_corpus_abstracts_All2 = pd.read_csv("data/human_aging_corpus/ah_corpus_abstracts_All2.csv", dtype={'year': 'object', 'pmid': 'object'})

ah_corpus_abstracts_All = ah_corpus_abstracts_All1.append(ah_corpus_abstracts_All2)
ah_corpus_abstracts_All = ah_corpus_abstracts_All.drop(['Unnamed: 0', 'Unnamed: 0_x', 'Unnamed: 0_y'], axis=1)
ah_corpus_abstracts_All_metadata = ah_corpus_abstracts_All[['pmid', 'year', 'meshHeadings', 'ptype']]
ah_corpus_abstracts_All_metadata.to_csv("data/human_aging_corpus/ah_corpus_abstracts_All_metadata.csv")
ah_corpus_abstracts_All = ah_corpus_abstracts_All.drop(['year'], axis=1)

abstract_frame1 = abstract_split_clean(ah_corpus_abstracts_All, 0, 1000000)
abstract_frame1.to_csv("data/human_aging_corpus/ah_sentences_frame_0_1000000.csv")
abstract_frame2 = abstract_split_clean(ah_corpus_abstracts_All, 1000000, 2000000)
abstract_frame2.to_csv("data/human_aging_corpus/ah_sentences_frame_1000000_1902967.csv")

abstract_frame1 = pd.read_csv("data/human_aging_corpus/ah_sentences_frame_0_1000000.csv")
abstract_frame2 = pd.read_csv("data/human_aging_corpus/ah_sentences_frame_1000000_1902967.csv")
abstract_frames = abstract_frame1.append(abstract_frame2, sort=False)
del abstract_frame1, abstract_frame2

abstract = abstract_frames.unedited_sent
abstract = [re.sub('\|', ' ', i) for i in abstract[0:]]
abstract = [re.sub("'", '', i) for i in abstract[0:]]
abstract = [sent_tokenize(i) for i in abstract[0:]]
pmid = list(pd.Series(abstract_frames.pmid_all[0:]))
new_frame = pd.DataFrame({'pmid_all': pmid, 'abstracts': abstract})
abstract_series = new_frame.apply(lambda x: pd.Series(x['abstracts']), axis=1).stack().reset_index(level=1, drop=True)
abstract_series.name = 'abstracts'

abstract_series2 = new_frame.drop('abstracts', axis=1).join(abstract_series)
abstract_series2['abstracts_decode'] = [unidecode.unidecode(i) for i in abstract_series2.abstracts[0:]]
abstract_series2['abstracts_decode'] = [tokenizer.tokenize(i) for i in abstract_series2['abstracts_decode'][0:]]
abstract_series2['abstracts_decode'] = [TreebankWordDetokenizer().detokenize(i) for i in abstract_series2['abstracts_decode'][0:]]
abstract_series2['abstracts_decode'] = [i.lower() for i in abstract_series2['abstracts_decode'][0:]]
abstract_series2 = abstract_series2.drop_duplicates()
abstract_series2['len'] = [len(x) for x in abstract_series2.abstracts[0:]]
abstract_series2 = abstract_series2[abstract_series2['len'] > 5]
abstract_series2 = abstract_series2.drop(['len'], axis=1)
abstract_series2.to_csv("data/human_aging_corpus/abstracts.csv")
abstract_series2 = pd.read_csv("data/human_aging_corpus/abstracts.csv")
print("There are", abstract_series2.pmid_all.nunique(), "unique PubMed IDs in the human aging corpus after processing")
print("There are", abstract_series2.abstracts_decode.count(), "unique sentences in the human aging corpus after processing")
