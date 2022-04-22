import pandas as pd
import os
os.chdir("/Aging_hallmark_ARD_repository-main/")

################################### EDA PubMed baseline ###################################

def read_csv_Append(number):
    """Loops through and selects the abstracts with relevant Pubmed IDs"""
    if (number >= 100):
        file_string = "pubmed_data/baseline_csv/pubmed19n0" + str(number) + ".csv"
    elif (number >= 10):
        file_string = "pubmed_data/baseline_csv/pubmed19n00" + str(number) + ".csv"
    else:
        file_string = "pubmed_data/baseline_csv/pubmed19n000" + str(number) + ".csv"
    PubmedDoc = pd.read_csv(file_string, dtype={'pmid': 'object'})
    return(PubmedDoc)


pubmed_baseline0_new = pd.DataFrame()
for q in list(range(1, 100)):
    pubmed_baseline0_new = pubmed_baseline0_new.append(read_csv_Append(q))
    print(pubmed_baseline0_new.shape)

pubmed_baseline1_new = pd.DataFrame()
for q in list(range(100, 200)):
    pubmed_baseline1_new = pubmed_baseline1_new.append(read_csv_Append(q))
    print(pubmed_baseline1_new.shape)

pubmed_baseline2_new = pd.DataFrame()
for q in list(range(200, 300)):
    pubmed_baseline2_new = pubmed_baseline2_new.append(read_csv_Append(q))
    print(pubmed_baseline2_new.shape)

pubmed_baseline3_new = pd.DataFrame()
for q in list(range(300, 400)):
    pubmed_baseline3_new = pubmed_baseline3_new.append(read_csv_Append(q))
    print(pubmed_baseline3_new.shape)

pubmed_baseline4_new = pd.DataFrame()
for q in list(range(400, 500)):
    pubmed_baseline4_new = pubmed_baseline4_new.append(read_csv_Append(q))
    print(pubmed_baseline4_new.shape)

pubmed_baseline5_new = pd.DataFrame()
for q in list(range(500, 600)):
    pubmed_baseline5_new = pubmed_baseline5_new.append(read_csv_Append(q))
    print(pubmed_baseline5_new.shape)

pubmed_baseline6_new = pd.DataFrame()
for q in list(range(600, 700)):
    pubmed_baseline6_new = pubmed_baseline6_new.append(read_csv_Append(q))
    print(pubmed_baseline6_new.shape)

pubmed_baseline7_new = pd.DataFrame()
for q in list(range(700, 800)):
    pubmed_baseline7_new = pubmed_baseline7_new.append(read_csv_Append(q))
    print(pubmed_baseline7_new.shape)

pubmed_baseline8_new = pd.DataFrame()
for q in list(range(800, 900)):
    pubmed_baseline8_new = pubmed_baseline8_new.append(read_csv_Append(q))
    print(pubmed_baseline8_new.shape)

pubmed_baseline9_new = pd.DataFrame()
for q in list(range(900, 973)):
    pubmed_baseline9_new = pubmed_baseline9_new.append(read_csv_Append(q))
    print(pubmed_baseline9_new.shape)

pubmed_baseline_All_new = pd.concat([pubmed_baseline0_new, pubmed_baseline1_new, pubmed_baseline2_new, 
                                     pubmed_baseline3_new, pubmed_baseline4_new, pubmed_baseline5_new, 
                                     pubmed_baseline6_new, pubmed_baseline7_new, pubmed_baseline8_new, 
                                     pubmed_baseline9_new])

pubmed_baseline_All_new = pubmed_baseline_All_new.drop(['Unnamed: 0'], axis=1)
print("There are ", pubmed_baseline_All_new.shape[0], "rows in the 2019 PubMed baseline")
