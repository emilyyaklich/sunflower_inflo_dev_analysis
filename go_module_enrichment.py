# Name: go module enrichment
# Author: EY
# Date: Sept. 20 2023
# Version: Python 3.10
# Description: Will calculate the network module enrichment using GOATOOLS

import pandas as pd
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
import pickle as pkl
from goatools import obo_parser
import os
import sys
from pandas import ExcelWriter
import glob

# read in the GO ontology
go_obo = 'go-basic.obo'
go = obo_parser.GODag(go_obo)

# read in the data for which genes belong to which modules in the network
# the format is 'column', 'mRNA', 'go_term', 'module'
#   # there is only column per GO term, so the same mRNA molecules can be in multiple lines
#go_module_data= pd.read_csv('/home/ely67071/sunflower_stress_analysis/wgcna/go_module_data.csv')
go_module_data= pd.read_csv('go_module_data.csv')

# turn the df into a list to parse through
go_module_data_list= go_module_data.values.tolist()

# create a nested dictionary where the key is the mrna molecule and the nested dictionaries are the GO terms
#   # and module color. Ex:
#   # 'mRNA:Ha412HOChr04g0148041': {'go_terms': ['GO:0005471', 'GO:0006862', 'GO:0016020'], 'module': ['brown']}, 'mRNA:Ha412HOChr13g0588161': {'go_terms': ['GO:0003676'], 'module': ['darkgreen']}
go_dict = {}
for i in go_module_data_list:
    go_dict[i[1]] = {}
    go_dict[i[1]]['go_terms'] = []
    go_dict[i[1]]['module'] = []
for j in go_module_data_list:
    if j[1] in go_dict.keys():
        go_dict[j[1]]['go_terms'].append(j[2])
        if j[3] not in go_dict[j[1]]['module']:
            go_dict[j[1]]['module'].append(j[3])

# create mrna_id_dict which is a dictionary where the key is the mrna and the values are a list of GO terms
#   # 'mRNA:Ha412HOChr04g0148041': ['GO:0005471', 'GO:0006862', 'GO:0016020']}
# also crease a unique term list with all of the GO terms in the entire dataset represented only once
# also create a list that contains module names
mrna_id_dict = {}
unique_term_list = []
modules=[]
# create the complex GO ID dict and the unique ID list
for mrna, terms in go_dict.items():
    mrna_id_dict[mrna] = terms['go_terms']
    unique_term_list.append(terms['go_terms'])
    # only add the module if it is not already in the name list
    if terms["module"][0] not in modules:
        modules.append(terms["module"][0])
# flatten the nested list and then get only unique terms
unique_term_list = [item for elem in unique_term_list for item in elem]
unique_term_list = list(set(unique_term_list))


# for each go term in our dataset, label whether it is CC, BP, or MF
CC = []
BP = []
MF = []
for go_id in unique_term_list:
    if go[go_id].namespace == 'cellular_component':
        CC.append(go_id)
    elif go[go_id].namespace == 'biological_process':
        BP.append(go_id)
    elif go[go_id].namespace == 'molecular_function':
        MF.append(go_id)

# create the reference dictionary that everything will be compared to
# for each GO namespace create a nested dictionary which contains the GO ids associated with each mRNA
#   # {'MF': {'mRNA:Ha412HOChr04g0148041': {'GO:0005471'}, 'mRNA:Ha412HOChr13g0588161': {'GO:0003676'}, 'mRNA:Ha412HOChr03g0141851': {'GO:0004869'}, 'mRNA:Ha412HOChr03g0142491': {'GO:0046983', 'GO:0003700'}
reference_dict = {}
reference_dict['MF'] = {}
reference_dict['BP'] = {}
reference_dict['CC'] = {}
for mrna, go_ids_associated_w_module in mrna_id_dict.items():
    x = []
    for go_id_mf in MF:
        if go_id_mf in go_ids_associated_w_module:
            x.append(go_id_mf)
    if len(x) > 0:
        reference_dict['MF'][mrna] = set(x)
    y = []
    for go_id_bp in BP:
        if go_id_bp in go_ids_associated_w_module:
            y.append(go_id_bp)
    if len(y) > 0:
        reference_dict['BP'][mrna] = set(y)
    z = []
    for go_id_cc in CC:
        if go_id_cc in go_ids_associated_w_module:
            z.append(go_id_cc)
    if len(z) > 0:
        reference_dict['CC'][mrna] = set(z)

# create the enrichment
goeaobj = GOEnrichmentStudyNS(
        mrna_id_dict.keys(),  # List of mrna molecules in the analysis
        reference_dict,  # mrna/GO associations
        go,  # Ontologies
        propagate_counts=True,  # this will include teh ancestors
        alpha=0.05,  # default significance cut-off
        methods=['fdr_bh'])

# for each module in the list of modules, run the enrichment analysis
for module in modules:
    # create a list of all the mrna molecules that belong to a specific module
    mrna_w_module = []
    for mrna in go_dict:
        if go_dict[mrna]["module"][0] == module:
            if mrna not in mrna_w_module:
                mrna_w_module.append(mrna)
    # run the enrichment analysis
    goea_results_all = goeaobj.run_study(mrna_w_module)
    # separate the significant results
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    # write the significant results to an excel file (will only do this if the results are significant)
    goeaobj.wr_xlsx(module+'_module_enrichment_results.xlsx', goea_results_sig)

# combine all of the spreadsheets into a single spreadsheet with different "sheets"
writer = ExcelWriter("module_significant_enrichment_results_inflo.xlsx")
with pd.ExcelWriter("module_significant_enrichment_results_inflo.xlsx") as writer:
    for filename in glob.glob("*_module_enrichment_results.xlsx"):
        excel_file = pd.ExcelFile(filename)
        (_, f_name) = os.path.split(filename)
        res = f_name.split("_", 2)
        f_short_name= res[0]
        for sheet_name in excel_file.sheet_names:
            df_excel = pd.read_excel(filename, sheet_name=sheet_name)
            df_excel.to_excel(writer, f_short_name, index=False)

# remove individual spreadsheets
for filename in glob.glob("*_module_enrichment_results.xlsx"):
    os.remove(filename)

