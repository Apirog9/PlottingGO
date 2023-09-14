# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 16:53:52 2023

@author: APirog
"""
import json
import pandas as pd
import requests, sys


def get_by_API(term_id):
    
    requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" +term_id
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if r.ok:
        print(r)
        data = json.loads(r.text)
        name = data["results"][0]["name"]
    else:
        print(r)
        r.raise_for_status()
        sys.exit()
        name = None
        
    return name
    
    


def create_tsv(gafpath,mapping_dict):
    '''
    Parameters
    ----------
    gafpath : path to gaf file


    Returns
    -------
    path to processed gaf file to tsv format
    json object with gene to BP mapping
    json object with gene to MF mapping
    json object with gene to CC mapping
    json object with protein to BP mapping
    json object with protein to MF mapping
    json object with protein to CC mapping

    '''
    
    def getmap(key,dictionary):
        try:
            value = dictionary[key]
        except KeyError:
                value = None  
        return value
    
    def check_dict(go_nums,mapping_dict):
        current = list(mapping_dict.keys())
        novel_ids = [x for x in go_nums if not x in current]
        for go_term in novel_ids:
            newname = get_by_API(go_term)
            if newname:
                mapping_dict[go_term] = get_by_API(go_term)
                
        return mapping_dict
        
        
    
    #write tsv and load it back as dataframe
    oldfile = open(gafpath,'r').read().splitlines()
    newfile = gafpath.split('/')[-1].removesuffix('.gaf')+'.tsv'
    with open(newfile,'w') as output:
        output.write('Type\tIdentifier\tGene\tGO_function\tGO_term\tcol6\tcol7\tcol8\tGO_category\tcol10\tcol11\tmolecule_type\tcol13\tcol14\tcol15\tcol16\tcol17\n')
        for line in oldfile:
            if line.startswith('!'):
                pass
            else:
               output.write(line+'\n')
                
    ontology_dataframe = pd.read_csv(newfile, sep='\t')
    # remove non-protein hits, thay need special treatment
    ontology_dataframe = ontology_dataframe[ontology_dataframe['molecule_type'] == 'protein']
    # remove negative mappings
    ontology_dataframe['GO_function'] = ontology_dataframe['GO_function'].fillna('')
    ontology_dataframe = ontology_dataframe[~ontology_dataframe['GO_function'].str.contains('NOT')]

    #translate to names    
    ontology_ids = list(ontology_dataframe['GO_term'].unique())
    mapping_dict = json.loads(mapping_dict)
    mapping_dict = check_dict(ontology_ids,mapping_dict)
    ontology_dataframe['GO_term_name'] =  ontology_dataframe['GO_term'].apply(getmap,args=[mapping_dict])
    
    #split to term types
    frame_BP = ontology_dataframe[ontology_dataframe['GO_category'] == 'P']
    frame_MF = ontology_dataframe[ontology_dataframe['GO_category'] == 'F']
    frame_CC = ontology_dataframe[ontology_dataframe['GO_category'] == 'C']
    
    #generate gene and uniprot id lists
    genes = ontology_dataframe['Gene'].unique()
    genes = [x for x in genes if isinstance(x,str)]
    proteins = ontology_dataframe['Identifier'].unique()
    proteins = [x for x in proteins if isinstance(x,str)]
    
    #prepare dictionaries and transform to jsons
    gene_BP_dict = {}
    for gene in genes:
        dataslice = frame_BP[frame_BP['Gene']==gene]
        terms = list(dataslice['GO_term_name'])
        gene_BP_dict[gene] = list(set(terms))
    gene_BP_dict = json.dumps(gene_BP_dict)
    
    gene_MF_dict = {}
    for gene in genes:
        dataslice = frame_MF[frame_MF['Gene']==gene]
        terms = list(dataslice['GO_term_name'])
        gene_MF_dict[gene] = list(set(terms))
    gene_MF_dict = json.dumps(gene_MF_dict)
    
    gene_CC_dict = {}
    for gene in genes:
        dataslice = frame_CC[frame_CC['Gene']==gene]
        terms = list(dataslice['GO_term_name'])
        gene_CC_dict[gene] = list(set(terms))
    gene_CC_dict = json.dumps(gene_CC_dict)
    
    protein_BP_dict = {}
    for protein in proteins:
        dataslice = frame_BP[frame_BP['Identifier']==protein]
        terms = list(dataslice['GO_term_name'])
        protein_BP_dict[protein] = list(set(terms))
    protein_BP_dict = json.dumps(protein_BP_dict)
    
    protein_MF_dict = {}
    for protein in proteins:
        dataslice = frame_MF[frame_MF['Identifier']==protein]
        terms = list(dataslice['GO_term_name'])
        protein_MF_dict[protein] = list(set(terms))
    protein_MF_dict = json.dumps(protein_MF_dict)
    
    protein_CC_dict = {}
    for protein in proteins:
        dataslice = frame_CC[frame_CC['Identifier']==protein]
        terms = list(dataslice['GO_term_name'])
        protein_CC_dict[protein] = list(set(terms))
    protein_CC_dict = json.dumps(protein_CC_dict)
    
    mapping_dict = json.dumps(mapping_dict)
    return newfile,gene_BP_dict,gene_MF_dict,gene_CC_dict,protein_BP_dict,protein_MF_dict,protein_CC_dict,mapping_dict



def get_namedict(name_dict_path):
    mapping = pd.read_csv(name_dict_path, sep = ' ')
    GO_identifiers = mapping["identifier"].unique()
    mapping_dict = {}
    for identifier in GO_identifiers:
        mapping_dict[identifier]= str(mapping[mapping['identifier'] == identifier]["name"].iloc[0])
    mapping_dict = json.dumps(mapping_dict)
    
    return mapping_dict
        
    
    
    
                
    
    

    
