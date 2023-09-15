# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 16:53:52 2023

@author: APirog
"""
import json
import pandas as pd
import requests, sys


def get_by_API(term_id):
    '''
    Accessory function
    Attempt to retrieve a name for a GO term using QuickGO REST API 
    '''
    # create URL and get response
    requestURL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" +term_id
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    # if response is ok, get name from json
    if r.ok:
        data = json.loads(r.text)
        name = data["results"][0]["name"]
    # if something went wrong, return name as None
    else:
        r.raise_for_status()
        sys.exit()
        name = None
        
    return name
    
    
def create_tsv(gafpath,mapping_dict):
    '''
    Trensforms raw gaf file to identifier:terms pairs. such data are easily accessible and
    can be stored as JSON strings in sqlite database
    
    
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
    json object with updated GO identifier : name dictionary

    '''
    
    def getmap(key,dictionary):
        '''
        Accessory function
        simple retrieve dictionary value or None if key not exist
        '''
        try:
            value = dictionary[key]
        except KeyError:
                value = None  
        return value
    
    def check_dict(go_nums,mapping_dict):
        '''
        Accessory function
        Check completness of current GO identifier : GO term name dictionary
        try to update dictionary using QuickGO REST API 
        Returns
        updated GO identifier : GO name dictionary
        '''
        # check for non-existing identifiers
        current = list(mapping_dict.keys())
        novel_ids = [x for x in go_nums if not x in current]
        for go_term in novel_ids:
        # update if name retireval succeed
            newname = get_by_API(go_term)
            if newname:
                mapping_dict[go_term] = get_by_API(go_term)
                
        return mapping_dict
        
        
    
    #write tsv-like form of gaf file and load it back as dataframe
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
    
    # remove non-protein hits, they need special treatment, not yet implemented
    ontology_dataframe = ontology_dataframe[ontology_dataframe['molecule_type'] == 'protein']
    # fill NAs with empty strings in GO_function (reference like 'colocalizes_to') column
    ontology_dataframe['GO_function'] = ontology_dataframe['GO_function'].fillna('')
    # remove negative mappings, they need special treatment, not yet implemented
    ontology_dataframe = ontology_dataframe[~ontology_dataframe['GO_function'].str.contains('NOT')]

    #translate to names
    # get unique numerical ontology identifiers
    ontology_ids = list(ontology_dataframe['GO_term'].unique())
    # load dictionary from JSON string
    mapping_dict = json.loads(mapping_dict)
    # check name dictionary completness and if possible update by QuickGO REST API
    mapping_dict = check_dict(ontology_ids,mapping_dict)
    # get names
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
        print(gene)
        dataslice = frame_BP[frame_BP['Gene']==gene]
        terms = list(dataslice['GO_term_name'])
        gene_BP_dict[gene] = list(set(terms))
    gene_BP_dict = json.dumps(gene_BP_dict)
    
    gene_MF_dict = {}
    for gene in genes:
        print(gene)
        dataslice = frame_MF[frame_MF['Gene']==gene]
        terms = list(dataslice['GO_term_name'])
        gene_MF_dict[gene] = list(set(terms))
    gene_MF_dict = json.dumps(gene_MF_dict)
    
    gene_CC_dict = {}
    for gene in genes:
        print(gene)
        dataslice = frame_CC[frame_CC['Gene']==gene]
        terms = list(dataslice['GO_term_name'])
        gene_CC_dict[gene] = list(set(terms))
    gene_CC_dict = json.dumps(gene_CC_dict)
    
    protein_BP_dict = {}
    for protein in proteins:
        print(protein)
        dataslice = frame_BP[frame_BP['Identifier']==protein]
        terms = list(dataslice['GO_term_name'])
        protein_BP_dict[protein] = list(set(terms))
    protein_BP_dict = json.dumps(protein_BP_dict)
    
    protein_MF_dict = {}
    for protein in proteins:
        print(protein)
        dataslice = frame_MF[frame_MF['Identifier']==protein]
        terms = list(dataslice['GO_term_name'])
        protein_MF_dict[protein] = list(set(terms))
    protein_MF_dict = json.dumps(protein_MF_dict)
    
    protein_CC_dict = {}
    for protein in proteins:
        print(protein)
        dataslice = frame_CC[frame_CC['Identifier']==protein]
        terms = list(dataslice['GO_term_name'])
        protein_CC_dict[protein] = list(set(terms))
    protein_CC_dict = json.dumps(protein_CC_dict)
    
    mapping_dict = json.dumps(mapping_dict)
    return newfile,gene_BP_dict,gene_MF_dict,gene_CC_dict,protein_BP_dict,protein_MF_dict,protein_CC_dict,mapping_dict



def get_namedict(name_dict_path):
    '''
    transform csv file to json object containing GO identifier:GO term name dictionary
    used to populate initial dictionary
    name_dict_path - path to csv file
    Returns:
    json object with initial GO identifier : name dictionary
    '''
    mapping = pd.read_csv(name_dict_path, sep = ' ')
    GO_identifiers = mapping["identifier"].unique()
    mapping_dict = {}
    for identifier in GO_identifiers:
        mapping_dict[identifier]= str(mapping[mapping['identifier'] == identifier]["name"].iloc[0])
    mapping_dict = json.dumps(mapping_dict)
    
    return mapping_dict
        
    
    
    
                
    
    

    
