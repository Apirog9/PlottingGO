# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 21:17:46 2023

@author: APirog
"""
import json


def make_annotation(ids,dictionary_list):
    '''
    ids - list of identifiers
    dictionary_list - list of JSON strings containing GO term annotation dictionaries to use
    Returns:
    annotation - Python dictionary containing identifier:list of annotated terms pairs
    '''
    annotation = {identifier:[] for identifier in ids}
    for dictionary in dictionary_list:
        dictionary = json.loads(dictionary)
        for identifier in ids:
            try:
                annotations = dictionary[identifier]
            # if identifier absent from annotation, return empty list
            except KeyError:
                annotations = []
            annotation[identifier] = annotation[identifier] + annotations
    return annotation
    
    
    