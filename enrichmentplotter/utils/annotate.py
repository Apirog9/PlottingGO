# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 21:17:46 2023

@author: APirog
"""
import json


def make_annotation(ids,dictionary_list):
    annotation = {identifier:[] for identifier in ids}
    for dictionary in dictionary_list:
        dictionary = json.loads(dictionary)
        for identifier in ids:
            try:
                annotations = dictionary[identifier]
            except KeyError:
                annotations = []
            annotation[identifier] = annotation[identifier] + annotations
    return annotation
    
    
    