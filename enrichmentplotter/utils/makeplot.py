# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 16:04:57 2023

@author: APirog
"""
import pandas as pd
import itertools as it
from scipy.stats import fisher_exact
from scipy.stats.contingency import association
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import plotnine as pl
import numpy as np
import json
import base64
from io import BytesIO



def calculate_fisher_exact(data,categorycolumn,annotationcolumn):
    '''
    calculate Fisher exact test results for contingency tables generated for
    query vs background term frequencies
    
    data - pandas dataframe, row names are identifiers/genes, but they are not used here yet
    categorycolumn - column with term names, separated by ';'
    annotationcolunm -  boolean column, True if row belongs to query, False if not
    
    Returns
    
    dictionary
    term : [contingency table,p-value,association,number of items,pseudoassociation(percent of query proteins in total)]
    '''
    # define query/background rowsets
    selected = data[data[annotationcolumn]==True]
    not_selected = data[data[annotationcolumn]==False]
    
    # split into list of lists
    no_selitems = list(not_selected[categorycolumn])
    selitems = list(selected[categorycolumn])
    selitems = [list(set(str(x).split(';'))) for x in selitems]
    no_selitems = [list(set(str(x).split(';'))) for x in no_selitems]
    
    # generate list of term names present in query
    selitems_names = list(it.chain.from_iterable(selitems))
    for_calc = list(set(selitems_names))
    for_calc = [x for x in for_calc if x!= '']

    # generate result dictionary
    datadict = {}
    for term in for_calc:
        #initiate and fill contingency table with counts
        table = np.zeros((2,2),dtype=np.int32())
        in_bkgd = len([x for x in no_selitems if term in x  ])
        in_sample = len([x for x in selitems if term in x  ])
        rest_bkgd = len([x for x in no_selitems if term not in x  ])
        rest_sample = len([x for x in selitems if term not in x  ])
        
        table[0,0] = in_bkgd
        table[0,1] = rest_bkgd
        table[1,0] = in_sample
        table[1,1] = rest_sample
        # calculate sum of identifiers in query and background, where term was present
        suma = in_bkgd + in_sample
        # calculate fisher exact test
        res = fisher_exact(table)
        # get p-value
        pval = res[1]
        # get association TODO check method, probably Cramer
        assoc = association(table)
        # calculate pseudoassociation (percentage of identifiers in query)
        pseudoassoc = (table[1][0] / (table[1][0]+table[0][0]))*100
        # add entry to dictionary
        datadict[term] = [table,pval,assoc,suma,pseudoassoc]

    return datadict


def nicebarplot_guess_parameters(frame,GO_name_column,p_val_column,association_column,protein_number_column,maxp,minnum,minfrac):
    '''
    Perform initial filtering to enable guesswork to adjust plot font size and plot size(TODO plot size changes)
    
    Returns
    dataframe for plotting, shape of dataframe
    '''
    # filter rows by pvalue, minimum protein number, minmum protein fraction in query
    frame = frame[frame[p_val_column] <= maxp]
    frame = frame[frame[protein_number_column] >=minnum]
    frame = frame[frame[association_column] >=minfrac]
    # sort values for proper barplot creation TODO reverse sorting of bars may be easily implemented
    frame = frame.sort_values(by=association_column)
    return frame,frame.shape
    
    
    
def nicebarplot(frame,GO_name_column,p_val_column,association_column,protein_number_column,
                x_t_s,
                y_t_s,
                x_l_s,
                y_l_s,
                p_t_s,
                f_x_s,
                f_y_s,
                cmap):
    '''
    Plot generation in plotnine
    
    frame - dataframe prepared by prepare barplot
    GO_name_column - column with GO term names
    p_val_column - column with corrected p value
    association_column - column with association measure (or pseudoasociation, like 'fraction_in_valid')
    protein_number_column - column with number of protein in query+background
    x_t_s - font size x axis title
    y_t_s - font size y axis title
    x_l_s - font size x axis labels
    y_l_s - font size y axis labels
    p_t_s - font size plot title
    f_x_s - figure size x
    f_y_s - figure size y
    cmap  - colormap (for p-value)
    
    
    Returns
    plotnine figure. (may be easily transformed to matplotlib figure or png)
    
    '''

    xlabels = frame[GO_name_column]
    plot = pl.ggplot(frame) + pl.geom_bar(pl.aes(x=GO_name_column,y=association_column,fill=p_val_column),stat = 'identity')+\
    pl.scale_x_discrete(limits=xlabels)+\
    pl.coord_flip()+\
    pl.labels.xlab('Gene ontology term')+\
    pl.labels.ylab('Percent of significantly\n regulated proteins')+\
    pl.theme(axis_text_x=pl.element_text(size=x_l_s,color='black'))+\
    pl.theme(axis_text_y=pl.element_text(size=y_l_s,color='black'))+\
    pl.theme(axis_title_y=pl.element_text(size=x_t_s,color='black'))+\
    pl.theme(axis_title_x=pl.element_text(size=y_t_s,color='black'))+\
    pl.theme(panel_grid_minor_y=pl.element_line(color="lightgray"))+\
    pl.theme(panel_grid_minor_x=pl.element_line(color="lightgray"))+\
    pl.theme(panel_grid_major_y=pl.element_line(color="lightgray",size=0.5))+\
    pl.theme(panel_grid_major_x=pl.element_line(color="lightgray",size=0.5))+\
    pl.theme(axis_ticks_length=0)+\
    pl.theme(legend_position = (0.87,0.25),
             legend_direction='vertical',
             legend_title=pl.element_text(size=x_t_s,color='black',weight='bold'),
             legend_text=pl.element_text(size=x_t_s),
             legend_title_align='left',
             legend_key_width=0.5*x_t_s,
             legend_key_height=0.5*x_t_s,
             figure_size=(f_x_s,f_y_s))+\
    pl.theme(panel_background = pl.element_rect(fill = 'white', colour = 'white'))+\
    pl.scales.scale_fill_continuous(cmap_name=cmap)+\
    pl.geom_text(pl.aes(x=GO_name_column,y=association_column,label=protein_number_column),
                 stat='identity',
                 nudge_y=x_l_s/9,
                 size=x_l_s-1,
                 ha='left')
    
    
    #plot.save('barplot.png',width=11,height=16,units='cm',dpi=1000)
    
    return plot
    
    
    



def prepare_barplot(query_annotation,background_annotation,maxp,minnum,minfrac):
    '''
    Prepare data for plotting from annotation
    query_annotation - dictionary containing identifier:GO term names list pairs
    background_annotation - dictionary containing identifier:GO term names list pairs
    maxp - maximum CORRECTED p-value for a term to be plotted
    minnum - minimum query+background protein/gene number for a term to be plotted
    minfrac -  minimum fraction of query proteins in all proteins for a term to be plotted
    
    Returns 
    dataframe for nicebarplot
    shape of this dataframe
    '''
    # initialize empty dataframe for plotting - defined for clarity, not needed!
    plotting_frame = pd.DataFrame(columns = ['name','p_value','association','protein_number','GO','fraction_in_valid'])
    #prepare fisher_exact input
    # initialize empty dataframe for fisher_exact - defined for clarity, not needed!
    frame = pd.DataFrame(columns = ['name','terms','is_query'])
    framedict = []
    # transform dictionary to dataframe part for query
    for ident in query_annotation.keys():
        line = {}
        line['name'] = ident
        line['terms'] = ';'.join(query_annotation[ident])
        line['is_query'] = True
        framedict.append(line)
    # transform dictionary to dataframe part for background 
    for ident in background_annotation.keys():
        # prevent adding query items to background if identical item exist in both, would decrease test sensitivity
        if ident not in list(query_annotation.keys()):
            line = {}
            line['name'] = ident
            line['terms'] = ';'.join(background_annotation[ident])
            line['is_query'] = False
            framedict.append(line)
    #  make dataframe for fisher_exact from lsit of dictionaries    
    result_frame = pd.DataFrame(framedict)
    frame = pd.concat([frame,result_frame])
    frame.to_csv('testframe.tsv',sep='\t',index=False)
    #calculate fisher_exact output
    datadict = calculate_fisher_exact(frame,'terms', 'is_query')

    # transform fisher_exact output to dataframe
    framedict_input = []
    for key in datadict.keys():
        line = {}
        line['name'] = key
        line['p_value'] = datadict[key][1]
        line['association'] = datadict[key][2]
        line['protein_number'] = datadict[key][3]
        line['fraction_in_valid'] = datadict[key][4]
        #TODO add process class later for annotation
        #line['GO'] = process_name
        framedict_input.append(line)
    result_frame_input = pd.DataFrame(framedict_input)
    plotting_frame = pd.concat([plotting_frame,result_frame_input])
    # make correction for multiple testing by Benjamini-Hochberg method
    corrected_bh = multipletests(plotting_frame['p_value'],method = 'fdr_bh')
    plotting_frame['corrected p value'] = corrected_bh[1]
    # perform initial filtering to enable guessing some plot sizing parameters
    frame,frameshape = nicebarplot_guess_parameters(plotting_frame,'name','corrected p value','fraction_in_valid','protein_number',maxp,minnum,minfrac)
    
    return frame,frameshape
        
        
def make_barplot(plotting_frame,
                 x_t_s,
                 y_t_s,
                 x_l_s,
                 y_l_s,
                 p_t_s,
                 f_x_s,
                 f_y_s,
                 cmap):
    '''
    Perform actual plotting and transform to form compatible with embedding in template
    '''
    # generate plot as plotnine figure
    plot = nicebarplot(plotting_frame,'name','corrected p value','fraction_in_valid','protein_number',
                       x_t_s,
                       y_t_s,
                       x_l_s,
                       y_l_s,
                       p_t_s,
                       f_x_s,
                       f_y_s,
                       cmap)
                       
    # save plot in png format to buffer                  
    buf = BytesIO()
    plot.save(buf, format='png',dpi=200)
    buf.seek(0)
    # puts the contents into correct format
    string = base64.b64encode(buf.read())
    string = string.decode('utf-8')
    return string
    
    


        

    
    
    
    
    
    
