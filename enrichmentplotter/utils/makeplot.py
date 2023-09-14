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





#def annotate(data,column,lista):
#    lista = open(lista,'r').read().splitlines()
#    data['selected'] = data[column].isin(lista)
#    return data



def calculate_fisher_exact(data,categorycolumn,annotationcolumn,pvalcutoff):

    selected = data[data[annotationcolumn]==True]
    not_selected = data[data[annotationcolumn]==False]
    no_selitems = list(not_selected[categorycolumn])
    selitems = list(selected[categorycolumn])
    selitems = [list(set(str(x).split(';'))) for x in selitems]
    no_selitems = [list(set(str(x).split(';'))) for x in no_selitems]
    selitems_names = list(it.chain.from_iterable(selitems))
    for_calc = list(set(selitems_names))
    for_calc = [x for x in for_calc if x!= '']

    
    datadict = {}
    for term in for_calc:
        table = np.zeros((2,2),dtype=np.int32())
        in_bkgd = len([x for x in no_selitems if term in x  ])
        in_sample = len([x for x in selitems if term in x  ])
        rest_bkgd = len([x for x in no_selitems if term not in x  ])
        rest_sample = len([x for x in selitems if term not in x  ])

        table[0,0] = in_bkgd
        table[0,1] = rest_bkgd
        table[1,0] = in_sample
        table[1,1] = rest_sample
        suma = in_bkgd + in_sample
        res = fisher_exact(table)
        pval = res[1]
        assoc = association(table)
        pseudoassoc = (table[1][0] / (table[1][0]+table[0][0]))*100
        datadict[term] = [table,pval,assoc,suma,pseudoassoc]

    
    return datadict


def nicebarplot_guess_parameters(frame,GO_name_column,p_val_column,association_column,protein_number_column,maxp,minnum,minfrac):
    frame = frame[frame[p_val_column] <= maxp]
    frame = frame[frame[protein_number_column] >=minnum]
    frame = frame[frame[association_column] >=minfrac]
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
    pl.theme(legend_position = (0.14,0.87),
             legend_direction='vertical',
             legend_title=pl.element_text(size=x_t_s,color='black',weight='bold'),
             legend_text=pl.element_text(size=p_t_s),
             legend_title_align='center',
             legend_key_width=2*x_t_s,
             legend_key_height=2*x_t_s,
             figure_size=(f_x_s,f_y_s))+\
    pl.theme(panel_background = pl.element_rect(fill = 'white', colour = 'white'))+\
    pl.scales.scale_fill_continuous(cmap_name=cmap)+\
    pl.geom_text(
     pl.aes(x=GO_name_column,y=association_column,label=protein_number_column),
     stat='identity',nudge_y=x_l_s/3,size=x_l_s,ha='left')
    
    
    
    #plot.save('barplot.png',width=11,height=16,units='cm',dpi=1000)
    
    return plot
    
    
    



def prepare_barplot(query_annotation,background_annotation,maxp,minnum,minfrac):
    plotting_frame = pd.DataFrame(columns = ['name','p_value','association','protein_number','GO','fraction_in_valid'])
    #prepare fisher_exact input
    frame = pd.DataFrame(columns = ['name','terms','is_query'])
    framedict = []
    for annotation in query_annotation.keys():
        line = {}
        line['name'] = annotation
        line['terms'] = ';'.join(query_annotation[annotation])
        line['is_query'] = True
        framedict.append(line)
        
    for annotation in background_annotation.keys():
        line = {}
        line['name'] = annotation
        line['terms'] = ';'.join(background_annotation[annotation])
        line['is_query'] = False
        framedict.append(line)
        
    result_frame = pd.DataFrame(framedict)
    frame = pd.concat([frame,result_frame])
    frame.to_csv('testframe.tsv',sep='\t',index=False)
    datadict = calculate_fisher_exact(frame,'terms', 'is_query', 0.01)

    
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
    corrected_bh = multipletests(plotting_frame['p_value'],method = 'fdr_bh')
    plotting_frame['corrected p value'] = corrected_bh[1]
    frame,frameshape = nicebarplot_guess_parameters(plotting_frame,'name','corrected p value','fraction_in_valid','protein_number',maxp,minnum,minfrac)
    print(frame)
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
    
    plot = nicebarplot(plotting_frame,'name','corrected p value','fraction_in_valid','protein_number',
                       x_t_s,
                       y_t_s,
                       x_l_s,
                       y_l_s,
                       p_t_s,
                       f_x_s,
                       f_y_s,
                       cmap)
    buf = BytesIO()
    plot.save(buf, format='png',dpi=200)
    buf.seek(0)
    string = base64.b64encode(buf.read())
    return string.decode('utf-8')
    
    


        

    
    
    
    
    
    
