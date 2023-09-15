from django.db import models
import uuid

from django.db.models.signals import post_save

# Create your models here.


class NameDictionary(models.Model):
    '''
    Model contains data about mapping numeric GO identifiers to readable names
    Use base table to fill the model with set of common names, stores path to table file
    Actual dictionary is stored in JSON and updated within /utilities/gaf_transform.py
    Only one model should exist and should be created in admin panel
    '''
    base_table_name = models.CharField(max_length=200, primary_key=True)
    base_table = models.FileField(upload_to='GO_tables', max_length=254)
    base_table_path = models.CharField(max_length=1000,null=True,blank=True,default='')
    name_dict = models.JSONField(null=True,blank=True)
    


class AnnotationTable(models.Model):
    '''
    Model contains data about actual GO annotations. Currently works well for human tables and other ones which use gene or uniprot identifiers as 
    an annotation bases
    annotation_name - should contain organism and term set name (e.g. complete, slim,fat)
    raw_gaf - .gaf file to create table and annotations
    raw_gaf_path - store path to file to be transformed with /utilities/gaf_transform.py
    processed_gaf_path -  path to processed tsv-like representation of GAF file
    *_dict - jsons containing actual mappings to use by /utilities/annotate.py
    
    '''
    anotation_name = models.CharField(max_length=200, primary_key=True)
    raw_gaf = models.FileField(upload_to='GO_tables', max_length=254)
    raw_gaf_path = models.CharField(max_length=1000,null=True,blank=True,default='')
    processed_gaf_path = models.CharField(max_length=1000,null=True,blank=True,default='')
    gene_BP_dict = models.JSONField(null=True,blank=True)
    gene_MF_dict = models.JSONField(null=True,blank=True)
    gene_CC_dict = models.JSONField(null=True,blank=True)
    protein_BP_dict = models.JSONField(null=True,blank=True)
    protein_MF_dict = models.JSONField(null=True,blank=True)
    protein_CC_dict = models.JSONField(null=True,blank=True)
    
    
    
    class Meta:
    # may be useful to allow logged user to add table, without admin permission
        permissions = (("modify_data", "Modify GO dataset"),)
        


class ListEnrichmentInput(models.Model):
    '''
    Model contains actual user input to perform enrichment analysis
    CHOICES - identifier choices, currently work for human and other GO data that use Uniprot ID or gene as an annotation base
    query_id - unique query id to trace whole transformation and give a semi-permanent link to results
    id_list - unlimited field containing one accepted ID per line
    background_list - unlimited field containing one accepted ID per line
    GO_set - gene ontology set, should match organism of the list&background_list
    GO_BP,GO_MF,GO_CC -  fields to check whether user decided to include those classes for plottins
    input type - select Uniprt Id or Gene
    query_items_annotation,background_items_annotation - annotations generated by /utilities/annotate.py
    precalculated_datapath_annotation - now unused, may be used to store media files created by other functions, if not compatible with database
    '''
    CHOICES = (('Gene', 'Gene'),('Uniprot', 'Uniprot ID'))


    query_id = models.UUIDField(primary_key=True, default=uuid.uuid4)
    id_list = models.TextField()
    background_list = models.TextField()
    GO_set = models.ForeignKey(AnnotationTable, on_delete=models.CASCADE,null=True)
    GO_BP = models.BooleanField()
    GO_MF = models.BooleanField()
    GO_CC = models.BooleanField()
    input_type = models.CharField(max_length=100,choices=CHOICES,default='Gene')
    query_items_annotation = models.JSONField(null=True)
    background_items_annotation = models.JSONField(null=True)
    precalculated_datapath_annotation = models.CharField(max_length=2000,blank=True)
    
    
class PlotParams(models.Model):
    '''
    Model contains default, precalculated or user-provided plotting parameters
    cmaps - matplotlib colormap choices
    params_id -  model identifier
    input_data - ForeginKey reference to single user input
    max_p_value - maximum CORRECTED p-value for a term to be plotted
    min_protein_number - minimum query+background protein/gene number for a term to be plotted
    min_fraction -  minimum fraction of query proteins in all proteins for a term to be plotted
    colormap - matplotlib colormap to reflect CORRECTED p value
    font_size_x_labels,font_size_y_labels,font_size_x_title,font_size_y_title,font_size_plot_title -font sizes for plot
    figure_size_x,figure_size_y - figure size in inches
    '''
    cmaps = (('viridis', 'viridis'),('plasma', 'plasma'))
    params_id = models.UUIDField(primary_key=True, default=uuid.uuid4)
    input_data =  models.ForeignKey('ListEnrichmentInput', on_delete=models.CASCADE, null=True)
    max_p_value = models.FloatField(default = 0.05)
    min_protein_number = models.IntegerField(default = 3)
    min_fraction = models.IntegerField(default = 10)
    colormap = models.CharField(max_length=100,choices=cmaps,default='viridis')
    font_size_x_labels = models.FloatField(null=True,blank=True)
    font_size_y_labels = models.FloatField(null=True,blank=True)
    font_size_x_title = models.FloatField(null=True,blank=True)
    font_size_y_title = models.FloatField(null=True,blank=True)
    font_size_plot_title = models.FloatField(null=True,blank=True)
    figure_size_x = models.FloatField(default=4)
    figure_size_y = models.FloatField(default=8)
    

