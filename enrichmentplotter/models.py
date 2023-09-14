from django.db import models
import uuid

from django.db.models.signals import post_save

# Create your models here.


class NameDictionary(models.Model):
    # this is used only to refer to and update
    # term:name mapping as new terms are added
    # will be possible to update names as necessary
    # table should be uploaded from admin , and only one should exist
    base_table_name = models.CharField(max_length=200, primary_key=True)
    base_table = models.FileField(upload_to='GO_tables', max_length=254)
    base_table_path = models.CharField(max_length=1000,null=True,blank=True,default='')
    name_dict = models.JSONField(null=True,blank=True)
    


class AnnotationTable(models.Model):
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
        permissions = (("modify_data", "Modify GO dataset"),)
        


class ListEnrichmentInput(models.Model):
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
    



#class PlotParameters(models.Model):
