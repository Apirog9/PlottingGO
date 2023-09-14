from django.shortcuts import render
from django.shortcuts import redirect
from django.http import HttpResponseRedirect
from django.contrib.auth.mixins import PermissionRequiredMixin
from django.views.generic.edit import CreateView
from .models import AnnotationTable
from django.db.models.signals import post_save
from .utils import gaf_transform
from .utils import annotate
from .utils import makeplot
from functools import wraps
import json
from django.urls import reverse_lazy
from django.contrib.staticfiles import finders

from enrichmentplotter.forms import InputData
from enrichmentplotter.forms import InputDataPlottingContext
from enrichmentplotter.models import ListEnrichmentInput
from enrichmentplotter.models import NameDictionary
from enrichmentplotter.models import PlotParams



# Create your views here.



def index(request):
    """View function for home page of site."""
    # Render the HTML template index.html with the data in the context variable
    return render(request, 'index.html')
    
    
    
class UploadGaf(PermissionRequiredMixin,CreateView):
    #Currently available only from admin!
    permission_required = 'enrichmentplotter.modify_data'
    model = AnnotationTable
    fields = ['anotation_name', 'raw_gaf']
    


    def model_created_or_updated(sender, created=False, **kwargs):
        print(sender)
        table = kwargs['instance']
        raw_path = table.raw_gaf.path
        name_dict = NameDictionary.objects.all()[0].name_dict
        
        new_path,gene_bp,gene_mf,gene_cc,protein_bp,protein_mf,protein_cc,updated_name_dict = gaf_transform.create_tsv(raw_path,name_dict)
        print(new_path)
        AnnotationTable.objects.filter(pk=table.pk).update(raw_gaf_path = raw_path)
        AnnotationTable.objects.filter(pk=table.pk).update(processed_gaf_path = new_path)
        AnnotationTable.objects.filter(pk=table.pk).update(gene_BP_dict = gene_bp)
        AnnotationTable.objects.filter(pk=table.pk).update(gene_MF_dict = gene_mf)
        AnnotationTable.objects.filter(pk=table.pk).update(gene_CC_dict = gene_cc)
        AnnotationTable.objects.filter(pk=table.pk).update(protein_BP_dict = protein_bp)
        AnnotationTable.objects.filter(pk=table.pk).update(protein_MF_dict = protein_mf)
        AnnotationTable.objects.filter(pk=table.pk).update(protein_CC_dict = protein_cc)
        NameDictionary.objects.all().update(name_dict = updated_name_dict)

    post_save.connect(model_created_or_updated, sender=AnnotationTable)
    success_url = reverse_lazy('index')
    
    
def input_data(request):
    #create directory from uuid
    #parse id field to list, save in directory, update precalculated_datapath_query
    #parse background field to list, save in directory, precalculated_datapath_background
    #call perform_annotation function
    #store annotation results in file, update precalculated_datapath_annotation
    #TODO update to use input from 
    
    def input_processing(sender, created=False, **kwargs):
        input = kwargs['instance']
        
        # get id type
        id_type = input.input_type
        # get correct AnnotationTable object
        table = input.GO_set
        # get correct set of term dictionaries
        dictionaries = []
        if id_type == 'Gene':
            if input.GO_BP:
                dictionaries.append(table.gene_BP_dict)
            if input.GO_MF:
                dictionaries.append(table.gene_MF_dict)
            if input.GO_CC:
                dictionaries.append(table.gene_CC_dict)
        if id_type == 'Uniprot':
            if input.GO_BP:
                dictionaries.append(table.protein_BP_dict)
            if input.GO_MF:
                dictionaries.append(table.protein_MF_dict)
            if input.GO_CC:
                dictionaries.append(table.protein_CC_dict)
        

        # get annotated id_list as json
        id_data = input.id_list
        id_data = id_data.split('\n')
        id_data = [x.strip('\r') for x in id_data]
        id_data = [x.strip(' ') for x in id_data]
        annotated_id_data = annotate.make_annotation(id_data,dictionaries)
        id_data_json = json.dumps(annotated_id_data)
        
        # get annotated background_list as json
        background_data = input.background_list
        background_data = background_data.split('\n')
        background_data = [x.strip('\r') for x in background_data]
        background_data = [x.strip(' ') for x in background_data]
        annotated_background_data = annotate.make_annotation(background_data,dictionaries)
        background_data_json = json.dumps(annotated_background_data)
        
        # update fields
        ListEnrichmentInput.objects.filter(pk=input.pk).update(query_items_annotation = id_data_json)
        ListEnrichmentInput.objects.filter(pk=input.pk).update(background_items_annotation = background_data_json)
        
    
       
    post_save.connect(input_processing, sender=ListEnrichmentInput)

    if request.method == 'POST':
        form = InputData(request.POST)
        
        if form.is_valid():
            data = form.save()
            return redirect('data_ready_to_plot',data.query_id)
    else:
         form = InputData()

    return render(request, 'enrichmentplotter/input_data.html',  {'form': form})
    
#load file with annotation
# get matched number
# get unmatched number
#
def data_ready(request,input_instance_id):
    input_instance = ListEnrichmentInput.objects.filter(pk=input_instance_id).values().get()
    
    query_annotation = json.loads(input_instance['query_items_annotation'])
    background_annotation = json.loads(input_instance['background_items_annotation'])
    query_length = len(list(query_annotation.keys()))
    background_length = len(list(background_annotation.keys()))
    not_annotated_query = [x for x in query_annotation.keys() if query_annotation[x] == []]
    not_annotated_background = [x for x in background_annotation.keys() if background_annotation[x] == []]
    not_annotated_query_num = len(not_annotated_query)
    not_annotated_background_num = len(not_annotated_query)
    context = {
    'query_length':query_length,
    'background_length':background_length,
    'not_annotated_query_num':not_annotated_query_num,
    'not_annotated_background_num':not_annotated_background_num,
    'input_instance_id':input_instance_id
    }



    return render(request, 'enrichmentplotter/data_ready.html',context=context)
    
def data_plotted(request,input_instance_id):

    input_instance = ListEnrichmentInput.objects.filter(pk=input_instance_id).values().get()
    query_annotation = json.loads(input_instance['query_items_annotation'])
    background_annotation = json.loads(input_instance['background_items_annotation'])
    query_length = len(list(query_annotation.keys()))
    background_length = len(list(background_annotation.keys()))
    not_annotated_query = [x for x in query_annotation.keys() if query_annotation[x] == []]
    not_annotated_background = [x for x in background_annotation.keys() if background_annotation[x] == []]
    not_annotated_query_num = len(not_annotated_query)
    not_annotated_background_num = len(not_annotated_query)
    context = {
    'query_length':query_length,
    'background_length':background_length,
    'not_annotated_query_num':not_annotated_query_num,
    'not_annotated_background_num':not_annotated_background_num,
    'input_instance_id':input_instance_id
    }
    # if the page is opened, create params model
    if request.method == 'GET' and not PlotParams.objects.filter(input_data=input_instance_id).exists():
        print('getmethod')
        plot_params = PlotParams(input_data=ListEnrichmentInput.objects.filter(pk=input_instance_id).get())
        maxp = plot_params.max_p_value
        minnum = plot_params.min_protein_number
        minfrac = plot_params.min_fraction
        frame,frameshape = makeplot.prepare_barplot(query_annotation,background_annotation,maxp,minnum,minfrac)
        nrows = frameshape[0]
        plot_params.font_size_x_title = 8 + (8/nrows) 
        plot_params.font_size_y_title = 8 + (8/nrows)
        plot_params.font_size_x_labels = 5 + (6/nrows)
        plot_params.font_size_y_labels = 5 + (6/nrows)
        plot_params.font_size_plot_title = 10 + (8/nrows)
        plot_params.save()
        form = InputDataPlottingContext(instance=plot_params)
    if request.method == 'POST':
        print('post')
        plot_params = PlotParams.objects.filter(input_data=input_instance_id).get()
        form = InputDataPlottingContext(request.POST,instance=plot_params)
        if form.is_valid():
            form.save()
            return redirect('data_plotted',input_instance_id)
            
        


    # get plotting parameters back from plotparams model
    plot_params = PlotParams.objects.filter(input_data=input_instance_id).get()
    maxp = plot_params.max_p_value
    minnum = plot_params.min_protein_number
    minfrac = plot_params.min_fraction
    frame,frameshape = makeplot.prepare_barplot(query_annotation,background_annotation,maxp,minnum,minfrac)
    nrows = frameshape[0]
    x_t_s =  plot_params.font_size_x_title
    y_t_s =  plot_params.font_size_y_title
    x_l_s =  plot_params.font_size_x_labels
    y_l_s =  plot_params.font_size_y_labels
    p_t_s =  plot_params.font_size_plot_title
    f_x_s =  plot_params.figure_size_x
    f_y_s =  plot_params.figure_size_y
    cmap =   plot_params.colormap
    if nrows>0:
        plotimg = makeplot.make_barplot(frame,x_t_s,y_t_s,x_l_s,y_l_s,p_t_s,f_x_s,f_y_s,cmap)
    else:
        plotimg = None
    
    form = InputDataPlottingContext(instance=plot_params)
    
   
    
    
    
    context['plot'] = plotimg
    context['form'] = form
    
    
    return render(request, 'enrichmentplotter/data_plotted.html',context=context)
    
    
class UploadNameDictionary(PermissionRequiredMixin,CreateView):

    permission_required = 'enrichmentplotter.modify_data'
    model = NameDictionary
    fields = ['base_table_name', 'base_table']
    


    def model_created_or_updated(sender, created=False, **kwargs):
        table = kwargs['instance']
        raw_path = table.base_table.path
        
        name_dictionary = gaf_transform.get_namedict(raw_path)
        NameDictionary.objects.filter(pk=table.pk).update(base_table_path = raw_path)
        NameDictionary.objects.filter(pk=table.pk).update(name_dict = name_dictionary)


    post_save.connect(model_created_or_updated, sender=NameDictionary)
    success_url = reverse_lazy('index')






    
    