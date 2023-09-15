from django.shortcuts import render
from django.shortcuts import redirect
from django.http import HttpResponseRedirect
from django.contrib.auth.mixins import PermissionRequiredMixin
from django.views.generic.edit import CreateView
from .models import AnnotationTable
from django.db.models.signals import post_save
# GAF file upload and annotation data model (AnnotationTable) input creation
from .utils import gaf_transform
# user input annotation with AnnotationTable and creating input for ListEnrichmentInput model creation
from .utils import annotate
# plotting data preparation and actual plotting
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
    '''
    Simple start page view
    '''
    # Render the HTML template index.html with the data in the context variable
    return render(request, 'index.html')
    
    
class UploadGaf(PermissionRequiredMixin,CreateView):
    '''
    Defines both page (experimental) and admin site behaviour for uploading gaf files
    while uploading only name and file should be provided by user or admin
    '''
    permission_required = 'enrichmentplotter.modify_data'
    model = AnnotationTable
    fields = ['anotation_name', 'raw_gaf']
    


    def model_created_or_updated(sender, created=False, **kwargs):
        '''
        defines transformations performed after model creation
        '''
        # get path to raw gaf file
        table = kwargs['instance']
        raw_path = table.raw_gaf.path
        # get NameDictionary model (only one should exist)
        name_dict = NameDictionary.objects.all()[0].name_dict
        
        # transform gaf file to .tsv and dictionaries used for annotation
        new_path,gene_bp,gene_mf,gene_cc,protein_bp,protein_mf,protein_cc,updated_name_dict = gaf_transform.create_tsv(raw_path,name_dict)
        # update fields 
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
    # after creation, go back to welcome page, TODO not implemented completely
    success_url = reverse_lazy('index')
    
    
def input_data(request):
    '''
    define processing of user input and user input form creation
    '''
    def input_processing(sender, created=False, **kwargs):
        input = kwargs['instance']
        # get id type
        id_type = input.input_type
        # get correct AnnotationTable object
        table = input.GO_set
        # get correct set of term dictionaries depending on id_type
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
        
        # get annotated id_list as json by annotate.make_annotation
        id_data = input.id_list
        id_data = id_data.split('\n')
        id_data = [x.strip('\r') for x in id_data]
        id_data = [x.strip(' ') for x in id_data]
        annotated_id_data = annotate.make_annotation(id_data,dictionaries)
        id_data_json = json.dumps(annotated_id_data)
        
        # get annotated background_list as json by annotate.make_annotation
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
    
    # form creation and management
    if request.method == 'POST':
        form = InputData(request.POST)
        if form.is_valid():
            # if form is valid, perform data annotation and proceed to data ready page - use query_id to keep the data id connected 
            # with other parameters!
            data = form.save()
            return redirect('data_ready_to_plot',data.query_id)
    else:
         form = InputData()

    return render(request, 'enrichmentplotter/input_data.html',  {'form': form})
    

def data_ready(request,input_instance_id):
    '''
    defines redy to plot view page, showing some annotation data and propmting to plot
    '''
    # get input_instance id
    input_instance = ListEnrichmentInput.objects.filter(pk=input_instance_id).values().get()
    # get annotations
    query_annotation = json.loads(input_instance['query_items_annotation'])
    background_annotation = json.loads(input_instance['background_items_annotation'])
    # get some numbers, query and backgroud item numbers, and number of not annotated items (no item in gaf file, or no terms in selected classes)
    # for this item
    query_length = len(list(query_annotation.keys()))
    background_length = len(list(background_annotation.keys()))
    not_annotated_query = [x for x in query_annotation.keys() if query_annotation[x] == []]
    not_annotated_background = [x for x in background_annotation.keys() if background_annotation[x] == []]
    not_annotated_query_num = len(not_annotated_query)
    not_annotated_background_num = len(not_annotated_background)
    # prepare page rendering context inpput_instance_id traces input
    context = {
    'query_length':query_length,
    'background_length':background_length,
    'not_annotated_query_num':not_annotated_query_num,
    'not_annotated_background_num':not_annotated_background_num,
    'input_instance_id':input_instance_id
    }


    return render(request, 'enrichmentplotter/data_ready.html',context=context)
    
def data_plotted(request,input_instance_id):
    '''
    defines plotted data view, together with form to update and re-plot data with updated parameters
    including form is experimental at this point
    '''
    # get data identical as for ready_to_plot , may or may not be used in template
    input_instance = ListEnrichmentInput.objects.filter(pk=input_instance_id).values().get()
    query_annotation = json.loads(input_instance['query_items_annotation'])
    background_annotation = json.loads(input_instance['background_items_annotation'])
    query_length = len(list(query_annotation.keys()))
    background_length = len(list(background_annotation.keys()))
    not_annotated_query = [x for x in query_annotation.keys() if query_annotation[x] == []]
    not_annotated_background = [x for x in background_annotation.keys() if background_annotation[x] == []]
    not_annotated_query_num = len(not_annotated_query)
    not_annotated_background_num = len(not_annotated_background)
    # generate context, at this point identical to ready_to_plot
    context = {
    'query_length':query_length,
    'background_length':background_length,
    'not_annotated_query_num':not_annotated_query_num,
    'not_annotated_background_num':not_annotated_background_num,
    'input_instance_id':input_instance_id
    }
    # if the page is opened, and PlotParams model connected to input data not yest exist create params model
    # there could be only one PlotParams model per input data
    if request.method == 'GET' and not PlotParams.objects.filter(input_data=input_instance_id).exists():
        # create model
        plot_params = PlotParams(input_data=ListEnrichmentInput.objects.filter(pk=input_instance_id).get())
        # get defaults for plotting 
        maxp = plot_params.max_p_value
        minnum = plot_params.min_protein_number
        minfrac = plot_params.min_fraction
        # create plotting dataframe and frameshape to guess some plotting parameters (font sie, plot size)
        frame,frameshape = makeplot.prepare_barplot(query_annotation,background_annotation,maxp,minnum,minfrac)
        nrows = frameshape[0]
        #guess font size
        plot_params.font_size_x_title = int(4 + (8/nrows)) 
        plot_params.font_size_y_title = int(4 + (8/nrows))
        plot_params.font_size_x_labels = int(4 + (6/nrows))
        plot_params.font_size_y_labels = int(4 + (6/nrows))
        plot_params.font_size_plot_title = int(4 + (8/nrows))
        plot_params.figure_size_x = 6
        plot_params.figure_size_y = 1 + 0.125*nrows
        # save model
        plot_params.save()
        # create InputDataPlottingContext form with defaults or initial precalculations
        form = InputDataPlottingContext(instance=plot_params)
    # if user posted form !TODO check if this explanation is right!
    if request.method == 'POST':
        # get the current PlotParams model instance
        plot_params = PlotParams.objects.filter(input_data=input_instance_id).get()
        # connect InputDataPlottingContext form to current PlotParams model
        form = InputDataPlottingContext(request.POST,instance=plot_params)
        # if form is valid, re-render page which will result in replotting
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
    # if any valid GO terms found, perform plotting
    if nrows>0:
        plotimg = makeplot.make_barplot(frame,x_t_s,y_t_s,x_l_s,y_l_s,p_t_s,f_x_s,f_y_s,cmap)
    # if no valid GO terms found, perform plotting, set plot image to None
    else:
        plotimg = None
    # re-create form with current parameters to show
    form = InputDataPlottingContext(instance=plot_params)
    
   # add plot image and updated form to context
    context['plot'] = plotimg
    context['form'] = form
    
    # plot the data
    return render(request, 'enrichmentplotter/data_plotted.html',context=context)
    
    
class UploadNameDictionary(PermissionRequiredMixin,CreateView):
    '''
    defines initial NameDictionary model creation
    probably should be done only once by admin, so no actaul upload page is yet defined
    '''
    permission_required = 'enrichmentplotter.modify_data'
    model = NameDictionary
    fields = ['base_table_name', 'base_table']
    


    def model_created_or_updated(sender, created=False, **kwargs):
        # get current instance and path to csv file
        table = kwargs['instance']
        raw_path = table.base_table.path
        # create name dictionary
        name_dictionary = gaf_transform.get_namedict(raw_path)
        # update file path field and actual dictionary as JSON
        NameDictionary.objects.filter(pk=table.pk).update(base_table_path = raw_path)
        NameDictionary.objects.filter(pk=table.pk).update(name_dict = name_dictionary)

    post_save.connect(model_created_or_updated, sender=NameDictionary)
    # probably will never be implemented as page
    success_url = reverse_lazy('index')






    
    