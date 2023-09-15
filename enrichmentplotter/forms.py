from django import forms
from django.forms import ModelForm
from enrichmentplotter.models import ListEnrichmentInput
from enrichmentplotter.models import PlotParams
from django.core.exceptions import ValidationError




class InputData(ModelForm):
    '''
    Form used for data input for ListEnrichmentInput model
    id_list - unlimited field containing one accepted ID per line
    background_list - unlimited field containing one accepted ID per line
    GO_set - gene ontology set, should match organism of the list&background_list
    input type - select Uniprt Id or Gene
    GO_BP,GO_MF,GO_CC - select at least one GO class to use on plot. 
    '''


    def clean(self):
        # check if at least one ontology class is selected
        cleaned_data = super().clean()
        BP = cleaned_data.get('GO_BP')
        MF = cleaned_data.get('GO_MF')
        CC = cleaned_data.get('GO_CC')
        if not any([BP,CC,MF]):
            raise ValidationError("Select some type of ontology")
        
    class Meta:
        model = ListEnrichmentInput
        fields = ['id_list','background_list','GO_set','input_type','GO_BP','GO_MF','GO_CC']
        
        
class InputDataPlottingContext(ModelForm):
    '''
    Form used to store and define plotting parameters for ProtParams model
    max_p_value - maximum CORRECTED p-value for a term to be plotted
    min_protein_number - minimum query+background protein/gene number for a term to be plotted
    min_fraction -  minimum fraction of query proteins in all proteins for a term to be plotted
    colormap - matplotlib colormap to reflect CORRECTED p value
    font_size_x_labels,font_size_y_labels,font_size_x_title,font_size_y_title,font_size_plot_title -font sizes for plot
    figure_size_x,figure_size_y - figure size in inches
    
    '''

    class Meta:
        model = PlotParams
        fields = ['max_p_value','min_protein_number','min_fraction','colormap','font_size_x_labels','font_size_y_labels','font_size_x_title','font_size_y_title','font_size_plot_title','figure_size_x','figure_size_y']
