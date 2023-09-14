from django import forms
from django.forms import ModelForm
from enrichmentplotter.models import ListEnrichmentInput
from enrichmentplotter.models import PlotParams
from django.core.exceptions import ValidationError



#accept text input or TODO text file
#parse text inputs to lists
#write lists in local text file
#perform annotation
#report number of ids found and missing
#here or in view
class InputData(ModelForm):

    def clean(self):
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

    class Meta:
        model = PlotParams
        fields = ['max_p_value','min_protein_number','min_fraction','colormap','font_size_x_labels','font_size_y_labels','font_size_x_title','font_size_y_title','font_size_plot_title','figure_size_x','figure_size_y']
