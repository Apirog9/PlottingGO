from django.urls import path
from . import views

# defines index view for welcome page
# defines GAF file input page, currenly experimental
# defines input data selection page
# defines data ready page, showing soe annotation performance and prompting to plot
# defines data plotted page, containing actual plot, and form for updating plot parameters
urlpatterns = [
    path('', views.index, name='index'),
    path('upload/gaf/', views.UploadGaf.as_view(), name='upload_gaf'),
    path('upload/input/', views.input_data, name='input_data'),
    path('upload/<uuid:input_instance_id>/ready/', views.data_ready, name='data_ready_to_plot'),
    path('upload/<uuid:input_instance_id>/plotted/', views.data_plotted, name='data_plotted')
]