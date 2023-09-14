from django.urls import path
from . import views


urlpatterns = [
    path('', views.index, name='index'),
    path('upload/gaf/', views.UploadGaf.as_view(), name='upload_gaf'),
    path('upload/gaf/', views.UploadGaf.as_view(), name='upload_gaf'),
    path('upload/input/', views.input_data, name='input_data'),
    path('upload/<uuid:input_instance_id>/ready/', views.data_ready, name='data_ready_to_plot'),
    path('upload/<uuid:input_instance_id>/plotted/', views.data_plotted, name='data_plotted')
]