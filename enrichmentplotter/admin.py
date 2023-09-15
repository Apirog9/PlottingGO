from django.contrib import admin
from .models import AnnotationTable
from .models import NameDictionary

# Register your models here.
# Enable Annotation Table files to be added and transformed from admin panel
admin.site.register(AnnotationTable)
# Enable Name dictionary file to be uploaded and transformed from admin panel
# Only one instance of the model should exis, its contensts will be updated within \utilities\gaf_transform.py
admin.site.register(NameDictionary)


