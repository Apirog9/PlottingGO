from django.contrib import admin
from .models import AnnotationTable
from .models import NameDictionary

# Register your models here.
admin.site.register(AnnotationTable)
admin.site.register(NameDictionary)


