# Generated by Django 4.2.2 on 2023-09-07 18:07

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        (
            "enrichmentplotter",
            "0019_rename_background_items_listenrichmentinput_background_items_annotation",
        ),
    ]

    operations = [
        migrations.CreateModel(
            name="NameDictionary",
            fields=[
                (
                    "base_table_name",
                    models.CharField(max_length=200, primary_key=True, serialize=False),
                ),
                ("base_table", models.FileField(max_length=254, upload_to="GO_tables")),
                (
                    "base_table_path",
                    models.CharField(
                        blank=True, default="", max_length=1000, null=True
                    ),
                ),
                ("name_dict", models.JSONField(blank=True, null=True)),
            ],
        ),
    ]
