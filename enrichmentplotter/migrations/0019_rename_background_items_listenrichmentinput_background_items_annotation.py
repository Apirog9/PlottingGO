# Generated by Django 4.2.2 on 2023-09-06 18:29

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ("enrichmentplotter", "0018_alter_listenrichmentinput_go_cc_and_more"),
    ]

    operations = [
        migrations.RenameField(
            model_name="listenrichmentinput",
            old_name="background_items",
            new_name="background_items_annotation",
        ),
    ]