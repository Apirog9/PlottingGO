# Generated by Django 4.2.2 on 2023-09-15 09:09

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        (
            "enrichmentplotter",
            "0022_plotparams_figure_size_x_plotparams_figure_size_y_and_more",
        ),
    ]

    operations = [
        migrations.AlterField(
            model_name="plotparams",
            name="max_p_value",
            field=models.FloatField(default=0.05),
        ),
    ]
