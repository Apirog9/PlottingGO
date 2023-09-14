# Generated by Django 4.2.2 on 2023-09-09 10:16

from django.db import migrations, models
import django.db.models.deletion
import uuid


class Migration(migrations.Migration):
    dependencies = [
        ("enrichmentplotter", "0020_namedictionary"),
    ]

    operations = [
        migrations.CreateModel(
            name="PlotParams",
            fields=[
                (
                    "params_id",
                    models.UUIDField(
                        default=uuid.uuid4, primary_key=True, serialize=False
                    ),
                ),
                (
                    "colormap",
                    models.CharField(
                        choices=[("viridis", "viridis"), ("plasma", "plasma")],
                        default="viridis",
                        max_length=100,
                    ),
                ),
                (
                    "input_data",
                    models.ForeignKey(
                        null=True,
                        on_delete=django.db.models.deletion.CASCADE,
                        to="enrichmentplotter.listenrichmentinput",
                    ),
                ),
            ],
        ),
    ]