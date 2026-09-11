from django.db import migrations

LEGACY_MODULE = "merzoug_assembly"
ESD_MODULE = "extended_sequence_design"


def normalize_esd_module(apps, schema_editor):
    project = apps.get_model("projects", "Project")
    project.objects.filter(module=LEGACY_MODULE).update(module=ESD_MODULE)


def reverse_module_normalization(apps, schema_editor):
    project = apps.get_model("projects", "Project")
    project.objects.filter(module=ESD_MODULE).update(module=LEGACY_MODULE)


class Migration(migrations.Migration):
    dependencies = [("projects", "0003_project_provenance")]

    operations = [
        migrations.RunPython(normalize_esd_module, reverse_module_normalization),
    ]
