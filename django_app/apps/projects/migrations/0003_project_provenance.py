from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [("projects", "0002_alter_project_module")]

    operations = [
        migrations.AddField(
            model_name="project",
            name="provenance",
            field=models.JSONField(blank=True, default=dict, editable=False),
        ),
    ]
