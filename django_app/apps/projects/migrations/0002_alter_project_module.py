from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("projects", "0001_initial"),
    ]

    operations = [
        migrations.AlterField(
            model_name="project",
            name="module",
            field=models.CharField(
                choices=[
                    ("general", "General"),
                    ("ssd", "Small Sequence Design"),
                    ("extended_sequence_design", "Extended Sequence Design"),
                    ("cloning", "Cloning"),
                    ("codon_optimization", "Codon Optimisation"),
                    ("plasmid_visualizer", "Imported Sequence"),
                    ("primer_generator", "Primer Generator"),
                    ("alignment_tools", "Alignment Tools"),
                ],
                default="general",
                max_length=64,
            ),
        ),
    ]
