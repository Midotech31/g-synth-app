from django.conf import settings
from django.db import models


class Project(models.Model):


    MODULES = (
        ("general",              "General"),
        ("ssd",                  "Small Sequence Design"),
        ("extended_sequence_design",     "Extended Sequence Design"),
        ("cloning",              "Cloning"),
        ("codon_optimization",   "Codon Optimisation"),
        ("plasmid_visualizer",   "Imported Sequence"),
        ("primer_generator",     "Primer Generator"),
        ("alignment_tools",      "Alignment Tools"),
    )

    user = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        related_name="projects",
    )
    name = models.CharField(max_length=255)
    module = models.CharField(max_length=64, choices=MODULES, default="general")
    sequence = models.TextField(blank=True, default="")
    notes = models.TextField(blank=True, default="")
    data = models.JSONField(default=dict, blank=True)
    provenance = models.JSONField(default=dict, blank=True, editable=False)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "gsynth_project"
        ordering = ("-updated_at",)
        indexes = [
            models.Index(fields=["user", "-updated_at"]),
            models.Index(fields=["user", "module"]),
        ]

    def __str__(self) -> str:
        return f"{self.name} ({self.module}) — {self.user.email}"
