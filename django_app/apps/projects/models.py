"""User-scoped project storage.

Mirrors the schema of the Streamlit `projects` table so migrating existing
users' data is a straight copy (same field names, same semantics).
"""
from django.conf import settings
from django.db import models


class Project(models.Model):
    """A saved item belonging to one user.

    `data` is a JSONField — any G-Synth module can persist arbitrary
    result payloads without schema migrations. `sequence` and `notes`
    are surfaced separately because they're queried on the list view.
    """

    # What the app actually saves. Kept in step with the views: a module
    # missing from here still stores, but shows as an unlabelled chip and
    # fails model validation, so the two drift silently.
    MODULES = (
        ("general",              "General"),
        ("ssd",                  "Small Sequence Design"),
        ("merzoug_assembly",     "Merzoug Assembly"),
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
