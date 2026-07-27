from django.contrib import admin

from apps.projects.models import Project


@admin.register(Project)
class ProjectAdmin(admin.ModelAdmin):
    list_display = ("name", "user", "module", "updated_at")
    list_filter = ("module", "created_at")
    search_fields = ("name", "user__email", "sequence")
    readonly_fields = ("created_at", "updated_at")
