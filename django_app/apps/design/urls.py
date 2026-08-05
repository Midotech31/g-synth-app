"""Design routes — mounted under /api/design/ in config.urls."""
from django.urls import path

from apps.design.views import (
    AlignView,
    CloneExportView,
    CloneView,
    ConstructExportView,
    EnzymeCatalogueView,
    LigationView,
    MerzougAssemblyView,
    OptimiseView,
    OrderSheetView,
    PrimerExportView,
    ProtocolView,
    SequencingPrimerView,
    SSDDesignView,
    VectorCatalogueView,
    VectorSequenceView,
    VerifyView,
)

urlpatterns = [
    path("enzymes/", EnzymeCatalogueView.as_view(), name="design-enzymes"),
    path("ssd/", SSDDesignView.as_view(), name="design-ssd"),
    path("assembly/", MerzougAssemblyView.as_view(), name="design-assembly"),
    path("assembly/order-sheet/", OrderSheetView.as_view(), name="design-order-sheet"),
    path("assembly/protocol/", ProtocolView.as_view(), name="design-protocol"),
    path("optimise/", OptimiseView.as_view(), name="design-optimise"),
    path("ligation/", LigationView.as_view(), name="design-ligation"),
    path("primers/", SequencingPrimerView.as_view(), name="design-primers"),
    path("primers/export/", PrimerExportView.as_view(), name="design-primers-export"),
    path("align/", AlignView.as_view(), name="design-align"),
    path("verify/", VerifyView.as_view(), name="design-verify"),
    path("clone/", CloneView.as_view(), name="design-clone"),
    path("clone/export/", CloneExportView.as_view(), name="design-clone-export"),
    path("assembly/export/", ConstructExportView.as_view(), name="design-construct-export"),
    path("vectors/", VectorCatalogueView.as_view(), name="design-vectors"),
    path("vectors/<str:key>/", VectorSequenceView.as_view(), name="design-vector"),
]
