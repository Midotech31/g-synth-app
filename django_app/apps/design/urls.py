"""Design routes — mounted under /api/design/ in config.urls."""
from django.urls import path

from apps.design.views import (
    CloneView,
    EnzymeCatalogueView,
    MerzougAssemblyView,
    OrderSheetView,
    ProtocolView,
    SSDDesignView,
    VectorCatalogueView,
    VectorSequenceView,
)

urlpatterns = [
    path("enzymes/", EnzymeCatalogueView.as_view(), name="design-enzymes"),
    path("ssd/", SSDDesignView.as_view(), name="design-ssd"),
    path("assembly/", MerzougAssemblyView.as_view(), name="design-assembly"),
    path("assembly/order-sheet/", OrderSheetView.as_view(), name="design-order-sheet"),
    path("assembly/protocol/", ProtocolView.as_view(), name="design-protocol"),
    path("clone/", CloneView.as_view(), name="design-clone"),
    path("vectors/", VectorCatalogueView.as_view(), name="design-vectors"),
    path("vectors/<str:key>/", VectorSequenceView.as_view(), name="design-vector"),
]
