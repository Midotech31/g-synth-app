"""Design routes — mounted under /api/design/ in config.urls."""
from django.urls import path

from apps.design.views import (
    EnzymeCatalogueView,
    MerzougAssemblyView,
    OrderSheetView,
    ProtocolView,
    SSDDesignView,
)

urlpatterns = [
    path("enzymes/", EnzymeCatalogueView.as_view(), name="design-enzymes"),
    path("ssd/", SSDDesignView.as_view(), name="design-ssd"),
    path("assembly/", MerzougAssemblyView.as_view(), name="design-assembly"),
    path("assembly/order-sheet/", OrderSheetView.as_view(), name="design-order-sheet"),
    path("assembly/protocol/", ProtocolView.as_view(), name="design-protocol"),
]
