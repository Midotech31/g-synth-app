"""Reflex configuration for the G-Synth POC."""
import reflex as rx

config = rx.Config(
    app_name="gsynth_poc",
    db_url="sqlite:///gsynth_poc.db",
    show_built_with_reflex=False,
)
