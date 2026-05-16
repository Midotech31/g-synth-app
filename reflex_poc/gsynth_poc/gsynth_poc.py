"""Main Reflex app — wires pages, theme, and routes."""
from __future__ import annotations
import reflex as rx

from gsynth_poc.pages.crispr import crispr_page
from gsynth_poc.pages.dashboard import dashboard_page
from gsynth_poc.pages.index import index_page
from gsynth_poc.pages.login import login_page
from gsynth_poc.pages.project_detail import project_detail_page
from gsynth_poc.pages.signup import signup_page


app = rx.App(
    theme=rx.theme(appearance="light", accent_color="cyan", radius="medium"),
    stylesheets=[],
)

app.add_page(index_page, route="/", title="G-Synth POC")
app.add_page(login_page, route="/login", title="Sign in — G-Synth")
app.add_page(signup_page, route="/signup", title="Create account — G-Synth")
app.add_page(dashboard_page, route="/dashboard", title="Dashboard — G-Synth")
app.add_page(crispr_page, route="/crispr", title="CRISPR designer — G-Synth")
app.add_page(project_detail_page, route="/projects/[project_id]",
             title="Project — G-Synth")
