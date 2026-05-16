# modules/settings.py

"""
Settings & Session Module
- JSON-based settings
- Recent files list
- Dark-mode toggle
- Logging settings (handled globally)
"""

import streamlit as st
import json
import os

SETTINGS_FILE = "settings.json"
default_settings = {
    "dark_mode": False,
    "last_used_directory": "",
    "recent_files": []
}

def main():
    st.title("Settings & Session")

    # Load or initialize
    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, "r") as f:
                current = json.load(f)
        except:
            st.warning("Error loading settings; using defaults.")
            current = default_settings.copy()
    else:
        current = default_settings.copy()

    # Dark mode toggle
    dark = st.checkbox("Dark Mode", value=current.get("dark_mode", False))
    current["dark_mode"] = dark

    # Recent files list (display only)
    st.write("### Recent Files")
    recents = current.get("recent_files", [])
    if recents:
        for f in recents:
            st.write(f"- {f}")
    else:
        st.write("No recent files.")

    # Save settings
    if st.button("Save Settings"):
        try:
            with open(SETTINGS_FILE, "w") as f:
                json.dump(current, f, indent=2)
            st.success("Settings saved.")
        except Exception as e:
            st.error(f"Error saving settings: {e}")
