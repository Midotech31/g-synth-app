# -*- coding: utf-8 -*-
"""
app.py

Main entry point for the G-Synth Streamlit application.
Enhanced with modern CSS, interactive sidebar, better cards, and workflow ordering.
FIXED: Sidebar text legibility and removed purple banners from module pages (except Home).
"""

import json
import logging
from pathlib import Path

import streamlit as st
from streamlit_option_menu import option_menu

# ───────────────────────────────────────────────────────────────────────────────
# 1. LOGGING CONFIGURATION
# ───────────────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    filename="g-synth.log",
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("G-Synth:APP")

# ───────────────────────────────────────────────────────────────────────────────
# 2. CONSTANTS & MAPPINGS (REORDERED FOR WORKFLOW)
# ───────────────────────────────────────────────────────────────────────────────
PAGE_TITLE = "G-Synth 2025.2.0"
SETTINGS_PATH = Path("settings.json")
CSS_PATH = Path("assets/styles.css")
LOGO_PATH = "assets/logo.png"

# WORKFLOW-ORDERED MODULE MAPPING
MODULE_MAPPING = {
    # Core Workflow (Main Pipeline)
    "🏠 Home": "home",
    "🔬 Small Sequence Design": "ssd",
    "🧬 Translation Tools": "translation",
    "⚙️ Codon Optimization": "codon_optimization",
    "🔗 Extended Synthesis": "extended_synthesis",
    "🧪 Hybridization": "hybridization",
    "🔗 Ligation Check": "ligation_check",
    "🎯 Primer Generator": "primer_generator",
    "🔄 Reverse Complement": "reverse_complement",
    
    # Advanced Analysis Tools
    "✂️ CRISPR sgRNA Designer": "crispr_designer",
    "🧫 Plasmid Visualizer": "plasmid_visualizer",
    "🧩 Ligation Calculator": "ligase_calculator",
    "📊 Alignment Tools": "alignment_tools",
    
    # AI/Computational Tools
    "🖥️ In Silico Docking": "docking_module",
    "🧠 Functional Prediction": "functional_prediction",
    
    # System Tools
    "⚡ Settings": "settings",
    "📖 Help & Guide": "help_guide",
}

# ───────────────────────────────────────────────────────────────────────────────
# 3. ENHANCED HELPER FUNCTIONS
# ───────────────────────────────────────────────────────────────────────────────
def inject_custom_css():
    """
    FIXED CSS injection with maximum specificity for sidebar buttons.
    """
    # Load external CSS file FIRST
    if CSS_PATH.exists():
        try:
            with open(CSS_PATH, "r", encoding="utf-8") as f:
                css = f.read()
            # Apply base CSS first
            st.markdown(f"<style>{css}</style>", unsafe_allow_html=True)
        except Exception as e:
            st.warning(f"Failed to load custom CSS: {e}")
            logger.error(f"CSS injection error: {e}")
    
    # CRITICAL FIX: Apply sidebar fixes AFTER base CSS with MAXIMUM specificity
    sidebar_fix_css = """
    /* ═══════════════════════════════════════════════════════════════════════════════
       CRITICAL SIDEBAR TEXT FIX - MAXIMUM SPECIFICITY AND NUCLEAR OPTION
    ═══════════════════════════════════════════════════════════════════════════════ */
    
    /* Enhanced sidebar background */
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%) !important;
        border-right: 3px solid #ff4b4b !important;
        box-shadow: 4px 0 20px rgba(0,0,0,0.3) !important;
    }
    
    [data-testid="stSidebar"] > div {
        background: transparent !important;
        padding-top: 1rem !important;
    }
    
    /* NUCLEAR OPTION 1: Maximum specificity sidebar button styling */
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] div[data-testid="stSidebarContent"] div .stButton > button {
        color: #000000 !important;
        font-weight: 700 !important;
        background: rgba(255, 255, 255, 0.98) !important;
        border: 2px solid rgba(255, 255, 255, 0.4) !important;
        border-radius: 12px !important;
        padding: 0.75rem 1.25rem !important;
        margin: 0.25rem 0.5rem !important;
        width: calc(100% - 1rem) !important;
        text-align: left !important;
        font-size: 0.95rem !important;
        line-height: 1.4 !important;
        box-shadow: 0 2px 8px rgba(0,0,0,0.15) !important;
        transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1) !important;
        cursor: pointer !important;
        position: relative !important;
        overflow: hidden !important;
        display: flex !important;
        align-items: center !important;
        justify-content: flex-start !important;
        font-family: 'Nunito Sans', sans-serif !important;
    }
    
    /* NUCLEAR OPTION 2: Force ALL text elements inside buttons to be black and bold */
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] div[data-testid="stSidebarContent"] div .stButton > button p,
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] div[data-testid="stSidebarContent"] div .stButton > button span,
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] div[data-testid="stSidebarContent"] div .stButton > button div,
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] div[data-testid="stSidebarContent"] div .stButton > button *,
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] .stButton button p,
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] .stButton button span,
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] .stButton button div,
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] .stButton button * {
        color: #000000 !important;
        font-weight: 700 !important;
        text-shadow: none !important;
        font-family: 'Nunito Sans', sans-serif !important;
    }
    
    /* NUCLEAR OPTION 3: Alternative selectors for different HTML structures */
    [data-testid="stSidebar"] .element-container .stButton button {
        color: #000000 !important;
        font-weight: 700 !important;
        background: rgba(255, 255, 255, 0.98) !important;
    }
    
    [data-testid="stSidebar"] button[kind="secondary"] {
        color: #000000 !important;
        background-color: rgba(255, 255, 255, 0.98) !important;
        font-weight: 700 !important;
    }
    
    [data-testid="stSidebar"] button[kind="primary"] {
        color: #000000 !important;
        background-color: rgba(255, 255, 255, 0.98) !important;
        font-weight: 700 !important;
    }
    
    [data-testid="stSidebar"] button {
        color: #000000 !important;
        background: rgba(255, 255, 255, 0.98) !important;
        font-weight: 700 !important;
    }
    
    /* NUCLEAR OPTION 4: Universal selector for any text in sidebar buttons */
    [data-testid="stSidebar"] * button * {
        color: #000000 !important;
        font-weight: 700 !important;
    }
    
    /* Hover effects with maximum specificity */
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] div[data-testid="stSidebarContent"] div .stButton > button:hover {
        color: #ffffff !important;
        background: linear-gradient(135deg, #ff4b4b, #ff6b6b) !important;
        border-color: #ff4b4b !important;
        transform: translateX(8px) scale(1.02) !important;
        box-shadow: 0 8px 25px rgba(255, 75, 75, 0.4) !important;
    }
    
    /* Hover text with maximum specificity */
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] div[data-testid="stSidebarContent"] div .stButton > button:hover p,
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] div[data-testid="stSidebarContent"] div .stButton > button:hover span,
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] div[data-testid="stSidebarContent"] div .stButton > button:hover div,
    html body div[data-testid="stApp"] div[data-testid="stSidebar"] div[data-testid="stSidebarContent"] div .stButton > button:hover *,
    [data-testid="stSidebar"] .stButton > button:hover p,
    [data-testid="stSidebar"] .stButton > button:hover span,
    [data-testid="stSidebar"] .stButton > button:hover div,
    [data-testid="stSidebar"] .stButton > button:hover * {
        color: #ffffff !important;
        text-shadow: 0 1px 2px rgba(0,0,0,0.1) !important;
    }
    
    /* Workflow section styling */
    .workflow-section {
        background: linear-gradient(135deg, #ff4444 0%, #cc0000 100%) !important;
        color: white !important;
        padding: 1rem 1.5rem !important;
        border-radius: 12px !important;
        margin: 1.5rem 0 1rem 0 !important;
        box-shadow: 0 6px 20px rgba(255, 68, 68, 0.3) !important;
        text-align: center !important;
    }
    
    .workflow-section h3 {
        margin: 0 !important;
        font-weight: 700 !important;
        text-shadow: 0 2px 4px rgba(0,0,0,0.2) !important;
        color: #ffffff !important;
    }
    
    /* Enhanced card styling */
    .feature-card, .card, .floating-card {
        background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%) !important;
        border: 2px solid #e2e8f0 !important;
        border-radius: 16px !important;
        padding: 1.5rem !important;
        margin: 0.75rem 0 !important;
        box-shadow: 0 4px 15px rgba(0,0,0,0.08) !important;
        transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1) !important;
        position: relative !important;
        overflow: hidden !important;
        cursor: pointer !important;
    }
    
    .feature-card:hover, .card:hover, .floating-card:hover {
        background: linear-gradient(135deg, #ff4b4b, #ff6b6b) !important;
        color: #ffffff !important;
        border-color: #ff4b4b !important;
        box-shadow: 0 12px 40px rgba(255, 75, 75, 0.4) !important;
        transform: translateY(-8px) scale(1.02) !important;
    }
    
    /* Status indicator styling */
    .status-indicator {
        position: fixed !important;
        bottom: 1.5rem !important;
        left: 1rem !important;
        background: linear-gradient(145deg, #22c55e, #16a34a) !important;
        color: white !important;
        padding: 0.5rem 1rem !important;
        border-radius: 20px !important;
        font-size: 0.85rem !important;
        font-weight: 600 !important;
        box-shadow: 0 4px 15px rgba(34, 197, 94, 0.4) !important;
        z-index: 1000 !important;
        border: 2px solid rgba(255, 255, 255, 0.2) !important;
        backdrop-filter: blur(10px) !important;
        animation: pulse 2s infinite !important;
    }
    
    @keyframes pulse {
        0%, 100% { opacity: 1; }
        50% { opacity: 0.7; }
    }
    """
    
    # Apply sidebar fixes with highest priority
    st.markdown(f"<style>{sidebar_fix_css}</style>", unsafe_allow_html=True)
    
    # ULTIMATE NUCLEAR OPTION: JavaScript force-apply if CSS still fails
    javascript_fix = """
    <script>
    function forceSidebarStyles() {
        // Wait for Streamlit to render
        setTimeout(function() {
            const sidebarButtons = document.querySelectorAll('[data-testid="stSidebar"] button');
            console.log('Found sidebar buttons:', sidebarButtons.length);
            
            sidebarButtons.forEach((btn, index) => {
                // Force styles directly via JavaScript
                btn.style.setProperty('color', '#000000', 'important');
                btn.style.setProperty('font-weight', '700', 'important');
                btn.style.setProperty('background', 'rgba(255, 255, 255, 0.98)', 'important');
                btn.style.setProperty('border', '2px solid rgba(255, 255, 255, 0.4)', 'important');
                btn.style.setProperty('border-radius', '12px', 'important');
                btn.style.setProperty('padding', '0.75rem 1.25rem', 'important');
                btn.style.setProperty('margin', '0.25rem 0.5rem', 'important');
                btn.style.setProperty('width', 'calc(100% - 1rem)', 'important');
                btn.style.setProperty('text-align', 'left', 'important');
                btn.style.setProperty('font-size', '0.95rem', 'important');
                btn.style.setProperty('line-height', '1.4', 'important');
                btn.style.setProperty('box-shadow', '0 2px 8px rgba(0,0,0,0.15)', 'important');
                btn.style.setProperty('cursor', 'pointer', 'important');
                btn.style.setProperty('font-family', 'Nunito Sans, sans-serif', 'important');
                
                // Force styles on all text elements inside button
                const textElements = btn.querySelectorAll('*');
                textElements.forEach(el => {
                    el.style.setProperty('color', '#000000', 'important');
                    el.style.setProperty('font-weight', '700', 'important');
                    el.style.setProperty('text-shadow', 'none', 'important');
                    el.style.setProperty('font-family', 'Nunito Sans, sans-serif', 'important');
                });
                
                console.log(`Styled button ${index}:`, {
                    color: btn.style.color,
                    backgroundColor: btn.style.background,
                    fontWeight: btn.style.fontWeight
                });
            });
        }, 500);
        
        // Repeat after 2 seconds to catch any late-loading elements
        setTimeout(function() {
            forceSidebarStyles();
        }, 2000);
    }
    
    // Apply styles immediately and on any DOM changes
    forceSidebarStyles();
    
    // Watch for DOM changes and reapply styles
    const observer = new MutationObserver(function(mutations) {
        mutations.forEach(function(mutation) {
            if (mutation.type === 'childList') {
                forceSidebarStyles();
            }
        });
    });
    
    // Start observing
    if (document.querySelector('[data-testid="stSidebar"]')) {
        observer.observe(document.querySelector('[data-testid="stSidebar"]'), {
            childList: true,
            subtree: true
        });
    }
    </script>
    """
    
    # Apply JavaScript fix
    st.markdown(javascript_fix, unsafe_allow_html=True)

def load_settings():
    """
    Loads JSON settings from settings.json, returning a dict.
    If loading fails, logs the error and returns an empty dict.
    """
    if SETTINGS_PATH.exists():
        try:
            with open(SETTINGS_PATH, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            st.warning(f"Error loading settings; using defaults: {e}")
            logger.error(f"Settings load error: {e}")
            return {}
    return {}

def load_module(module_key: str):
    """
    Dynamically imports a module from the 'modules' package by its key.
    Returns the imported module or None if import fails.
    """
    try:
        module = __import__(f"modules.{module_key}", fromlist=[module_key])
        return module
    except UnicodeDecodeError as e:
        st.error(f"Encoding error in module '{module_key}': {e}")
        st.warning(f"Please check that modules/{module_key}.py is saved with UTF-8 encoding")
        logger.error(f"Unicode decode error for modules.{module_key}: {e}")
        return None
    except Exception as e:
        st.error(f"Error loading module '{module_key}': {e}")
        logger.error(f"Import failure for modules.{module_key}: {e}")
        return None

def display_enhanced_sidebar():
    """
    Enhanced sidebar with workflow organization and modern design.
    Returns the selected menu option as a string.
    """
    with st.sidebar:
        # Enhanced logo section
        if Path(LOGO_PATH).exists():
            st.image(LOGO_PATH, use_container_width=True)
        else:
            st.markdown("""
            <div style="
                text-align: center; 
                padding: 1.5rem 0; 
                margin-bottom: 2rem;
                background: linear-gradient(135deg, rgba(255, 75, 75, 0.2), rgba(255, 107, 107, 0.2));
                border-radius: 16px;
                border: 2px solid rgba(255, 255, 255, 0.2);
            ">
                <h2 style="
                    color: #ffffff; 
                    margin: 0; 
                    font-size: 1.8rem; 
                    font-weight: 800;
                    text-shadow: 0 2px 4px rgba(0,0,0,0.3);
                ">🧬 G-Synth</h2>
                <p style="
                    color: #ffffff; 
                    margin: 0.5rem 0 0 0; 
                    font-size: 1rem; 
                    opacity: 0.9;
                    font-weight: 500;
                ">Gene Synthesis Toolkit</p>
            </div>
            """, unsafe_allow_html=True)

        # Workflow sections
        workflow_sections = {
            "Core Workflow": [
                "🏠 Home", "🔬 Small Sequence Design", "🧬 Translation Tools", 
                "⚙️ Codon Optimization", "🔗 Extended Synthesis", "🧪 Hybridization", "🔗 Ligation Check"
            ],
            "Design Tools": [
                "🎯 Primer Generator", "🔄 Reverse Complement"
            ],
            "Advanced Analysis": [
                "✂️ CRISPR sgRNA Designer", "🧫 Plasmid Visualizer", 
                "🧩 Ligation Calculator", "📊 Alignment Tools"
            ],
            "AI Tools": [
                "🖥️ In Silico Docking", "🧠 Functional Prediction"
            ],
            "System": [
                "⚡ Settings", "📖 Help & Guide"
            ]
        }

        # Create expandable sections
        selected = None
        
        for section_name, modules in workflow_sections.items():
            st.markdown(f"""
            <div class="workflow-section">
                <h3>{section_name}</h3>
            </div>
            """, unsafe_allow_html=True)
            
            # Create buttons for each module in the section
            for module in modules:
                if module in MODULE_MAPPING:
                    # Create a styled button for each module
                    if st.button(
                        module, 
                        key=f"btn_{module}",
                        use_container_width=True,
                        help=f"Navigate to {module.split(' ', 1)[-1]}"
                    ):
                        selected = module
                        st.session_state.selected_module = module
        
        # Get the current selection
        if not selected and 'selected_module' in st.session_state:
            selected = st.session_state.selected_module
        elif not selected:
            selected = "🏠 Home"
            st.session_state.selected_module = selected
        
        # Enhanced status indicator
        st.markdown("""
        <div class="status-indicator">
            🟢 System Online
        </div>
        """, unsafe_allow_html=True)
        
        # Module information
        if selected and selected != "🏠 Home":
            module_name = selected.split(' ', 1)[-1]
            st.markdown(f"""
            <div style="
                margin-top: 2rem;
                padding: 1rem;
                background: rgba(255, 255, 255, 0.1);
                border-radius: 12px;
                border: 1px solid rgba(255, 255, 255, 0.2);
            ">
                <h4 style="color: #ffffff; margin: 0 0 0.5rem 0; font-size: 0.9rem;">
                    Current Module:
                </h4>
                <p style="color: #ffffff; margin: 0; font-size: 0.8rem; opacity: 0.9;">
                    {module_name}
                </p>
            </div>
            """, unsafe_allow_html=True)
        
    return selected

def render_enhanced_header(selected: str):
    """
    Enhanced page header with modern design.
    FIXED: Remove purple banners from module pages (except Home).
    """
    # ONLY show the enhanced header/banner on the Home page
    # All other modules get clean interfaces without purple banners
    if selected == "🏠 Home":
        # Home page handles its own welcome banner in home.py
        # No header needed here as home.py has its own design
        pass
    # For all other modules: NO PURPLE BANNER
    # This removes the purple gradient headers from module pages as requested

def run_selected_module(selected: str):
    """
    Enhanced module runner with better error handling.
    """
    module_key = MODULE_MAPPING.get(selected)
    if not module_key:
        st.error("Invalid selection.")
        return

    module = load_module(module_key)
    if not module:
        return

    # Add loading indicator
    with st.spinner(f"Loading {selected.split(' ', 1)[-1]}..."):
        if hasattr(module, "main"):
            try:
                module.main()
            except Exception as e:
                st.error(f"Error running '{selected}': {e}")
                logger.error(f"Runtime error in modules.{module_key}.main(): {e}")
                
                # Enhanced error display
                with st.expander("🔍 Debug Information", expanded=False):
                    st.exception(e)
                    st.code(f"""
Module: {module_key}
Function: main()
Error: {str(e)}
                    """)
                    
        elif hasattr(module, "app"):
            try:
                module.app()
            except Exception as e:
                st.error(f"Error running '{selected}': {e}")
                logger.error(f"Runtime error in modules.{module_key}.app(): {e}")
                
                # Enhanced error display
                with st.expander("🔍 Debug Information", expanded=False):
                    st.exception(e)
                    st.code(f"""
Module: {module_key}
Function: app()
Error: {str(e)}
                    """)
        else:
            st.warning(f"Module '{selected}' has no entrypoint (main/app).")
            st.info("Please ensure your module has either a main() or app() function.")

def render_enhanced_footer():
    """
    Enhanced footer with modern design.
    """
    st.markdown("---")
    st.markdown("""
    <div style="
        text-align: center; 
        font-size: 0.9rem; 
        color: #666; 
        padding: 2.5rem 2rem;
        background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
        border-radius: 16px;
        margin-top: 3rem;
        border: 2px solid #e2e8f0;
        box-shadow: 0 4px 15px rgba(0,0,0,0.05);
    ">
        <div style="margin-bottom: 1rem;">
            <h4 style="
                margin: 0 0 0.5rem 0; 
                font-weight: 700; 
                color: #374151;
                font-size: 1.1rem;
            ">
                🧬 G-Synth Toolkit
            </h4>
            <p style="
                margin: 0; 
                font-weight: 600;
                color: #6b7280;
            ">
                Advanced Gene Synthesis Platform for Therapeutic Peptides
            </p>
        </div>
        <div style="
            display: flex; 
            justify-content: center; 
            gap: 2rem; 
            margin: 1.5rem 0;
            flex-wrap: wrap;
        ">
            <span style="
                background: linear-gradient(135deg, #667eea, #764ba2); 
                color: white; 
                padding: 0.4rem 0.8rem; 
                border-radius: 20px; 
                font-size: 0.8rem;
                font-weight: 600;
            ">Version 2025.2.0</span>
            <span style="
                background: linear-gradient(135deg, #ff4b4b, #ff6b6b); 
                color: white; 
                padding: 0.4rem 0.8rem; 
                border-radius: 20px; 
                font-size: 0.8rem;
                font-weight: 600;
            ">Powered by Streamlit</span>
            <span style="
                background: linear-gradient(135deg, #10b981, #059669); 
                color: white; 
                padding: 0.4rem 0.8rem; 
                border-radius: 20px; 
                font-size: 0.8rem;
                font-weight: 600;
            ">Workflow Optimized</span>
        </div>
        <p style="
            margin: 0.5rem 0 0 0; 
            font-size: 0.8rem; 
            opacity: 0.8;
            color: #9ca3af;
        ">
            © 2025 Dr. Mohamed Merzoug • Genomics Technology Platform • Higher School of Biological Sciences of Oran
        </p>
    </div>
    """, unsafe_allow_html=True)

# ───────────────────────────────────────────────────────────────────────────────
# 4. MAIN APPLICATION
# ───────────────────────────────────────────────────────────────────────────────
def main():
    # Streamlit page configuration
    st.set_page_config(
        page_title=PAGE_TITLE,
        layout="wide",
        initial_sidebar_state="expanded",
        page_icon="🧬"
    )

    # Inject enhanced CSS with FIXED sidebar styling
    inject_custom_css()

    # Load user settings
    settings = load_settings()

    # Enhanced sidebar navigation
    selected = display_enhanced_sidebar()

    # FIXED: Enhanced page header - removed purple banners from module pages
    render_enhanced_header(selected)

    # Load and run the selected module
    run_selected_module(selected)

    # Enhanced footer
    render_enhanced_footer()

if __name__ == "__main__":
    main()