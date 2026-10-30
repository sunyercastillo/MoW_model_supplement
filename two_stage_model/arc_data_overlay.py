#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Arc Data Overlay Module for Bali et al. (2012) Dashboard
=========================================================

Provides functionality to overlay real arc geochemical data on model plots.
Reads arc data from Excel file in the project folder (arc_data.xlsx).
Computes element ratios internally from raw ppm concentrations.

Author: Arc Data Extension
Date: January 2026
"""

import os
import contextlib
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from typing import Dict, List, Tuple, Optional, Any
import warnings

# ============================================================================
# DEFAULT STYLING
# ============================================================================

# Default color palette for arcs (Plotly categorical colors)
DEFAULT_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5'
]

# Default symbols
DEFAULT_SYMBOLS = [
    'circle', 'square', 'diamond', 'triangle-up', 'triangle-down',
    'cross', 'x', 'star', 'hexagon', 'pentagon'
]

# Required columns in ArcData sheet
REQUIRED_COLUMNS = ['SampleID', 'Arc']

# Numeric columns that should be coerced (raw concentrations and isotopes only)
NUMERIC_COLUMNS = [
    'Mo_ppm', 'Ce_ppm', 'W_ppm', 'U_ppm', 'Th_ppm', 'Nb_ppm', 'La_ppm',
    'd98Mo', 'd98Mo_2SE',
    'Latitude', 'Longitude', 'ShowLabel', 'Size', 'Opacity', 'OutlineWidth', 'Visible'
]

# Default arc data filename
ARC_DATA_FILENAME = 'arc_data.xlsx'


# ============================================================================
# DATA LOADING AND VALIDATION
# ============================================================================

def load_arc_data_from_file(file_path: str = None) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], List[str]]:
    """
    Load arc data from Excel file in the project folder.
    
    Parameters
    ----------
    file_path : str, optional
        Path to Excel file. If None, looks for 'arc_data.xlsx' in the app directory.
    
    Returns
    -------
    tuple
        (arc_df, style_df, warnings_list)
        arc_df: DataFrame with arc sample data (with computed ratios)
        style_df: DataFrame with styling info (or None if missing)
        warnings_list: List of warning/info messages
    """
    warnings_list = []
    arc_df = None
    style_df = None
    
    # Determine file path
    if file_path is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_dir, ARC_DATA_FILENAME)
    
    # Check if file exists
    if not os.path.exists(file_path):
        # Not an error - just no arc data file present
        return None, None, []
    
    try:
        # Read Excel file
        excel_file = pd.ExcelFile(file_path)
        sheet_names = excel_file.sheet_names
        
        # Load ArcData sheet (required)
        if 'ArcData' in sheet_names:
            arc_df = pd.read_excel(excel_file, sheet_name='ArcData')
        elif len(sheet_names) > 0:
            # Try first sheet as fallback
            arc_df = pd.read_excel(excel_file, sheet_name=0)
            warnings_list.append(f"'ArcData' sheet not found. Using first sheet: '{sheet_names[0]}'")
        else:
            warnings_list.append("No sheets found in Excel file.")
            return None, None, warnings_list
        
        # Load Style sheet (optional)
        if 'Style' in sheet_names:
            style_df = pd.read_excel(excel_file, sheet_name='Style')
        else:
            warnings_list.append("'Style' sheet not found. Using default styling.")
        
        # Validate arc_df
        arc_df, validation_warnings = validate_arc_data(arc_df)
        warnings_list.extend(validation_warnings)
        
        # Validate style_df
        if style_df is not None:
            style_df, style_warnings = validate_style_data(style_df)
            warnings_list.extend(style_warnings)
        
    except Exception as e:
        warnings_list.append(f"Error loading arc data file: {str(e)}")
        return None, None, warnings_list
    
    return arc_df, style_df, warnings_list


def validate_arc_data(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    """
    Validate and clean arc data DataFrame.
    
    Returns
    -------
    tuple
        (cleaned_df, warnings_list)
    """
    warnings_list = []
    
    if df is None or df.empty:
        warnings_list.append("Arc data is empty.")
        return df, warnings_list
    
    # Check required columns
    missing_required = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_required:
        warnings_list.append(f"Missing required columns: {missing_required}")
    
    # Coerce numeric columns
    for col in NUMERIC_COLUMNS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Fill missing optional string columns
    for col in ['Group', 'Subgroup', 'Label', 'Volcano', 'Reference', 'Notes']:
        if col not in df.columns:
            df[col] = ''
        else:
            df[col] = df[col].fillna('').astype(str)
    
    # Set ShowLabel default
    if 'ShowLabel' not in df.columns:
        df['ShowLabel'] = 0
    else:
        df['ShowLabel'] = df['ShowLabel'].fillna(0).astype(int)
    
    # Compute ratios if missing
    df = compute_ratios(df)
    
    # Count valid samples
    n_samples = len(df)
    n_with_arc = df['Arc'].notna().sum() if 'Arc' in df.columns else 0
    
    if n_samples > 0:
        warnings_list.append(f"Loaded {n_samples} samples ({n_with_arc} with Arc designation)")
    
    return df, warnings_list


def validate_style_data(df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
    """
    Validate and clean style data DataFrame.
    """
    warnings_list = []
    
    if df is None or df.empty:
        return df, warnings_list
    
    # Check required columns
    required_style_cols = ['KeyType', 'Key']
    missing = [col for col in required_style_cols if col not in df.columns]
    if missing:
        warnings_list.append(f"Style sheet missing columns: {missing}")
        return None, warnings_list
    
    # Set defaults
    if 'Color' not in df.columns:
        df['Color'] = '#1f77b4'
    if 'Symbol' not in df.columns:
        df['Symbol'] = 'circle'
    if 'Size' not in df.columns:
        df['Size'] = 8
    if 'Opacity' not in df.columns:
        df['Opacity'] = 0.9
    if 'OutlineColor' not in df.columns:
        df['OutlineColor'] = '#000000'
    if 'OutlineWidth' not in df.columns:
        df['OutlineWidth'] = 1
    if 'Visible' not in df.columns:
        df['Visible'] = 1
    
    # Coerce numeric columns
    df['Size'] = pd.to_numeric(df['Size'], errors='coerce').fillna(8).astype(int)
    df['Opacity'] = pd.to_numeric(df['Opacity'], errors='coerce').fillna(0.9)
    df['OutlineWidth'] = pd.to_numeric(df['OutlineWidth'], errors='coerce').fillna(1).astype(int)
    df['Visible'] = pd.to_numeric(df['Visible'], errors='coerce').fillna(1).astype(int)
    
    return df, warnings_list


def compute_ratios(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute element ratios from raw ppm concentrations.
    
    Ratios computed:
    - Mo_Ce = Mo_ppm / Ce_ppm
    - U_Th = U_ppm / Th_ppm
    - W_Th = W_ppm / Th_ppm
    - Mo_W = Mo_ppm / W_ppm
    - Nb_La = Nb_ppm / La_ppm
    
    Parameters
    ----------
    df : pd.DataFrame
        Arc data with raw concentrations
    
    Returns
    -------
    pd.DataFrame
        DataFrame with computed ratio columns added
    """
    ratio_formulas = {
        'Mo_Ce': ('Mo_ppm', 'Ce_ppm'),
        'Ce_Mo': ('Ce_ppm', 'Mo_ppm'),
        'U_Th': ('U_ppm', 'Th_ppm'),
        'W_Th': ('W_ppm', 'Th_ppm'),
        'Mo_W': ('Mo_ppm', 'W_ppm'),
        'Nb_La': ('Nb_ppm', 'La_ppm'),
        'W_Ce': ('W_ppm', 'Ce_ppm'),
    }
    
    for ratio_col, (num_col, den_col) in ratio_formulas.items():
        # Always compute ratios from concentrations
        df[ratio_col] = np.nan
        
        if num_col in df.columns and den_col in df.columns:
            valid = df[num_col].notna() & df[den_col].notna() & (df[den_col] != 0)
            df.loc[valid, ratio_col] = df.loc[valid, num_col] / df.loc[valid, den_col]
    
    return df


# ============================================================================
# STYLING FUNCTIONS
# ============================================================================

def get_style_for_sample(sample: pd.Series, style_df: Optional[pd.DataFrame], 
                         color_by: str = 'Arc', symbol_by: str = 'Arc',
                         default_styles: Dict = None) -> Dict:
    """
    Get styling for a single sample based on style_df or defaults.
    
    Parameters
    ----------
    sample : pd.Series
        Single sample row
    style_df : pd.DataFrame or None
        Style definitions
    color_by : str
        'Arc' or 'Group'
    symbol_by : str
        'Arc', 'Group', or 'None'
    default_styles : dict
        Pre-computed default styles for arcs/groups
    
    Returns
    -------
    dict
        Style dictionary with color, symbol, size, etc.
    """
    style = {
        'color': '#1f77b4',
        'symbol': 'circle',
        'size': 8,
        'opacity': 0.9,
        'outline_color': '#000000',
        'outline_width': 1,
        'visible': True
    }
    
    # Get keys for lookup
    arc_key = sample.get('Arc', '')
    group_key = sample.get('Group', '')
    
    # Try to find style in style_df
    if style_df is not None and not style_df.empty:
        # First try exact match for color_by key
        color_key = arc_key if color_by == 'Arc' else group_key
        symbol_key = arc_key if symbol_by == 'Arc' else (group_key if symbol_by == 'Group' else None)
        
        # Look for color match
        for key_type in [color_by, 'Arc', 'Group']:
            key_val = arc_key if key_type == 'Arc' else group_key
            match = style_df[(style_df['KeyType'] == key_type) & (style_df['Key'] == key_val)]
            if not match.empty:
                row = match.iloc[0]
                style['color'] = row.get('Color', style['color'])
                style['size'] = int(row.get('Size', style['size']))
                style['opacity'] = float(row.get('Opacity', style['opacity']))
                style['outline_color'] = row.get('OutlineColor', style['outline_color'])
                style['outline_width'] = int(row.get('OutlineWidth', style['outline_width']))
                style['visible'] = bool(row.get('Visible', 1))
                break
        
        # Look for symbol match
        if symbol_key:
            for key_type in [symbol_by, 'Arc', 'Group']:
                key_val = arc_key if key_type == 'Arc' else group_key
                match = style_df[(style_df['KeyType'] == key_type) & (style_df['Key'] == key_val)]
                if not match.empty:
                    style['symbol'] = match.iloc[0].get('Symbol', style['symbol'])
                    break
    
    # Use default styles if no match found
    if default_styles:
        if color_by == 'Arc' and arc_key in default_styles:
            style['color'] = default_styles[arc_key].get('color', style['color'])
        elif color_by == 'Group' and group_key in default_styles:
            style['color'] = default_styles[group_key].get('color', style['color'])
        
        if symbol_by == 'Arc' and arc_key in default_styles:
            style['symbol'] = default_styles[arc_key].get('symbol', style['symbol'])
        elif symbol_by == 'Group' and group_key in default_styles:
            style['symbol'] = default_styles[group_key].get('symbol', style['symbol'])
    
    return style


def build_default_styles(arc_df: pd.DataFrame, style_df: Optional[pd.DataFrame]) -> Dict:
    """
    Build default style mapping for all unique arcs and groups.
    """
    default_styles = {}
    
    # Get unique arcs and groups
    arcs = arc_df['Arc'].dropna().unique().tolist() if 'Arc' in arc_df.columns else []
    groups = arc_df['Group'].dropna().unique().tolist() if 'Group' in arc_df.columns else []
    groups = [g for g in groups if g and str(g).strip()]
    
    # Assign colors and symbols
    for i, arc in enumerate(arcs):
        color_idx = i % len(DEFAULT_COLORS)
        symbol_idx = i % len(DEFAULT_SYMBOLS)
        default_styles[arc] = {
            'color': DEFAULT_COLORS[color_idx],
            'symbol': DEFAULT_SYMBOLS[symbol_idx]
        }
    
    for i, group in enumerate(groups):
        if group not in default_styles:
            color_idx = (i + len(arcs)) % len(DEFAULT_COLORS)
            symbol_idx = (i + len(arcs)) % len(DEFAULT_SYMBOLS)
            default_styles[group] = {
                'color': DEFAULT_COLORS[color_idx],
                'symbol': DEFAULT_SYMBOLS[symbol_idx]
            }
    
    # Override with style_df if provided
    if style_df is not None:
        for _, row in style_df.iterrows():
            key = row.get('Key', '')
            if key and key in default_styles:
                if pd.notna(row.get('Color')):
                    default_styles[key]['color'] = row['Color']
                if pd.notna(row.get('Symbol')):
                    default_styles[key]['symbol'] = row['Symbol']
            elif key:
                default_styles[key] = {
                    'color': row.get('Color', DEFAULT_COLORS[0]),
                    'symbol': row.get('Symbol', 'circle')
                }
    
    return default_styles


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def add_arc_overlay_to_figure(
    fig: go.Figure,
    arc_df: pd.DataFrame,
    style_df: Optional[pd.DataFrame],
    x_col: str,
    y_col: str,
    x_err_col: Optional[str] = None,
    y_err_col: Optional[str] = None,
    color_by: str = 'Arc',
    symbol_by: str = 'Arc',
    show_error_bars: bool = True,
    show_labels: bool = False,
    legend_by: str = 'Arc',
    filters: Dict = None
) -> go.Figure:
    """
    Add arc data overlay to an existing Plotly figure.
    
    Parameters
    ----------
    fig : go.Figure
        Existing Plotly figure
    arc_df : pd.DataFrame
        Arc sample data
    style_df : pd.DataFrame or None
        Style definitions
    x_col : str
        Column name for x-axis values
    y_col : str
        Column name for y-axis values
    x_err_col : str, optional
        Column name for x error bars
    y_err_col : str, optional
        Column name for y error bars
    color_by : str
        'Arc' or 'Group'
    symbol_by : str
        'Arc', 'Group', or 'None'
    show_error_bars : bool
        Whether to show error bars
    show_labels : bool
        Whether to show sample labels
    legend_by : str
        'Arc' or 'Group' for legend grouping
    filters : dict, optional
        Filter criteria {column: [values]}
    
    Returns
    -------
    go.Figure
        Figure with arc data overlay added
    """
    if arc_df is None or arc_df.empty:
        return fig

    # Required ratio columns must exist (compute_ratios normally creates them, but a
    # malformed arc_data.xlsx may be missing the underlying concentrations) — skip
    # the overlay gracefully rather than raising KeyError downstream.
    if x_col not in arc_df.columns or y_col not in arc_df.columns:
        return fig

    # Apply filters
    df = arc_df.copy()
    if filters:
        for col, values in filters.items():
            if col in df.columns and values:
                df = df[df[col].isin(values)]
    
    # Filter to samples with valid x and y values
    df = df[df[x_col].notna() & df[y_col].notna()].copy()
    
    if df.empty:
        return fig
    
    # Build default styles
    default_styles = build_default_styles(arc_df, style_df)
    
    # Group by legend key
    legend_col = 'Arc' if legend_by == 'Arc' else 'Group'
    if legend_col not in df.columns:
        legend_col = 'Arc'
    
    # Get unique legend keys
    legend_keys = df[legend_col].dropna().unique()
    
    # Track which legend entries have been added
    legend_added = set()
    
    for legend_key in legend_keys:
        if not legend_key or pd.isna(legend_key):
            continue
        
        group_df = df[df[legend_col] == legend_key]
        
        if group_df.empty:
            continue
        
        # Get style for this group (use first sample's style)
        sample_style = get_style_for_sample(
            group_df.iloc[0], style_df, color_by, symbol_by, default_styles
        )
        
        if not sample_style['visible']:
            continue
        
        x_vals = group_df[x_col].values
        y_vals = group_df[y_col].values
        
        # Error bars
        error_x = None
        error_y = None
        if show_error_bars:
            if x_err_col and x_err_col in group_df.columns:
                x_err_vals = group_df[x_err_col].values
                if np.any(np.isfinite(x_err_vals)):
                    error_x = dict(type='data', array=x_err_vals, visible=True)
            if y_err_col and y_err_col in group_df.columns:
                y_err_vals = group_df[y_err_col].values
                if np.any(np.isfinite(y_err_vals)):
                    error_y = dict(type='data', array=y_err_vals, visible=True, thickness=0.8)
        
        # Hover text
        hover_text = []
        for _, row in group_df.iterrows():
            text = f"<b>{row.get('SampleID', 'N/A')}</b><br>"
            text += f"Arc: {row.get('Arc', 'N/A')}<br>"
            if row.get('Group'):
                text += f"Group: {row['Group']}<br>"
            if row.get('Volcano'):
                text += f"Volcano: {row['Volcano']}<br>"
            text += f"{x_col}: {row[x_col]:.4g}<br>"
            text += f"{y_col}: {row[y_col]:.4g}"
            hover_text.append(text)
        
        # Add trace
        show_legend = legend_key not in legend_added
        legend_added.add(legend_key)
        
        fig.add_trace(go.Scatter(
            x=x_vals,
            y=y_vals,
            mode='markers',
            name=str(legend_key),
            marker=dict(
                color=sample_style['color'],
                symbol=sample_style['symbol'],
                size=sample_style['size'],
                opacity=sample_style['opacity'],
                line=dict(
                    color=sample_style['outline_color'],
                    width=sample_style['outline_width']
                )
            ),
            error_x=error_x,
            error_y=error_y,
            hovertemplate='%{hovertext}<extra></extra>',
            hovertext=hover_text,
            showlegend=show_legend,
            legendgroup=str(legend_key)
        ))
        
        # Add labels if requested
        if show_labels and 'ShowLabel' in group_df.columns:
            label_df = group_df[group_df['ShowLabel'] == 1]
            if not label_df.empty:
                _label_src = label_df['Label'] if 'Label' in label_df.columns else label_df.get('SampleID')
                labels = _label_src.fillna(label_df.get('SampleID', '')).values
                fig.add_trace(go.Scatter(
                    x=label_df[x_col].values,
                    y=label_df[y_col].values,
                    mode='text',
                    text=labels,
                    textposition='top center',
                    textfont=dict(size=8, color=sample_style['color']),
                    showlegend=False,
                    hoverinfo='skip'
                ))
    
    return fig


def get_plot_summary(arc_df: pd.DataFrame) -> Dict:
    """
    Get summary of which plots can show arc data.
    """
    if arc_df is None or arc_df.empty:
        return {}
    
    plots = {
        'Mo/Ce vs U/Th': ('Mo_Ce', 'U_Th'),
        'W/Th vs U/Th': ('W_Th', 'U_Th'),
        'Mo/W vs U/Th': ('Mo_W', 'U_Th'),
        'Mo/W vs Nb/La': ('Mo_W', 'Nb_La'),
        'δ98Mo vs Mo/Ce': ('Mo_Ce', 'd98Mo')
    }
    
    summary = {}
    for plot_name, (x_col, y_col) in plots.items():
        valid_count = arc_df[x_col].notna() & arc_df[y_col].notna()
        summary[plot_name] = {
            'x_col': x_col,
            'y_col': y_col,
            'valid_samples': valid_count.sum(),
            'total_samples': len(arc_df)
        }
    
    return summary


# ============================================================================
# STREAMLIT UI COMPONENTS
# ============================================================================

def render_arc_data_upload_section(expanded: bool = False, key_prefix: str = "",
                                   use_expander: bool = True) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], Dict]:
    """
    Render the arc data overlay section in Streamlit.
    Automatically loads arc data from 'arc_data.xlsx' in the project folder.

    Parameters
    ----------
    expanded : bool
        Whether the expander should be open by default (default False)
    key_prefix : str
        Prefix applied to every widget key and to the settings session-state key, so the
        panel can be rendered independently in more than one place (e.g. Part 1 and Part 2).
    use_expander : bool
        Wrap the controls in a "Lava data overlay" expander (True) or render them inline
        with a markdown header (False — use inside another expander to avoid nesting errors).

    Returns
    -------
    tuple
        (arc_df, style_df, settings_dict)
    """
    _sk = f"{key_prefix}arc_settings"  # per-panel settings key (data itself is shared)

    # Initialize session state for arc data (arc_df / style_df / loaded flag are shared)
    if 'arc_df' not in st.session_state:
        st.session_state.arc_df = None
    if 'style_df' not in st.session_state:
        st.session_state.style_df = None
    if _sk not in st.session_state:
        st.session_state[_sk] = {
            'color_by': 'Group',
            'symbol_by': 'Group',
            'show_error_bars': True,
            'show_labels': False,
            'legend_by': 'Group',
            'filters': {}
        }
    if 'arc_data_loaded' not in st.session_state:
        st.session_state.arc_data_loaded = False

    arc_df = None
    style_df = None
    settings = st.session_state[_sk].copy()

    # Load arc data from file (only once per session or on reload)
    if not st.session_state.arc_data_loaded:
        arc_df, style_df, warnings = load_arc_data_from_file()

        if arc_df is not None and not arc_df.empty:
            st.session_state.arc_df = arc_df
            st.session_state.style_df = style_df
            st.session_state.arc_data_loaded = True
        else:
            st.session_state.arc_data_loaded = True  # Mark as checked even if no file

    # Get data from session state
    arc_df = st.session_state.arc_df
    style_df = st.session_state.style_df

    # Container: an expander, or an inline block (to avoid nested-expander errors)
    _box = st.expander("Lava data overlay", expanded=expanded) if use_expander else contextlib.nullcontext()
    with _box:
        if not use_expander:
            st.markdown("**Lava data overlay**")
        # Reload button
        col_reload, col_info = st.columns([1, 3])
        with col_reload:
            if st.button("🔄 Reload", key=f"{key_prefix}reload_arc_data"):
                st.session_state.arc_data_loaded = False
                st.rerun()
        with col_info:
            st.caption("Reads from `arc_data.xlsx`")

        if arc_df is None or arc_df.empty:
            st.info(
                "📁 **No arc data file found.**\n\n"
                "Create `arc_data.xlsx` with ArcData sheet."
            )
            return None, None, settings

        # Show data summary
        st.success(f"✅ {len(arc_df)} arc samples loaded")

        # Settings in compact layout
        col1, col2 = st.columns(2)

        with col1:
            settings['color_by'] = st.radio(
                "Color by:",
                options=['Arc', 'Group'],
                index=0 if settings['color_by'] == 'Arc' else 1,
                key=f'{key_prefix}arc_color_by',
                horizontal=True,
                label_visibility='collapsed'
            )
            st.caption("Color by")
            settings['show_error_bars'] = st.checkbox(
                "Show error bars",
                value=settings['show_error_bars'],
                key=f'{key_prefix}arc_show_errors'
            )

        with col2:
            symbol_options = ['Arc', 'Group', 'None']
            settings['symbol_by'] = st.radio(
                "Symbol by:",
                options=symbol_options,
                index=symbol_options.index(settings['symbol_by']) if settings['symbol_by'] in symbol_options else 1,
                key=f'{key_prefix}arc_symbol_by',
                horizontal=True,
                label_visibility='collapsed'
            )
            st.caption("Symbol by")
            settings['show_labels'] = st.checkbox(
                "Show labels",
                value=settings['show_labels'],
                key=f'{key_prefix}arc_show_labels'
            )

        settings['legend_by'] = st.radio(
            "Legend by:",
            options=['Arc', 'Group'],
            index=0 if settings['legend_by'] == 'Arc' else 1,
            key=f'{key_prefix}arc_legend_by',
            horizontal=True,
            label_visibility='collapsed'
        )
        st.caption("Legend by")

        # Filters in nested expander
        with st.expander("🔍 Filter Data", expanded=False):
            filter_col1, filter_col2 = st.columns(2)

            with filter_col1:
                arcs = arc_df['Arc'].dropna().unique().tolist()
                selected_arcs = st.multiselect(
                    "Filter by Arc:",
                    options=arcs,
                    default=arcs,
                    key=f'{key_prefix}arc_filter_arc'
                )
                if selected_arcs and len(selected_arcs) < len(arcs):
                    settings['filters']['Arc'] = selected_arcs
                elif 'Arc' in settings['filters']:
                    del settings['filters']['Arc']

            with filter_col2:
                groups = [g for g in arc_df['Group'].dropna().unique().tolist() if g]
                if groups:
                    selected_groups = st.multiselect(
                        "Filter by Group:",
                        options=groups,
                        default=groups,
                        key=f'{key_prefix}arc_filter_group'
                    )
                    if selected_groups and len(selected_groups) < len(groups):
                        settings['filters']['Group'] = selected_groups
                    elif 'Group' in settings['filters']:
                        del settings['filters']['Group']

    # Update session state
    st.session_state[_sk] = settings

    return arc_df, style_df, settings
