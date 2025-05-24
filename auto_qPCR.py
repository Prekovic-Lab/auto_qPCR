import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO

st.set_page_config(layout='wide', page_title='Prekovic Lab qPCR Analysis')
st.title('🧬 Prekovic Lab qPCR Data Analysis')

uploaded_file = st.file_uploader("📁 Upload qPCR Excel file (.xlsx)", type=['xlsx'])

if uploaded_file:
    xl = pd.ExcelFile(uploaded_file)

    # Properly loading sheets with correct header rows
    results = xl.parse('Results', skiprows=21)
    amplification_data = xl.parse('Amplification Data', skiprows=21)
    melt_curve_raw = xl.parse('Melt Curve Raw Data', skiprows=21)

    # Dynamic housekeeping gene selection from actual data
    available_genes = results['Target Name'].unique().tolist()
    
    st.sidebar.header("Housekeeping Gene Settings")
    housekeeping_genes = st.sidebar.multiselect(
        "Select Housekeeping Genes:",
        options=available_genes,
        default=[gene for gene in available_genes if gene.lower() in ['ubc', 'actb', 'gapdh', 'b-actin']]
    )

    norm_method = st.sidebar.radio("Normalization Method:", ['Individual', 'Geometric Mean'], index=1)

    # Genes of interest selection dynamically
    genes_of_interest = st.sidebar.multiselect(
        "Genes of Interest:", 
        options=[gene for gene in available_genes if gene not in housekeeping_genes],
        default=[gene for gene in available_genes if gene not in housekeeping_genes][:3]
    )

    # Ensure housekeeping genes were selected
    if not housekeeping_genes:
        st.error("Please select at least one housekeeping gene.")
        st.stop()

    # Convert CT to numeric explicitly
    results['CT'] = pd.to_numeric(results['CT'], errors='coerce')

    # Calculate housekeeping reference
    housekeeping_cts = results[results['Target Name'].isin(housekeeping_genes)].pivot_table(
        index='Sample Name', columns='Target Name', values='CT', aggfunc='mean'
    )

    if norm_method == 'Individual':
        hk_reference = housekeeping_cts.mean(axis=1)
    else:
        hk_reference = housekeeping_cts.apply(lambda x: np.exp(np.mean(np.log(x.dropna()))), axis=1)

    # Normalization (ΔCt calculation)
    results['ΔCt'] = results.apply(lambda row: row['CT'] - hk_reference.loc[row['Sample Name']], axis=1)

    # Melt Curve Visualization (individual well)
    st.header('🔥 Melt Curve Analysis')
    well_selected = st.text_input("Enter Well Position (e.g., A1, B12):", value='A1').upper()
    melt_curve_selected = melt_curve_raw[melt_curve_raw['Well Position'] == well_selected]

    if not melt_curve_selected.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.lineplot(x='Temperature', y='Derivative', data=melt_curve_selected, ax=ax)
        ax.set_title(f"Melt Curve for Well {well_selected}")
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Derivative")
        st.pyplot(fig)
    else:
        st.warning(f"No melt curve data found for well {well_selected}.")

    # Amplification Curves
    st.header('📈 Amplification Curves')
    amp_gene = st.selectbox("Select Gene for Amplification Curve:", available_genes)
    amp_curve_data = amplification_data[amplification_data['Target Name'] == amp_gene]

    fig2, ax2 = plt.subplots(figsize=(7, 4))
    for well in amp_curve_data['Well Position'].unique():
        well_data = amp_curve_data[amp_curve_data['Well Position'] == well]
        ax2.plot(well_data['Cycle'], well_data['ΔRn'], alpha=0.3)

    ax2.set_title(f"Amplification Curves for {amp_gene}")
    ax2.set_xlabel("Cycle")
    ax2.set_ylabel("ΔRn")
    st.pyplot(fig2)

    # Barplot of Normalized Gene Expression
    st.header("📊 Normalized Gene Expression")
    plot_data = results[results['Target Name'].isin(genes_of_interest)].dropna(subset=['ΔCt'])

    fig3, ax3 = plt.subplots(figsize=(8, 5))
    sns.barplot(
        x='Target Name', y='ΔCt', data=plot_data, ci='sd', capsize=.2, palette='viridis', ax=ax3
    )
    sns.stripplot(
        x='Target Name', y='ΔCt', data=plot_data, color='black', alpha=0.7, ax=ax3
    )
    ax3.set_xlabel("Gene")
    ax3.set_ylabel("Normalized Expression (ΔCt)")
    ax3.set_title("Normalized Expression of Selected Genes")
    st.pyplot(fig3)

    # Export barplot as PDF
    buffer = BytesIO()
    fig3.savefig(buffer, format='pdf')
    buffer.seek(0)
    st.download_button(
        label="⬇️ Download Barplot (PDF)",
        data=buffer,
        file_name="qPCR_Gene_Expression.pdf",
        mime="application/pdf"
    )

    # Export normalized results CSV
    normalized_csv = results.to_csv(index=False).encode('utf-8')
    st.download_button(
        "⬇️ Download Normalized Data (CSV)",
        data=normalized_csv,
        file_name="normalized_qPCR_results.csv",
        mime="text/csv"
    )

else:
    st.info('Please upload a valid qPCR Excel (.xlsx) file to proceed.')

# Footer
st.markdown("---")
st.markdown("🔗 [Prekovic Lab](https://www.prekovic-lab.org) | ✉️ [Contact](mailto:s.prekovic@umcutrecht.nl)")
