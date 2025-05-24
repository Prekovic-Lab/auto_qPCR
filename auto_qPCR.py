import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO

st.set_page_config(layout='wide', page_title='Prekovic Lab qPCR Analysis')

st.title('🧬 Prekovic Lab qPCR Data Analysis')

# File upload
uploaded_file = st.file_uploader("📁 Upload qPCR Excel file (.xlsx)", type=['xlsx'])

if uploaded_file:
    xl = pd.ExcelFile(uploaded_file)

    # Load specific sheets
    sample_setup = xl.parse('Sample Setup', skiprows=6)
    amplification_data = xl.parse('Amplification Data', skiprows=7)
    results = xl.parse('Results', skiprows=20)
    melt_curve_raw = xl.parse('Melt Curve Raw Data', skiprows=7)
    melt_curve_result = xl.parse('Melt Curve Result', skiprows=7)

    # Housekeeping gene selection
    st.sidebar.header("Housekeeping Gene Settings")
    housekeeping_genes = st.sidebar.multiselect(
        "Select Housekeeping Genes:",
        options=['UBC', 'ACTB', 'GAPDH'],
        default=['UBC', 'ACTB']
    )
    norm_method = st.sidebar.radio(
        "Normalization Method:",
        options=['Individual', 'Geometric Mean'],
        index=1
    )

    # Genes of interest selection
    gene_options = results['Target Name'].unique().tolist()
    genes_of_interest = st.sidebar.multiselect("Genes of Interest:", gene_options, default=gene_options[:2])

    # Normalization calculations
    def geometric_mean(ct_values):
        return np.exp(np.mean(np.log(ct_values), axis=0))

    normalized_results = results.copy()
    housekeeping_cts = results[results['Target Name'].isin(housekeeping_genes)].pivot_table(
        index='Sample Name', columns='Target Name', values='CT')

    if norm_method == 'Individual':
        hk_reference = housekeeping_cts.mean(axis=1)
    else:  # Geometric Mean
        hk_reference = housekeeping_cts.apply(geometric_mean, axis=1)

    normalized_results['ΔCt'] = normalized_results.apply(
        lambda row: row['CT'] - hk_reference.loc[row['Sample Name']], axis=1
    )

    # Melt curve visualization
    st.header('🔥 Melt Curve Analysis')
    well_selected = st.text_input("Enter Well Position (e.g., A1, B12):", value='A1').upper()
    melt_curve_selected = melt_curve_raw[melt_curve_raw['Well'] == well_selected]

    if not melt_curve_selected.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.lineplot(x='Temperature', y='Derivative', data=melt_curve_selected, ax=ax)
        ax.set_title(f"Melt Curve for Well {well_selected}")
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Derivative")
        st.pyplot(fig)
    else:
        st.warning(f"No melt curve data found for well {well_selected}.")

    # Amplification Curve Visualization
    st.header('📈 Amplification Curves')
    amp_gene = st.selectbox("Select Gene for Amplification Curve:", gene_options)

    amp_curve_data = amplification_data[amplification_data['Target Name'] == amp_gene]

    fig2, ax2 = plt.subplots(figsize=(7, 4))
    for well in amp_curve_data['Well'].unique():
        well_data = amp_curve_data[amp_curve_data['Well'] == well]
        ax2.plot(well_data['Cycle'], well_data['ΔRn'], alpha=0.3)

    ax2.set_title(f"Amplification Curves for {amp_gene}")
    ax2.set_xlabel("Cycle")
    ax2.set_ylabel("ΔRn")
    st.pyplot(fig2)

    # Barplot with replicates and standard deviation (Genes of interest)
    st.header("📊 Normalized Gene Expression")

    plot_data = normalized_results[normalized_results['Target Name'].isin(genes_of_interest)]

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

    # Export barplot as vectorized PDF
    buffer = BytesIO()
    fig3.savefig(buffer, format='pdf')
    buffer.seek(0)
    st.download_button(
        label="⬇️ Download Barplot (PDF)",
        data=buffer,
        file_name="qPCR_Gene_Expression.pdf",
        mime="application/pdf"
    )

    # Export normalized results
    csv_data = normalized_results.to_csv(index=False).encode('utf-8')
    st.download_button(
        "⬇️ Download Normalized Data (CSV)",
        data=csv_data,
        file_name="normalized_qPCR_results.csv",
        mime="text/csv"
    )

else:
    st.info('Please upload a valid qPCR Excel (.xlsx) file to proceed.')

# Footer
st.markdown("---")
st.markdown("🔗 [Prekovic Lab](https://www.prekovic-lab.org) | ✉️ [Contact](mailto:s.prekovic@umcutrecht.nl)")
