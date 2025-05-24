import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO

st.set_page_config(layout='wide', page_title='Prekovic Lab qPCR Analysis')
st.title('🧬 Prekovic Lab qPCR Data Analysis')

uploaded_file = st.file_uploader("📁 Upload qPCR Excel file (.xlsx)", type=['xlsx'])

analysis_type = st.sidebar.radio("Analysis Type:", ['Amplification Results', 'Melt Curve Analysis'])

if uploaded_file:
    xl = pd.ExcelFile(uploaded_file)
    results = xl.parse('Results', skiprows=21)
    melt_curve_raw = xl.parse('Melt Curve Raw Data', skiprows=21)

    results['CT'] = pd.to_numeric(results['CT'].replace('Undetermined', np.nan), errors='coerce')
    available_genes = results['Target Name'].dropna().unique().tolist()

    if analysis_type == 'Amplification Results':
        st.sidebar.header("Normalization Settings")
        housekeeping_genes = st.sidebar.multiselect(
            "Housekeeping Genes:", options=available_genes,
            default=[gene for gene in available_genes if gene.lower() in ['ubc', 'actb', 'b-actin', 'gapdh']]
        )

        norm_method = st.sidebar.radio("Normalization method:", ['Individual Mean', 'Geometric Mean'], index=1)

        genes_of_interest = st.sidebar.multiselect(
            "Genes of Interest:",
            options=[gene for gene in available_genes if gene not in housekeeping_genes],
            default=[gene for gene in available_genes if gene not in housekeeping_genes][:3]
        )

        if not housekeeping_genes:
            st.error("Select at least one housekeeping gene.")
            st.stop()

        hk_data = results[results['Target Name'].isin(housekeeping_genes)]
        hk_means = hk_data.groupby('Sample Name')['CT'].agg(lambda x: np.exp(np.mean(np.log(x.dropna()))) if norm_method=='Geometric Mean' else np.mean(x))

        results = results.join(hk_means.rename('HK_mean'), on='Sample Name')
        results['ΔCt'] = results['CT'] - results['HK_mean']

        plot_data = results[results['Target Name'].isin(genes_of_interest)].dropna(subset=['ΔCt'])

        st.header("📊 Normalized Gene Expression (ΔCt)")

        plt.figure(figsize=(8, 5))
        sns.barplot(x='Target Name', y='ΔCt', data=plot_data, ci='sd', capsize=.2, palette='viridis')
        sns.stripplot(x='Target Name', y='ΔCt', data=plot_data, color='black', alpha=0.6)
        plt.xlabel("Gene")
        plt.ylabel("Normalized Expression (ΔCt)")
        plt.title("Normalized Expression of Selected Genes")
        st.pyplot(plt.gcf())

        # PDF Export
        buffer = BytesIO()
        plt.savefig(buffer, format='pdf')
        buffer.seek(0)
        st.download_button("⬇️ Download Barplot PDF", buffer, "qPCR_Gene_Expression.pdf", "application/pdf")

        # CSV Export
        csv = results[['Sample Name', 'Target Name', 'CT', 'HK_mean', 'ΔCt']].to_csv(index=False).encode()
        st.download_button("⬇️ Download Normalized Data (CSV)", csv, "normalized_qPCR_results.csv", "text/csv")

    elif analysis_type == 'Melt Curve Analysis':
        st.header('🔥 Melt Curve Visualization')
        well_selected = st.text_input("Well Position (e.g., A1, B12):", 'A1').upper()
        melt_curve_selected = melt_curve_raw[melt_curve_raw['Well Position'] == well_selected]

        if not melt_curve_selected.empty:
            plt.figure(figsize=(6,4))
            sns.lineplot(x='Temperature', y='Derivative', data=melt_curve_selected)
            plt.title(f"Melt Curve for Well {well_selected}")
            plt.xlabel("Temperature (°C)")
            plt.ylabel("Derivative")
            st.pyplot(plt.gcf())
        else:
            st.warning(f"No data found for well {well_selected}.")

else:
    st.info('Upload a valid qPCR Excel (.xlsx) file.')

st.markdown("---")
st.markdown("[Prekovic Lab](https://www.prekovic-lab.org) | [Contact](mailto:s.prekovic@umcutrecht.nl)")
