import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from io import BytesIO

st.set_page_config(layout='wide', page_title='Prekovic Lab qPCR Analysis')
st.title('Prekovic Lab qPCR Data Analysis')

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
        hk_means = hk_data.groupby('Sample Name')['CT'].agg(
            lambda x: np.exp(np.mean(np.log(x.dropna()))) if norm_method == 'Geometric Mean' else np.mean(x)
        )

        results = results.join(hk_means.rename('HK_mean'), on='Sample Name')
        results['ΔCt'] = results['CT'] - results['HK_mean']
        results['Expression (2^-ΔCt)'] = 2 ** (-results['ΔCt'])

        # Extract condition and replicate
        results[['Condition', 'Replicate']] = results['Sample Name'].str.extract(r'(.*)\s+(\d+)$')

        plot_data = results[results['Target Name'].isin(genes_of_interest)].dropna(subset=['Expression (2^-ΔCt)'])

        st.header("Normalized Gene Expression")

        for gene in genes_of_interest:
            gene_data = plot_data[plot_data['Target Name'] == gene]
            summary = gene_data.groupby('Condition')['Expression (2^-ΔCt)'].agg(['mean', 'std']).reset_index()

            y_max = summary['mean'].max() + summary['std'].max() * 1.2

            fig = go.Figure()

            # Add barplot with error bars
            fig.add_trace(go.Bar(
                x=summary['Condition'],
                y=summary['mean'],
                error_y=dict(type='data', array=summary['std']),
                marker_color=px.colors.qualitative.Pastel,
                name='Mean ± SD'
            ))

            # Add individual replicate points
            for idx, row in gene_data.iterrows():
                fig.add_trace(go.Scatter(
                    x=[row['Condition']],
                    y=[row['Expression (2^-ΔCt)']],
                    mode='markers',
                    marker=dict(size=8, color='black', opacity=0.6),
                    hovertemplate=(
                        f"Condition: {row['Condition']}<br>"
                        f"Replicate: {row['Replicate']}<br>"
                        f"Well: {row['Well Position']}<br>"
                        f"Expression: {row['Expression (2^-ΔCt)']:.2f}<br>"
                        f"ΔCt: {row['ΔCt']:.2f}"
                    ),
                    showlegend=False
                ))

            # Force the y-axis to start exactly at 0
            fig.update_yaxes(range=[0, summary['mean'].max() + summary['std'].max()*1.2])


            fig.update_layout(
                title=f'Normalized Expression of {gene}',
                xaxis_title='Condition',
                yaxis_title='Expression (2^-ΔCt)',
                template='simple_white',
                width=800,
                height=500,
                    yaxis=dict(
        range=[0, y_max],   # force baseline at zero
        autorange=False,    # disable autorange to respect range
        zeroline=True,      # draw zero baseline
        zerolinewidth=2     # thickness of zero line
    )
            )

            st.plotly_chart(fig, use_container_width=True)

        # Export normalized results CSV
        normalized_csv = results[['Sample Name', 'Condition', 'Replicate', 'Well Position', 'Target Name', 'CT', 'HK_mean', 'ΔCt', 'Expression (2^-ΔCt)']].to_csv(index=False).encode('utf-8')
        st.download_button(
            "⬇️ Download Normalized Data (CSV)",
            data=normalized_csv,
            file_name="normalized_qPCR_results.csv",
            mime="text/csv"
        )

    elif analysis_type == 'Melt Curve Analysis':
        st.header('Melt Curve Visualization')
        well_selected = st.text_input("Well Position (e.g., A1, B12):", 'A1').upper()
        melt_curve_selected = melt_curve_raw[melt_curve_raw['Well Position'] == well_selected]

        if not melt_curve_selected.empty:
            fig = px.line(melt_curve_selected, x='Temperature', y='Derivative',
                          title=f"Melt Curve for Well {well_selected}",
                          template='simple_white')
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning(f"No data found for well {well_selected}.")

else:
    st.info('Upload a valid qPCR Excel (.xlsx) file.')

st.markdown("---")
st.markdown("[Prekovic Lab](https://www.prekovic-lab.org) | [Contact](mailto:s.prekovic@umcutrecht.nl)")
