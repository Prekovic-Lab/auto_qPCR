import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

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
        hk_means = hk_data.groupby('Sample Name')['CT'].agg(
            lambda x: np.exp(np.mean(np.log(x.dropna()))) if norm_method == 'Geometric Mean' else np.mean(x)
        )

        results = results.join(hk_means.rename('HK_mean'), on='Sample Name')
        results['ΔCt'] = results['CT'] - results['HK_mean']
        results['Expression (2^-ΔCt)'] = 2 ** (-results['ΔCt'])

        results[['Condition', 'Replicate']] = results['Sample Name'].str.extract(r'(.*)\s+(\d+)$')

        plot_data = results[results['Target Name'].isin(genes_of_interest)].dropna(subset=['Expression (2^-ΔCt)'])

        st.header("📊 Normalized Gene Expression")

        conditions = plot_data['Condition'].unique()
        condition_colors = {}
        
        st.subheader("🎨 Choose Condition Colors:")
        cols = st.columns(len(conditions))
        default_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
        
        for idx, cond in enumerate(conditions):
            with cols[idx]:
                condition_colors[cond] = st.color_picker(f"Color for {cond}", default_colors[idx % len(default_colors)])

        for gene in genes_of_interest:
            gene_data = plot_data[plot_data['Target Name'] == gene]
            summary = gene_data.groupby('Condition')['Expression (2^-ΔCt)'].agg(['mean', 'std']).reset_index()

            fig = go.Figure()

            for cond in conditions:
                cond_summary = summary[summary['Condition'] == cond]
                cond_data = gene_data[gene_data['Condition'] == cond]

                fig.add_trace(go.Bar(
                    x=[cond],
                    y=cond_summary['mean'],
                    error_y=dict(type='data', array=cond_summary['std']),
                    width=0.7,
                    marker=dict(color=condition_colors[cond], line=dict(color='black', width=1.5)),
                    name=cond
                ))

                # Explicit individual data points
                fig.add_trace(go.Scatter(
                    x=np.repeat(cond, len(cond_data)),
                    y=cond_data['Expression (2^-ΔCt)'],
                    mode='markers',
                    marker=dict(size=7, color='black', opacity=0.9, line=dict(width=1)),
                    hovertemplate=(
                        f"Condition: {cond}<br>" +
                        "Replicate: %{customdata[0]}<br>" +
                        "Well: %{customdata[1]}<br>" +
                        "Expression: %{y:.2f}<br>" +
                        "ΔCt: %{customdata[2]:.2f}<extra></extra>"
                    ),
                    customdata=cond_data[['Replicate', 'Well Position', 'ΔCt']],
                    showlegend=False
                ))

            fig.update_layout(
                title=f'Normalized Expression of {gene}',
                xaxis_title='Condition',
                yaxis_title='Expression (2^-ΔCt)',
                template='plotly_white',
                plot_bgcolor='white',
                paper_bgcolor='white',
                width=500,
                height=450,
                font=dict(color='black'),
                margin=dict(l=40, r=40, t=40, b=40),
                xaxis=dict(
                    showline=True, linewidth=1, linecolor='black', mirror=True,
                    ticks='outside', showgrid=False
                ),
                yaxis=dict(
                    showline=True, linewidth=1, linecolor='black', mirror=True,
                    ticks='outside', showgrid=True, gridcolor='lightgray', rangemode='tozero'
                )
            )

            st.plotly_chart(fig, use_container_width=False)

            # PDF vectorized plot download
            try:
                pdf_bytes = pio.to_image(fig, format='pdf', engine='kaleido')
                st.download_button(
                    label=f"⬇️ Download {gene} Plot (PDF)",
                    data=pdf_bytes,
                    file_name=f"{gene}_Expression_Plot.pdf",
                    mime='application/pdf'
                )
            except Exception as e:
                st.error(f"PDF Error: {e}. Ensure 'kaleido' is installed.")

        normalized_csv = results[['Sample Name', 'Condition', 'Replicate', 'Well Position', 'Target Name', 'CT', 'HK_mean', 'ΔCt', 'Expression (2^-ΔCt)']].to_csv(index=False).encode('utf-8')
        st.download_button(
            "⬇️ Download Normalized Data (CSV)",
            data=normalized_csv,
            file_name="normalized_qPCR_results.csv",
            mime="text/csv"
        )

    elif analysis_type == 'Melt Curve Analysis':
        st.header('🔥 Melt Curve Visualization')
        well_selected = st.text_input("Well Position (e.g., A1, B12):", 'A1').upper()
        melt_curve_selected = melt_curve_raw[melt_curve_raw['Well Position'] == well_selected]

        if not melt_curve_selected.empty:
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=melt_curve_selected['Temperature'],
                y=melt_curve_selected['Derivative'],
                mode='lines',
                line=dict(color='royalblue', width=2)
            ))

            fig.update_layout(
                title=f"Melt Curve for Well {well_selected}",
                xaxis_title='Temperature (°C)',
                yaxis_title='Derivative',
                template='plotly_white',
                plot_bgcolor='white',
                paper_bgcolor='white',
                width=700,
                height=450,
                font=dict(color='black'),
                margin=dict(l=40, r=40, t=40, b=40),
                xaxis=dict(
                    showline=True, linewidth=1, linecolor='black', mirror=True,
                    ticks='outside', showgrid=False
                ),
                yaxis=dict(
                    showline=True, linewidth=1, linecolor='black', mirror=True,
                    ticks='outside', showgrid=True, gridcolor='lightgray'
                )
            )

            st.plotly_chart(fig, use_container_width=False)

            try:
                pdf_bytes_melt = pio.to_image(fig, format='pdf', engine='kaleido')
                st.download_button(
                    label=f"⬇️ Download Melt Curve {well_selected} (PDF)",
                    data=pdf_bytes_melt,
                    file_name=f"MeltCurve_{well_selected}.pdf",
                    mime='application/pdf'
                )
            except Exception as e:
                st.error(f"PDF Error: {e}. Ensure 'kaleido' is installed.")
        else:
            st.warning(f"No data found for well {well_selected}.")

else:
    st.info('Upload a valid qPCR Excel (.xlsx) file.')

st.markdown("---")
st.markdown("[Prekovic Lab](https://www.prekovic-lab.org) | [Contact](mailto:s.prekovic@umcutrecht.nl)")
