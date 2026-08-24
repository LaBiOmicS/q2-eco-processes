import os
import json
import pandas as pd
import biom
try:
    import qiime2
except ImportError:
    qiime2 = None

def summarize_processes(
    output_dir: str,
    table: biom.Table,
    metadata: qiime2.Metadata = None
) -> None:
    """
    QIIME 2 Visualizer Action: Generates interactive Bootstrap 5 HTML summary dashboard (.qzv) for ecological assembly processes,
    featuring Plotly.js 2D quadrant scatterplots (BetaNTI vs RCbray), donut charts, and DataTables.
    """
    # Convert biom.Table to pandas DataFrame (samples x processes)
    df_proc = pd.DataFrame(
        table.matrix_data.toarray().T,
        index=table.ids(axis='sample'),
        columns=table.ids(axis='observation')
    )

    non_proc_cols = ["Normalized_Stochasticity_Ratio", "Mean_BetaNTI", "Mean_RCbray"]
    proc_cols = [c for c in df_proc.columns if c not in non_proc_cols]
    proc_totals = df_proc[proc_cols].mean(axis=0).sort_values(ascending=False)
    mean_nst = round(df_proc["Normalized_Stochasticity_Ratio"].mean(), 2) if "Normalized_Stochasticity_Ratio" in df_proc.columns else 50.0

    total_samples = len(df_proc)

    # Coordinates for Plotly 2D Scatterplot
    x_coords = df_proc["Mean_RCbray"].round(3).tolist() if "Mean_RCbray" in df_proc.columns else [0.0] * total_samples
    y_coords = df_proc["Mean_BetaNTI"].round(3).tolist() if "Mean_BetaNTI" in df_proc.columns else [0.0] * total_samples
    sample_labels = df_proc.index.tolist()

    # Export TSV table inside visualization bundle
    tsv_path = os.path.join(output_dir, "ecological_assembly_processes.tsv")
    df_proc.to_csv(tsv_path, sep="\t")

    # Generate HTML rows
    proc_rows_html = ""
    for proc, mean_val in proc_totals.items():
        val_fmt = round(mean_val, 2)
        proc_clean = proc.replace("_", " ")
        badge_cls = "bg-primary" if "Selection" in proc else ("bg-warning text-dark" if "Dispersal" in proc else "bg-secondary")
        proc_rows_html += f"""
        <tr>
            <td><strong>{proc_clean}</strong></td>
            <td><code>{proc}</code></td>
            <td>
                <div class="progress" style="height: 22px;">
                    <div class="progress-bar {badge_cls}" role="progressbar" style="width: {val_fmt}%;" aria-valuenow="{val_fmt}" aria-valuemin="0" aria-valuemax="100">{val_fmt}%</div>
                </div>
            </td>
        </tr>
        """

    # HTML Dashboard Template with Plotly.js
    index_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>q2-eco-processes: Ecological Assembly Dashboard</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    <style>
        body {{ background-color: #f8f9fa; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }}
        .header-box {{ background: linear-gradient(135deg, #2b5876 0%, #4e4376 100%); color: white; padding: 2rem; border-radius: 12px; margin-bottom: 2rem; }}
        .card-custom {{ border: none; border-radius: 12px; box-shadow: 0 4px 12px rgba(0,0,0,0.05); margin-bottom: 2rem; background: white; }}
        .alert-custom {{ border-left: 5px solid #17a2b8; background-color: #e9f7ef; border-radius: 8px; }}
    </style>
</head>
<body>
<div class="container-fluid py-4 px-4">
    <div class="header-box d-flex justify-content-between align-items-center">
        <div>
            <h1 class="display-5 fw-bold mb-1">🌲 q2-eco-processes Dashboard</h1>
            <p class="lead mb-0">Quantification of Ecological Assembly Processes (Stegen et al. & Ning et al.)</p>
        </div>
        <a href="ecological_assembly_processes.tsv" class="btn btn-light fw-bold" download>📥 Download TSV Table</a>
    </div>

    <div class="row text-center mb-4">
        <div class="col-md-4">
            <div class="card card-custom p-3">
                <h6 class="text-muted text-uppercase">Total Samples / Groups</h6>
                <h2 class="fw-bold text-success">{total_samples:,}</h2>
            </div>
        </div>
        <div class="col-md-4">
            <div class="card card-custom p-3">
                <h6 class="text-muted text-uppercase">Processes Quantified</h6>
                <h2 class="fw-bold text-primary">5 Stegen Processes</h2>
            </div>
        </div>
        <div class="col-md-4">
            <div class="card card-custom p-3">
                <h6 class="text-muted text-uppercase">Mean NST Stochasticity Ratio</h6>
                <h2 class="fw-bold text-info">{mean_nst}%</h2>
            </div>
        </div>
    </div>

    <div class="row mb-4">
        <div class="col-md-6">
            <div class="card card-custom p-3">
                <h5 class="fw-bold mb-3">🍩 Ecological Assembly Composition</h5>
                <div id="donutPlot" style="height: 350px;"></div>
            </div>
        </div>
        <div class="col-md-6">
            <div class="card card-custom p-3">
                <h5 class="fw-bold mb-3">📈 2D Boundary Scatterplot (&beta;NTI vs RC<sub>bray</sub>)</h5>
                <div id="scatterPlot" style="height: 350px;"></div>
            </div>
        </div>
    </div>

    <div class="card card-custom">
        <div class="card-header bg-white py-3">
            <h5 class="m-0 fw-bold">📊 Relative Influence of Ecological Processes</h5>
        </div>
        <div class="card-body p-0">
            <div class="table-responsive">
                <table class="table table-hover align-middle mb-0">
                    <thead class="table-light">
                        <tr>
                            <th>Ecological Process</th>
                            <th>Process Code</th>
                            <th style="width: 50%;">Mean Relative Influence (%)</th>
                        </tr>
                    </thead>
                    <tbody>
                        {proc_rows_html}
                    </tbody>
                </table>
            </div>
        </div>
    </div>

    <div class="alert alert-custom p-3 mb-4">
        <h6 class="fw-bold text-dark m-0">💡 Methodological Note on Undominated Drift:</h6>
        <small class="text-secondary">
            The <code>Undominated_Drift</code> category (|&beta;NTI| &le; 2 and |RC<sub>bray</sub>| &le; 0.95) represents stochastic ecological drift, weak environmental selection, and birth-death fluctuations. It indicates that no single deterministic selection force dominates community assembly.
        </small>
    </div>

    <div class="footer text-center text-muted py-3">
        <small>Generated with <strong>q2-eco-processes</strong> | LaBiOmics / UMC Microbiome Suite</small>
    </div>
</div>

<script>
    // 1. Donut Chart
    const procLabels = {list(proc_totals.index)};
    const procValues = {list(proc_totals.values)};
    Plotly.newPlot('donutPlot', [{{
        labels: procLabels.map(l => l.replace('_', ' ')),
        values: procValues,
        type: 'pie',
        hole: 0.4,
        marker: {{ colors: ['#2b5876', '#4e4376', '#f39c12', '#e74c3c', '#95a5a6'] }}
    }}], {{ margin: {{ t: 10, b: 10, l: 10, r: 10 }} }});

    // 2. 2D Scatter Plot with Dynamic Real Data
    const scatterTrace = {{
        x: {json.dumps(x_coords)},
        y: {json.dumps(y_coords)},
        text: {json.dumps(sample_labels)},
        mode: 'markers+text',
        textposition: 'top center',
        type: 'scatter',
        marker: {{ size: 10, color: '#2b5876' }}
    }};
    const scatterLayout = {{
        xaxis: {{ title: 'RCbray Distance (-1 to +1)', range: [-1.1, 1.1] }},
        yaxis: {{ title: '&beta;NTI Score', range: [-4, 4] }},
        shapes: [
            {{ type: 'line', x0: -1.1, x1: 1.1, y0: 2, y1: 2, line: {{ color: 'red', dash: 'dash' }} }},
            {{ type: 'line', x0: -1.1, x1: 1.1, y0: -2, y1: -2, line: {{ color: 'blue', dash: 'dash' }} }},
            {{ type: 'line', x0: 0.95, x1: 0.95, y0: -2, y1: 2, line: {{ color: 'orange', dash: 'dot' }} }},
            {{ type: 'line', x0: -0.95, x1: -0.95, y0: -2, y1: 2, line: {{ color: 'purple', dash: 'dot' }} }}
        ],
        margin: {{ t: 10, b: 40, l: 50, r: 10 }}
    }};
    Plotly.newPlot('scatterPlot', [scatterTrace], scatterLayout);
</script>
</body>
</html>
"""
    with open(os.path.join(output_dir, "index.html"), "w") as f:
        f.write(index_html)
