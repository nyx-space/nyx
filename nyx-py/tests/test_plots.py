import polars as pl
from nyx_space.plots.od import cr_cd, od_dashboard, residuals, uncertainty


def test_residuals_and_dashboard_ground_station():
    # Build a synthetic DataFrame simulating a ground station OD output
    epochs = [
        "2025-01-01T00:00:00Z",
        "2025-01-01T00:01:00Z",
        "2025-01-01T00:02:00Z",
        "2025-01-01T00:03:00Z",
        "2025-01-01T00:04:00Z",
        "2025-01-01T00:05:00Z",
    ]

    df = pl.DataFrame(
        {
            "Epoch (UTC)": epochs,
            "Tracker": ["DSS34", "DSS34", "DSS34", "DSS65", "DSS65", "DSS65"],
            "Prefit residual: Range (km)": [0.01, -0.02, 0.015, -0.01, 0.005, 0.02],
            "Postfit residual: Range (km)": [
                0.001,
                -0.002,
                0.0015,
                -0.001,
                0.0005,
                0.002,
            ],
            "Measurement noise: Range (km)": [0.005] * 6,
            "Whitened residual #0": [0.2, -0.4, 0.3, -0.2, 0.1, 0.4],
            "Residual Rejected": [False] * 6,
            "Sigma X (RIC) (km)": [0.01] * 6,
            "Sigma Y (RIC) (km)": [0.02] * 6,
            "Sigma Z (RIC) (km)": [0.03] * 6,
            "Sigma Vx (RIC) (km/s)": [1e-5] * 6,
            "Sigma Vy (RIC) (km/s)": [2e-5] * 6,
            "Sigma Vz (RIC) (km/s)": [3e-5] * 6,
            "cr": [1.2, 1.25, 1.23, 1.24, 1.22, 1.21],
            "Sigma Cr (Earth J2000) (unitless)": [0.05] * 6,
            "cd": [2.1] * 6,
        }
    )

    res_fig = residuals(df)
    assert res_fig is not None
    assert "Range (m) Residuals" in res_fig.layout.annotations[0].text

    dash_figs = od_dashboard(df)
    assert len(dash_figs) == 1
    assert "Range (m)" in dash_figs[0].layout.title.text

    unc_fig = uncertainty(df)
    assert unc_fig is not None

    crcd_fig = cr_cd(df)
    assert crcd_fig is not None


def test_residuals_and_dashboard_gnss():
    # Build a synthetic DataFrame simulating a GNSS OD output (without Tracker column, X/Y/Z residual columns)
    epochs = [f"2025-01-01T00:0{i}:00Z" for i in range(6)]

    df = pl.DataFrame(
        {
            "Epoch (UTC)": epochs,
            # No Tracker column provided, simulating GPS OD solution
            "Prefit residual: X (km)": [0.01, -0.02, 0.015, -0.01, 0.005, 0.02],
            "Postfit residual: X (km)": [0.001, -0.002, 0.0015, -0.001, 0.0005, 0.002],
            "Measurement noise: X (km)": [0.005] * 6,
            "Prefit residual: Y (km)": [0.02, -0.01, 0.025, -0.02, 0.015, 0.01],
            "Postfit residual: Y (km)": [0.002, -0.001, 0.0025, -0.002, 0.0015, 0.001],
            "Measurement noise: Y (km)": [0.005] * 6,
            "Prefit residual: Z (km)": [0.03, -0.03, 0.035, -0.03, 0.025, 0.03],
            "Postfit residual: Z (km)": [0.003, -0.003, 0.0035, -0.003, 0.0025, 0.003],
            "Measurement noise: Z (km)": [0.005] * 6,
            "Whitened residual #0": [0.2, -0.4, 0.3, -0.2, 0.1, 0.4],
            "Whitened residual #1": [0.4, -0.2, 0.5, -0.4, 0.3, 0.2],
            "Whitened residual #2": [0.6, -0.6, 0.7, -0.6, 0.5, 0.6],
            "Residual Rejected": [False] * 6,
            "Sigma X (RIC) (km)": [0.01] * 6,
            "Sigma Y (RIC) (km)": [0.02] * 6,
            "Sigma Z (RIC) (km)": [0.03] * 6,
            "Sigma Vx (RIC) (km/s)": [1e-5] * 6,
            "Sigma Vy (RIC) (km/s)": [2e-5] * 6,
            "Sigma Vz (RIC) (km/s)": [3e-5] * 6,
        }
    )

    res_fig = residuals(df)
    assert res_fig is not None
    annotation_titles = [
        a.text for a in res_fig.layout.annotations if "Residuals" in a.text
    ]
    assert len(annotation_titles) == 3
    assert "X (m) Residuals" in annotation_titles
    assert "Y (m) Residuals" in annotation_titles
    assert "Z (m) Residuals" in annotation_titles

    dash_figs = od_dashboard(df)
    assert len(dash_figs) == 3
    titles = [fig.layout.title.text for fig in dash_figs]
    assert any("X (m)" in t for t in titles)
    assert any("Y (m)" in t for t in titles)
    assert any("Z (m)" in t for t in titles)
    # Ensure none of the dashboard titles fell back to 'Whitened residual #0'
    assert not any("Whitened residual" in t for t in titles)
