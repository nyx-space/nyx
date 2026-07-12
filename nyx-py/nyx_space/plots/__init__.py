from os.path import basename

import polars as pl

from nyx_space.time import Epoch

TEMPLATE = "seaborn"

__all__ = ["od", "md"]


def convert_units(df):
    rename_dict = {}
    exprs = []
    for col in df.columns:
        if "(km/s)" in col:
            new_col = col.replace("(km/s)", "(m/s)")
            rename_dict[col] = new_col
            exprs.append(pl.col(col) * 1000)
        elif "(km^2)" in col:
            new_col = col.replace("(km^2)", "(m^2)")
            rename_dict[col] = new_col
            exprs.append(pl.col(col) * 1_000_000)
        elif "(km)" in col:
            new_col = col.replace("(km)", "(m)")
            rename_dict[col] = new_col
            exprs.append(pl.col(col) * 1000)

    if exprs:
        df = df.with_columns(exprs).rename(rename_dict)
    return df


def watermark(fig, filepath: str):
    """
    Adds a vertical text watermark on the right side of the plot containing
    'Nyx Space', the filename, and the current datetime.
    """
    now_str = Epoch.system_now().strftime("%Y-%m-%d %H:%M:%S %T")
    if filepath is None:
        text = f"{now_str} | Nyx Space"
    else:
        text = f"{basename(filepath)} | {now_str} | Nyx Space"

    fig.add_annotation(
        text=text,
        x=1.0,
        y=0.5,
        xref="paper",
        yref="paper",
        xanchor="left",
        yanchor="middle",
        showarrow=False,
        textangle=90,
        font=dict(size=12, color="gray"),
        xshift=15,
    )
    fig.update_layout(margin=dict(r=80))
    return fig
