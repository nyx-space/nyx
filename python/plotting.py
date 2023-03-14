import pyarrow.parquet as pq
import pandas as pd
import plotly.express as px
from sys import argv

def plot_measurements(path):
    '''
    Plots the provided measurement parquet file
    '''
    arc_pq = pq.read_table(path)

    df = arc_pq.to_pandas()
    # Set the epoch to a Python object
    df['Epoch (UTC)'] = pd.to_datetime(df['Epoch:Gregorian UTC'])

    # Plot all of the columns apart from the default ones.
    default = ['Epoch:Gregorian UTC', 'Epoch:Gregorian TAI', 'Epoch:TAI (s)', 'Tracking device']

    columns = [col for col in arc_pq.column_names if col not in default]

    for col in columns:
        # Create Plotly trace with units as axis label
        fig = px.scatter(df, x='Epoch (UTC)', y=col, color='Tracking device')

        # See if there is any metadata
        try:
            axis_title = '{} ({})'.format(col, ", ".join([f"{k.decode('utf8')}={v.decode('utf8')}" for k, v in arc_pq.field(col).metadata.items()]))
        except AttributeError:
            axis_title = col

        fig.update_layout(
            yaxis_title=axis_title,
            title=f"{col} over time"
        )

        fig.show()

if __name__ == '__main__':
    plot_measurements(argv[-1])