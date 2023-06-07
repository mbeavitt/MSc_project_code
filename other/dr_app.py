# Import statements
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP
import plotly.graph_objs as go

import pickle
import dash
from dash import dcc, html, Dash
from dash.dependencies import Input, Output, State

"""
Performs dimensionality reduction (DR) on a dataframe and visualizes the result. 

Arguments:
features_cols (list): A list of strings specifying the columns in 'df' to be considered as features for dimensionality reduction.
label_col (str): A string specifying the column in 'df' to be used as the label.
df (pandas.DataFrame): The dataframe on which to perform dimensionality reduction.
technique (str, optional): The dimensionality reduction technique to use. Can be 'pca' (default), 'tsne', or 'umap'.
title (str, optional): The title for the plot. If None, no title will be displayed.
plot_3d (bool, optional): If True, a 3D plot will be created. If False (default), a 2D plot will be created.
hover_cols (list, optional): A list of strings specifying additional columns in 'df' to include in hover data of the plotly graph. 

Returns:
Nothing. This function performs the operation and displays the resulting plot but does not return any value.

Raises:
Error: If a specified column in 'features_cols', 'label_col', or 'hover_cols' is not in 'df'.
Error: If the specified 'technique' is not recognized.
"""

# Dimensionality Reduction Function
def dr_labels(features_cols, label_col, df, technique='pca', title=None, plot_3d=False, hover_cols=None, point_size=3):

    scale_factor = 1.3

    for col in features_cols + [label_col] + (hover_cols if hover_cols else []):
        if col not in df.columns:
            print(f'Error: The column {col} is not in the dataframe')
            return

    df_filtered = df[df[label_col].notnull() & df[features_cols].notnull().all(axis=1)]
    df_standardized = (df_filtered[features_cols] - df_filtered[features_cols].mean()) / df_filtered[features_cols].std()

    num_components = 3 if plot_3d else 2

    if technique == 'pca':
        dr = PCA(n_components=num_components)
    elif technique == 'tsne':
        dr = TSNE(n_components=num_components)
    elif technique == 'umap':
        dr = UMAP(n_components=num_components, n_neighbors=20, min_dist=0.6, spread=2, n_jobs=1)
    else:
        print(f'Error: The technique {technique} is not recognized')
        return

    dr_result = dr.fit_transform(df_standardized)

    if plot_3d:
        dr_df = pd.DataFrame(data=dr_result, columns=['Component 1', 'Component 2', 'Component 3'])
    else:
        dr_df = pd.DataFrame(data=dr_result, columns=['Component 1', 'Component 2'])

    dr_df[label_col] = df_filtered[label_col].astype('category').values  # cast the label_col to categorical

    # Sort the unique labels in ascending order
    unique_labels = sorted(dr_df[label_col].unique())
    color_palette = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']

    # Then generate color_dict
    color_dict = {str(label): color_palette[i % len(color_palette)] for i, label in enumerate(unique_labels)}


    if hover_cols:
        hover_data = df_filtered[hover_cols].reset_index(drop=True)
        for col in hover_cols:
            dr_df[col] = hover_data[col].values
    else:
        hover_data = None

    if plot_3d:
        fig = go.Figure() # create a new figure
        for label in unique_labels:
            df_filtered = dr_df[dr_df[label_col] == label]
            hover_data_filtered = hover_data[dr_df[label_col] == label]
            
            hovertemplate_parts = ["<br>".join([f"{col}: %{{customdata[{i}]}}" for i, col in enumerate(hover_data.columns)])]
            
            fig.add_trace(go.Scatter3d(x=df_filtered['Component 1'], y=df_filtered['Component 2'], z=df_filtered['Component 3'],
                                mode='markers', name=str(label),
                                marker=dict(color=color_dict[str(label)], size=point_size),
                                customdata=hover_data_filtered.values,
                                hovertemplate=hovertemplate_parts[0]))

        fig.update_layout(width=800*scale_factor, height=600*scale_factor)
        fig.update_layout(title=title, autosize=True)

    else:
        # generate a plotly scatter plot figure
        fig = go.Figure()
        for label in unique_labels:
            df_filtered = dr_df[dr_df[label_col] == label]
            fig.add_trace(go.Scatter(x=df_filtered['Component 1'], y=df_filtered['Component 2'],
                                    mode='markers', name=str(label),
                                    marker_color=color_dict[str(label)], size=point_size)) # specify marker color
            
        fig.update_layout(title=title, autosize=True,
                          xaxis_title="Component 1", 
                          yaxis_title="Component 2",
                          width=800*scale_factor, height=600*scale_factor)
    # return the figure object instead of showing it
    return fig


# Loading the dataframe from a pickle file
with open('/Users/michaelbeavitt/Desktop/Code/Python/msc_proj/master_df.pickle', 'rb') as f:
    master_df = pickle.load(f)

# Creating a table of contents!
column_dict = {
    "ungated_flc": [],
    "pop1_flc": [],
    "pop2_flc": [],
    "pop3_flc": [],
    "ROS": [],
    "cell_count": [],
    "cytokine": [],
    "biometric": [],
    "netosis": []
}
# Define the start and end column for each category in a list of tuples
range_list = [
    ("1-Unstim-Cells-FoP", "2-PMA-Pop3-PMNs-Hmox1pos-Hmox1-MFI"),
    ("1-FoP", "2-%ROS-lo"),
    ("RBC", "Gra%"),
    ("CD163 (BR28) (28)", "TFR(BR13) (13)"),
    ("Initial_weight", "End_haem"),
    ("media_netosis", "nts_netosis")
]
# For each start and end column, get all column names in this range and store in the dictionary
for i, (start, end) in enumerate(range_list):
    start_idx = list(master_df.columns).index(start)
    end_idx = list(master_df.columns).index(end)

    if i == 0:  # The first range contains multiple categories
        column_dict["ungated_flc"] = [col for col in master_df.columns[start_idx : end_idx+1] if 'Pop' not in col]
        column_dict["pop1_flc"] = [col for col in master_df.columns[start_idx : end_idx+1] if 'Pop1' in col]
        column_dict["pop2_flc"] = [col for col in master_df.columns[start_idx : end_idx+1] if 'Pop2' in col]
        column_dict["pop3_flc"] = [col for col in master_df.columns[start_idx : end_idx+1] if 'Pop3' in col]
    else:  # All other ranges correspond to one category
        key = list(column_dict.keys())[i+3]  # Skip the first four keys
        column_dict[key] = [col for col in master_df.columns[start_idx : end_idx+1]]

label_colnames = [
    'geo_cluster',      # geographic clusters of samples
    'Group',             # one of three; Infected, Resolved, Control
    'Initial_RDT',       # initial rapid diagnostic malaria test (2nd visit)
    'Final_RDT',         # final rapid diagnostic malaria test (2nd visit)
    'Village',           # Village that patient lived in at time of sampling
    'Age',               # Age in years
    'Sex',               # male or female
    'killing_label',     # One of four; 'NonKiller', 'Killer-NonKiller', 'NonKiller-Killer', or 'Killer', based on Salmonella killing assay
    'days_after_start',   # days since start of sampling
    'visit_1_parasites_pa', # presence/absence of parasites on visit 1
    'visit_2_parasites_pa', # presence/absence of parasites on visit 2
    'ID' # Careful using this one to label clusters..!
]


'''
### Guide: ###

fetch column names (colnames) using "column_dict[key]"

fetch actual columns using master_df[colnames]

## keys: ##

"ungated_flc": Columns representing ungated flow cytometry data.gated 
"pop1_flc": Columns representing gated flow cytometry data for population 1 (inactivated(?)).
"pop2_flc": Columns representing gated flow cytometry data for population 2. (activated(?))
"pop3_flc": Columns representing gated flow cytometry data for population 3 (dead). 
"ROS": Columns representing reactive oxygen species (ROS) assay data.
"cell_count": Columns representing cell count data (RBC, granulocytes, plus some others).
"cytokine": Columns representing cytokine data.
"biometric": Columns representing biometric data, which includes weight, height, oxygen saturation and hemoglobin levels..
"netosis": Columns representing netosis data, a type of cell death specific to neutrophils. These are all columns from the sixth range.

You can add them together using '+' e.g.:

master_df[column_dict['pop1_flc'] + column_dict['pop2_flc']]


## Labels and hover data: ##

fetch label columns using master_df[<label_colname>]

label colnames are as follows:

'geo_cluster': geographic clusters of samples
'Group': one of three; Infected, Resolved, Control
'Initial_RDT': initial rapid diagnostic malaria test (2nd visit)
'Final_RDT': final rapid diagnostic malaria test (2nd visit)
'Village': Village that patient lived in at time of sampling
'Age': Age in years
'Sex': male or female
'killing_label': One of four; 'NonKiller', 'Killer-NonKiller', 'NonKiller-Killer', or 'Killer', based on Salmonella killing assay
'days_after_start': days since start of sampling

you can also iterate through them, they are stored in the list label_colnames

'''
# Initialize the app
app = dash.Dash(__name__)

# Define the layout
app.layout = html.Div([
    dcc.Dropdown(id='features_cols', options=[{'label': k, 'value': k} for k in column_dict.keys()], multi=True),
    dcc.Dropdown(id='label_col', options=[{'label': l, 'value': l} for l in label_colnames], multi=False),
    dcc.Checklist(id='hover_cols', options=[{'label': l, 'value': l} for l in label_colnames], value=[]),
    dcc.Slider(
        id='point_size',
        min=1,
        max=10,
        value=3,
        step=0.5,
        marks={i: str(i) for i in range(1, 11)}
    ),
    dcc.Graph(id='dr_plot', config={'displayModeBar': False})  # Save the plot object in the layout
])

# Create the initial plot object
fig = go.Figure()

# Define the callback
@app.callback(
    Output('dr_plot', 'figure'),
    [Input('features_cols', 'value'),
    Input('label_col', 'value'),
    Input('hover_cols', 'value'),
    Input('point_size', 'value')]
)
def update_plot(features_cols, label_col, hover_cols, point_size):
    global fig  # Use the global plot object

    try:
        if features_cols is None or label_col is None:
            return fig  # Return the existing plot object if either features_cols or label_col is None

        # Determine the trigger input
        triggered_input = dash.callback_context.triggered[0]['prop_id'].split('.')[0]

        # Update the plot based on the trigger
        if triggered_input in ['features_cols', 'label_col', 'hover_cols']:
            cols = []
            for feature in features_cols:
                cols += column_dict[feature]

            print("About to perform DR")
            updated_fig = dr_labels(features_cols=cols, label_col=label_col, df=master_df, technique='umap', plot_3d=True, hover_cols=hover_cols, point_size=point_size)
            print("DR done")

            fig = updated_fig  # Update the global plot object
        elif triggered_input == 'point_size':
            fig.update_traces(marker=dict(size=point_size))  # Update the point size in the plot object
            fig.update_layout(uirevision='constant')  # Disable plot regeneration on update

        return fig

    except Exception as e:
        print(f"Error occurred: {e}")
        return fig  # Return the existing plot object if an error occurs


# Run the app
if __name__ == '__main__':
    app.run_server(debug=True, use_reloader=False)


