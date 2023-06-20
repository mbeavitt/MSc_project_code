# Import statements
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP
import plotly.graph_objs as go

import pickle
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate

def performDR(features_cols, df, technique='umap'):

    df_filtered = df[df[features_cols].notnull().all(axis=1)]
    df_standardized = (df_filtered[features_cols] - df_filtered[features_cols].mean()) / df_filtered[features_cols].std()

    dr = {'pca': PCA(n_components=3), 'tsne': TSNE(n_components=3), 'umap': UMAP(n_components=3, n_neighbors=20, min_dist=0.6, spread=2)}[technique]
    dr_result = dr.fit_transform(df_standardized)
    dr_df = pd.DataFrame(data=dr_result, columns=['Component 1', 'Component 2', 'Component 3'])

    return dr_df

def createFigure(title, point_size, dr_df, master_df, features_cols, label_col, plot_3d=True, hover_data=None):

    # Removing NA values
    valid_rows = master_df[features_cols].notnull().all(axis=1)

    # Get valid labels
    labels_df = master_df[label_col][valid_rows]

    # Fill NA values
    labels_df = labels_df.fillna("not defined")

    # Convert to list
    labels = labels_df.tolist()

    # Concatenate with IDs
    labels_df = pd.concat([master_df['ID'][valid_rows], labels_df], axis=1)

    # Generate a sorted list of unique labels for factorizing
    sorted_unique_labels = sorted(list(set(labels)))

    # Create a mapping of sorted labels to their indices (factors)
    label_mapping = {label: i+1 for i, label in enumerate(sorted_unique_labels)}

    # Factorize the labels based on the sorted label mapping
    labels_df[label_col] = labels_df[label_col].map(label_mapping)

    # Add 'count' column
    labels_df['count'] = labels_df.groupby(label_col).cumcount() + 1


    # Adding labels to the dimension reduction dataframe and resetting the index so the hover data works too
    dr_df['labels'] = labels
    dr_df = dr_df.reset_index(drop=True)


    if hover_data is not None:

        # Grabbing the hover data from the master df and resetting the indexes to align with dr_df (why am I using dr_df?)
        hover_df = master_df[hover_data][valid_rows]
        hover_df = hover_df.reset_index(drop=True)

    else:
        hover_df = None

    unique_labels = sorted(set(labels))

    fig = go.Figure()

    color_palette = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']
    color_dict = {str(label): color_palette[i % len(color_palette)] for i, label in enumerate(unique_labels)}

    for _, row in labels_df.iterrows():
        print(row.values)

    #print(labels_df)

    for label in unique_labels:

        label_data = dr_df[dr_df['labels'] == label]
        color = color_dict[str(label)]
        if hover_df is not None:
            hover_data_filtered = hover_df.loc[label_data.index]
            hovertemplate_parts = ["<br>".join([f"{col}: %{{customdata[{i}]}}" for i, col in enumerate(hover_data)])]
        else:
            hover_data_filtered = None
            hovertemplate_parts = []

        if plot_3d:
            fig.add_trace(go.Scatter3d(x=label_data['Component 1'], y=label_data['Component 2'], z=label_data['Component 3'],
                                       mode='markers', name=str(label),
                                       marker=dict(color=color, size=point_size),
                                       customdata=hover_data_filtered.values if hover_data_filtered is not None else [],
                                       hovertemplate=hovertemplate_parts[0] if hovertemplate_parts else ''))
            
            fig.update_layout(width=800, height=600, title=title, autosize=True, dragmode = 'select')

        else:
            fig.add_trace(go.Scatter(x=label_data['Component 1'], y=label_data['Component 2'],
                                     mode='markers', name=str(label),
                                     marker=dict(color=color, size=point_size),
                                     customdata=hover_data_filtered.values if hover_data_filtered is not None else [],
                                     hovertemplate=hovertemplate_parts[0] if hovertemplate_parts else '',
                                     selectedpoints = [],
                                     selected = dict(marker = dict(color = 'orange')),
                                     unselected = dict(marker = dict(opacity = 0.9))))

            fig.update_layout(width=800, height=600, title=title, autosize=True, dragmode = 'select')

    return fig, labels_df

# Loading the dataframe from a pickle file
with open('master_df.pickle', 'rb') as f:
    master_df = pickle.load(f)

# Creating a table of contents!
column_dict = {
    "Ungated Flow Cytometry Data": [],
    "Population 1 Flow Data (Inactive Neutrophils)": [],
    "Population 2 Flow Data (Activated Neutrophils)": [],
    "Population 3 Flow Data (Dead Neutrophils)": [],
    "ROS assay Data": [],
    "Cell Count Data": [],
    "Cytokine Quantification Assay Data": [],
    "Biometric Data (Height, weight, etc.)": [],
    "Netosis Assay Data": [],
    "Netosis Assay Data (fold change)": []
}
# Define the start and end column for each category in a list of tuples
range_list = [
    ("1-Unstim-Cells-FoP", "2-PMA-Pop3-PMNs-Hmox1pos-Hmox1-MFI"),
    ("1-FoP", "2-%ROS-lo"),
    ("RBC", "Gra%"),
    ("CD163 (BR28) (28)", "TFR(BR13) (13)"),
    ("Initial_weight", "End_haem"),
    ("media_netosis", "nts_netosis"),
    ("PMA_fc", "NTS_fc")
]
# For each start and end column, get all column names in this range and store in the dictionary
for i, (start, end) in enumerate(range_list):
    start_idx = list(master_df.columns).index(start)
    end_idx = list(master_df.columns).index(end)

    if i == 0:  # The first range contains multiple categories
        column_dict["Ungated Flow Cytometry Data"] = [col for col in master_df.columns[start_idx : end_idx+1] if 'Pop' not in col]
        column_dict["Population 1 Flow Data (Inactive Neutrophils)"] = [col for col in master_df.columns[start_idx : end_idx+1] if 'Pop1' in col]
        column_dict["Population 2 Flow Data (Activated Neutrophils)"] = [col for col in master_df.columns[start_idx : end_idx+1] if 'Pop2' in col]
        column_dict["Population 3 Flow Data (Dead Neutrophils)"] = [col for col in master_df.columns[start_idx : end_idx+1] if 'Pop3' in col]
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
    'ID', # Careful using this one to label clusters..!
    'RNAseq_done',
    'Netosis_responder'
]

colour_colnames = [
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
    'RNAseq_done',
    'Netosis_responder'
]

label_to_display_map = {
    'geo_cluster': 'Geo Cluster',
    'Group': 'Infection Status',
    'Initial_RDT': 'Initial Malaria RDT',
    'Final_RDT': 'Final Malaria RDT',
    'Village': 'Village',
    'Age': 'Age',
    'Sex': 'Sex',
    'killing_label': 'Killing Assay Label',
    'days_after_start': 'Days Since Start',
    'visit_1_parasites_pa': 'Parasites at Start? (y/n)',
    'visit_2_parasites_pa': 'Parasites at End? (y/n)',
    'ID': 'Patient ID',
    'RNAseq_done': 'RNA-seq Data? (y/n)',
    'Netosis_responder': 'Netosis_responder'
}

# Initialize the app
app = dash.Dash(__name__)
application = app.server

# Define the layout
app.layout = html.Div([

    html.H1("UMAP dimensionality reduction of the AMAN (Asymptomatic Malaria, Anemia and Neutrophils) dataset", style={'textAlign': 'center'}),

    html.Div([  # Div for hover_cols
        html.Label('Hover Data: (optional)'),
        dcc.Checklist(id='hover_cols', options=[{'label': label_to_display_map[l], 'value': l} for l in label_colnames], value=[],
                      style={'float': 'left', 'width': '100%'})],
        style={'width': '13%', 'float': 'left'}
    ),

    html.Div([  # Div for features_cols and label_col
        
        html.Div([  # Div for plot_3d
            html.Label('Choose 2D/3D: (optional)'),
            dcc.Checklist(id='plot_3d', options=[{'label': '3D Plot', 'value': True}], value=[])],
            style={'width': '13%', 'float': 'left'}
        ),
        
        html.Div([  # Subdiv for features_cols
            html.Label('Feature selection:'),
            dcc.Dropdown(id='features_cols', options=[{'label': k, 'value': k} for k in column_dict.keys()], multi=True)],
            style={'width': '40%', 'float': 'left', 'margin-right': '20px'}
        ),
        html.Div([  # Div for update_button
            html.Button('Update Plot', id='update_button')],
            style={'float': 'left', 'margin-top': '25px', 'margin-right': '50px'}
        ),
        html.Div([  # Subdiv for label_col
            html.Label('Colour by:'),
            dcc.Dropdown(id='label_col', options=[{'label': label_to_display_map[l], 'value': l} for l in colour_colnames], multi=False)],
            style={'width': '23%', 'float': 'left', 'margin-right': '40px'}
        ),

        html.Div([  # Div for dr_plot
            dcc.Graph(id='dr_plot', config={'displayModeBar': False}, style={'width': '100%', 'height': '100vh', 'margin': '0 auto'},
                    responsive=True)],
            style={'clear': 'both'}
        )
    ], style={'width': '70%', 'float': 'left', 'margin-bottom': '20px'}),  # Parent Div for features_cols and label_col

    html.Div([  # Div for point_size
        html.Label('Point Size:'),
        dcc.Slider(id='point_size', min=1, max=10, value=3, step=0.5, marks={i: str(i) for i in range(1, 11)})],
        style={'width': '40%', 'margin': '0 auto'}
    ),

    dcc.Store(id='intermediate_data'),  # the new Store component

    html.Div(id='selected-data',
             style={
                'position': 'fixed',
                'top': '0',
                'right': '0',
                'bottom': '0',
                'width': '150px',
                'padding': '20px',
                'overflow-y': 'auto'
            })

], style={'width': '100%'})


# Create the initial plot object
fig = go.Figure()

@app.callback(
    Output('intermediate_data', 'data'),
    [Input('features_cols', 'value'),
     Input('update_button', 'n_clicks')]
)
def update_intermediate_data(features_cols, n_clicks):
    if n_clicks is None:
        raise PreventUpdate

    cols = []
    for feature in features_cols:
        cols += column_dict[feature]

    plot_dataframe = performDR(cols, master_df)
    
    data = plot_dataframe.to_dict('records')
    
    # Store both the data and features_cols in the dcc.Store
    return {"data": data, "features_cols": cols}


@app.callback(
    Output('selected-data', 'children'),
    Input('dr_plot', 'selectedData')
)
def display_selected_data(selectedData):
    global labels_df

    if selectedData is not None:
        label_ids = [point['curveNumber']+1 for point in selectedData['points']]
        counts = [point['pointIndex']+1 for point in selectedData['points']]

        #print(label_ids)
        #print(counts)

        id_list = []

        for label_id, count in zip(label_ids, counts):
            id_list.extend(labels_df[(labels_df.iloc[:, 1] == label_id) & (labels_df['count'] == count)]['ID'])

        print(id_list)
        return [html.P(f"ID: {str(id)}") for id in id_list]
    return []


@app.callback(
    Output('dr_plot', 'figure'),
    [Input('intermediate_data', 'data'),
     Input('label_col', 'value'),
     Input('hover_cols', 'value'),
     Input('point_size', 'value'),
     Input('plot_3d', 'value')]
)
def update_plot(intermediate_data, label_col, hover_cols, point_size, plot_3d):

    global fig
    global labels_df

    if intermediate_data is None or label_col is None:
        return fig

    data = intermediate_data['data']  # Retrieve the data from intermediate_data
    features_cols = intermediate_data['features_cols']  # Retrieve the features_cols from intermediate_data

    plot_dataframe = pd.DataFrame(data)

    updated_fig, labels_df = createFigure(plot_3d=plot_3d, label_col=label_col, title=None, point_size=point_size, dr_df=plot_dataframe, master_df=master_df, features_cols=features_cols, hover_data=hover_cols)
    fig = updated_fig

    return fig

# Run the app
if __name__ == '__main__':
    application.run(host='0.0.0.0', debug=True, use_reloader=False, port=8080)

