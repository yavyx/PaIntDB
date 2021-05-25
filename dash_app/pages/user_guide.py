import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from dash_app.app import app

data_upload_tab = dbc.Card(
    dbc.CardBody(
        [
            dcc.Markdown('You can create protein-protein interaction networks in 4 sequential steps, detailed below.'),
            dbc.CardImg(src=app.get_asset_url('tutorial_1.png'),
                        style={'width': '60vw',
                               'padding-bottom': '10px'}),
            html.Br(),
            dcc.Markdown(
                '#### 1. Data upload\n'
                'PaIntDB takes a list of  of *Pseudomonas aeruginosa* locus tags, which must be in the first column'
                ' of a .csv file. The three upload options work as follows:\n'
                '- Gene list: any additional columns are ignored.\n'
                '- Differentially-expressed gene list: the column with fold changes must be named '
                '\'log2FoldChange\' and the column with p-values must be named \'padj\'. Any additional columns are '
                'ignored. \n'
                '- Combined: A differentially-expressed gene list, with the same requirements as above, and '
                'another .csv file with locus tags in the first column.\n'
                '#### 2. Network generation options\n'
                'Choose the parameters to generate the network.\n'
                '- Strain: *P. aeruginosa* PAO1 or  PA14, depending on your data.\n'
                '- Network order:\n'
                '   - Zero-order: Maps direct interactions between the queried genes. '
                'Recommended for long lists (>200 genes).\n'
                '   - First-order: Uses your queried genes as "seed" genes and finds any interaction between them '
                'and the other genes in the database. Recommended for short lists (<200 genes). Warning: When used '
                'with longer lists this option results in very large networks and performance is very slow.\n'
                '- Interaction detection method:\n' 
                '   - All: includes interactions predicted computationally.\n'
                '   - Experimental:  only includes interactions verified experimentally, resulting in smaller '
                'networks.\n'
                '#### 3. GO term enrichment\n'
                'You can select to run the enrichment using all the input genes or only the genes that were mapped to '
                'the network. Fischer\'s exact test is used with a 0.05 p-value cutoff. The full enrichment results '
                'can be downloaded as a .csv file.\n'
                '#### 4. Explore network\n'
                'After enrichment is run, it is ready for visualization and analysis. Click the explore'
                ' button to continue.',
                style={'width': '70vw',
                       'font-size': '16px'}
            ),
        ],
    ),
    style={'background-color': '#ededed'}
)

explore_tab = dbc.Card(
    dbc.CardBody(
        [
            dbc.CardImg(src=app.get_asset_url('tutorial_2.png'),
                        style={'width': '60vw',
                               'padding-bottom': '10px'}),
            html.Br(),
            dcc.Markdown(
                '#### A. Color mapping\n'
                'If differential expression data is included, the node color indicates up- or down-regulation '
                'to identify co-regulated gene modules in the network. If Tn-Seq genes are included, the color mapping '
                'can be changed to indicate the experiment in which the genes were identified. Node labels can be '
                'toggled on or off.\n'
                '#### B. Select nodes\n'
                'Nodes can be selected according to their cellular location or the enriched GO terms containing them. '
                'If differential expression data is included, then it is possible to select up-regulated or '
                'down-regulated genes, and if a Tn-Seq dataset is included, there is another filter to select genes '
                'according to the experiment in which they were identified. All filters can be combined to fine-tune '
                'the selected nodes as desired. Individual genes of '
                'interest can also be added by name or locus tag to the query.\n\n'
                'Once you select your genes of interest with the filters, you can click on the button to generate a '
                'subnetwork connecting these genes.\n'
                'You can download the network as a '
                '.graphml file for use in [NetworkAnalyst](https://www.networkanalyst.ca) or '
                '[Cytoscape](https://cytoscape.org), '
                'or export it as a .png image.\n'
                '##### Subnetwork Menu\n'
                'After generating a subnetwork, you have the following options:\n\n'
                '- Including additional genes, called Steiner nodes, that are not part '
                'of your original data but have interactions with your genes. '
                'Useful to connect subnetworks with many smaller components.\n'
                '- Including low-confidence interactions in the subnetwork solution, instead of prioritizing '
                'experimental interactions.\n'
                '#### C. Network view\n'
                'You can zoom, pan and select nodes using the mouse. Box selection is possible with Shift + drag. '
                'This selection is independent from the filter selection and cannot be used to build subnetworks.\n'
                '#### D. Node details\n'
                'This table shows the names and descriptions of the selected genes, and differential expression data, '
                'if included. This table can be downloaded as a .csv file.',
                style={'width': '70vw',
                       'font-size': '16px'}
            ),
        ]
    ),
    style={'background-color': '#ededed'}
)

text_style = {'text-align': 'center', 'font-weight': 'bold'}

gallery_tab = dbc.Card(
    dbc.CardBody(
        [
            dbc.Row(
                dbc.Col(
                    dcc.Markdown(
                        'Here are some examples of subnetworks that can be created with PaIntDB using its filtering '
                        'system. These networks were made using the example data of differentially-expressed genes '
                        'identified through RNA-Seq in a *P. aeruginosa* PAO1 relA/spoT knockout mutant vs. wild type. '
                        'The Tn-Seq example data is a random subset of genes, used to illustrate how the '
                        'RNA-Seq/Tn-Seq integration works. The selected genes are '
                        'connected using the using the smallest number of additional nodes possible to generate the '
                        'subnetworks. The image title indicates the filter used.')
                ),
            ),
            html.Br(),
            dbc.Row([
                dbc.Col(
                    [
                        html.P('Experiment', style=text_style),
                        html.Img(src=app.get_asset_url('sig_source_legend.svg'), width='100%',
                                 style={'padding-bottom': '10px'}),
                        html.P('Regulation', style=text_style),
                        html.Img(src=app.get_asset_url('de_legend.svg'), width='100%')
                    ],
                    width=1),
                dbc.Col(
                    [
                        html.P('GO term: Protein Secretion', style=text_style),
                        html.Img(src=app.get_asset_url('protein_secretion_network.png'), width='100%')
                    ],
                    width=4),
                dbc.Col(
                    [
                        html.P('Experiment: Tn-Seq', style=text_style),
                        html.Img(src=app.get_asset_url('tn_seq_network.png'), width='100%')
                    ],
                    width=5),
            ],
                justify='center'
            ),
            dbc.Row([
                dbc.Col(
                    [
                        html.P('Localization: Outer Membrane Vesicle', style=text_style),
                        html.Img(src=app.get_asset_url('outer_membrane_vesicle.png'), width='100%')
                    ],
                    width=5),
                dbc.Col(
                    [
                        html.P('Localization: Cytoplasmic Membrane', style=text_style),
                        html.Img(src=app.get_asset_url('membrane_network.png'), width='100%')
                    ],
                    width=5)
            ],
                justify='center'
            )
        ]
    ),
    style={'background-color': '#ededed'}
)


layout = html.Div(
    [
        dbc.Jumbotron(
            [
                html.H1('User Guide'),
                dbc.Tabs(
                    [
                        dbc.Tab(data_upload_tab, label='Build Network'),
                        dbc.Tab(explore_tab, label='Explore Network'),
                        dbc.Tab(gallery_tab, label='Gallery'),
                    ],
                ),
            ],
            style={'margin': '10px',
                   'backgroundColor': '#a6edff'}
        )
    ]
)