import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

layout = html.Div(
    [
        dbc.Jumbotron(
            [
                dcc.Markdown(
                    '## Usage\n'
                    '***\n'
                    '#### 1. Data upload\n'
                    'PaIntDB takes a list of  of *Pseudomonas aeruginosa* locus tags, which must be in the first column'
                    ' of a .csv file. The three upload options work as follows:\n'
                    '- Gene list: any extra columns are ignored.\n'
                    '- Differentially-expressed gene list: the column with fold changes must be named '
                    '\'log2FoldChange\' and the column with p-values must be named \'padj\'. Any extra columns are '
                    'ignored. \n'
                    '- Combined: A differentially-expressed gene list, with the same requirements as above, and '
                    'another .csv file with locus tags in the first column.\n'
                    '#### 2. Network generation options:\n'
                    'Choose the parameters to generate the network.\n'
                    '- Strain: *P. aeruginosa* PAO1 or  PA14, depending on your data.\n'
                    '- Network order:\n'
                    '   - Zero-order: Maps direct interactions between the queried genes.'
                    'Recommended for long lists (>200 genes).\n'
                    '   - First-order: Uses your queried genes as "seed" genes and finds any interaction between them '
                    'and the other genes in the database. Recommended for short lists (<200 genes).\n'
                    '- Interaction detection method:\n'
                    '   - All: includes interactions predicted computationally.'
                    '   - Experimental:  only includes interactions verified experimentally, resulting in smaller '
                    'networks.\n'
                    '#### 3. GO term enrichment\n'
                    'You can select to run the enrichment using all the input genes or the genes that were mapped to '
                    'the network. Fischer\'s exact test is used with a 0.05 p-value cutoff.'
                    'none',
                    style={'width': '80vw',
                           'font-size': '16px'}
                ),
                dbc.Button('Get Started', color='primary', id='start', href='/menu', size='lg')
            ],
            style={'margin': '10px',
                   'backgroundColor': '#d1d1d1'}
        ),
    ]
)