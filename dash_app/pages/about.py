import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from dash_app.app import app

layout = html.Div(
    [
        dbc.Jumbotron(
            [
                html.H1('About'),
                html.Hr(),
                dcc.Markdown(
                    'PaIntDB contains more than 157,000 protein-protein and protein-metabolite interactions in '
                    '*Pseudomonas aeruginosa* strains PAO1 and PA14, compiled from various sources, including '
                    'interactions derived from orthologous proteins in *E. coli*.\n\n'
                    'The source code and database file are available in [Github](https://github.com/yavyx/PaIntDB).\n'
                    'If you encounter any bugs, please open an issue in Github.\n\n'
                    'List of compiled interaction databases and studies:\n'
                    '[APID](http://cicblade.dep.usal.es:8080/APID/init.action), '
                    '[BindingDB](https://www.bindingdb.org/bind/index.jsp), '
                    '[DIP](https://dip.doe-mbi.ucla.edu/dip/Main.cgi), '
                    '[EcoCyc](https://ecocyc.org/), '
                    '[Galan-Vazquez (2011)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3348663/), '
                    '[IMEx](https://www.imexconsortium.org/), '
                    '[IntAct](https://www.ebi.ac.uk/intact/), '
                    '[iRefIndex](https://irefindex.vib.be/wiki/index.php/iRefIndex), '
                    '[KEGG](https://www.genome.jp/kegg/), '
                    '[mentha](https://mentha.uniroma2.it/), '
                    '[MINT](https://mint.bio.uniroma2.it/), '
                    '[MPIDB](https://europepmc.org/article/pmc/pmc2638870), '
                    '[RegulonDB](http://regulondb.ccg.unam.mx/), '
                    '[UniProt](https://www.uniprot.org/), '
                    '[XLinkDB](http://xlinkdb.gs.washington.edu/xlinkdb/), '
                    '[Zhang (2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3404098/).',
                    style={'width': '60vw',
                           'font-size': '20px'}
                )
            ],
            style={'margin': '10px',
                   'backgroundColor': '#a6edff'}
        ),
        html.Div(
            [
                html.A(
                    html.Img(src=app.get_asset_url('hancock-lab-logo.svg'),
                             style={'width': '175px'}),
                    href='http://cmdr.ubc.ca/bobh/'
                ),
                dcc.Markdown('PaIntDB is being developed by the [Hancock Laboratory](http://cmdr.ubc.ca/bobh/) at the ' 
                             'University of British Columbia.\n\nFunding is currently provided by the Canadian '
                             'Institutes of Health Research FDN-154287.',
                             style={'margin-top': '10px'})
            ],
            style={'font-size': '14px',
                   'margin-left': '10px',
                   'padding': '7px'}
        )
    ],
    style={'background-color': '#ededed',
           'height': 'calc(100vh - 76px)',
           'overflow-y': 'auto'}
)
