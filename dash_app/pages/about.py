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
                    'PaIntDB is a tool to aid in the interpretation and visualization of high-throughput omics results '
                    'and the generation of new hypotheses implied by your data. It takes a list of genes with optional '
                    'expression data and maps the interactions between them using our database of more than '
                    '157,000 protein-protein interactions (PPI) in '
                    '*Pseudomonas aeruginosa* strains PAO1 and PA14 that was compiled from various sources, including '
                    'interactions derived from orthologous proteins in *E. coli*. You can visualize, explore and filter '
                    'the resulting PPI networks to find interesting groups of genes involved in the conditions of study.\n\n'
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
                    '[Zhang (2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3404098/).\n\n'
                    'If you use PaIntDB to analyze your data please cite this publication:\n\n'
                    'Castillo-Arnemann JJ, Solodova O, Dhillon BK, Hancock REW (2021). PaIntDB: Network-based omics integration'
                    ' and visualization using protein-protein interactions in *Pseudomonas aeruginosa*. *Bioinformatics*.'
                    ' doi: [10.1093/bioinformatics/btab363](10.1093/bioinformatics/btab363)\n\n'
                    'The example data is a set of differentially-expressed genes identified through RNA-Seq in a '
                    '*P. aeruginosa* PAO1 relA/spoT knockout mutant vs. wild type. '
                    'The Tn-Seq example data is a random subset of genes, used to illustrate how the '
                    'RNA-Seq/Tn-Seq integration works.\n\n'
                    'For the best experience we recommend using Chrome or Firefox.'
                    'The source code and database file are available in [Github](https://github.com/yavyx/PaIntDB).\n'
                    'If you encounter any bugs, please open an issue in Github.\n\n',
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
