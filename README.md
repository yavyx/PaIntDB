# PaIntDB
*Pseudomonas aeruginosa* Interactions Database

PaIntDB contains more than 150K protein-protein (PPI) and protein-metabolite interactions in *P. aeruginosa*. It allows 
the integration and visualization of high-throughput (proteomics, RNASeq and/or TnSeq) experimental results, using a 
list of genes as input and mapping them onto a network using the PPI data. It is available at www.paintdb.ca.

You need Python 3 and `pip` installed to install and run the application locally. To install the app on Windows, you also need Microsoft Visual C++ Build Tools, this is a good [installation tutorial](https://www.scivision.co/python-windows-visual-c-14-required/).

## Installation
1. Clone repo.
2. Navigate to the project root.
3. Create a new Python virtual environment:  `python -m venv env`
4. Activate virtual environment: `source env/bin/activate` on MacOS/Linux, `.\env\Scripts\activate` on Windows.
5. Download and install required libaries: `pip install -r requirements.txt`
6. Run app: `python -m dash_app.index`, and go to [http://127.0.0.1:8050/home](http://127.0.0.1:8050/home).
