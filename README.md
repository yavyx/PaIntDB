# PaIntDB
Pseudomonas aeruginosa Interactions Database

PaIntDB takes a list of genes as input and generates a network of protein-protein interactions with these genes. You need Python 3 and `pip` installed to install and run the application locally. To install the app on Windows, you also need Microsoft Visual C++ Build Tools, this is a good [installation tutorial](https://www.scivision.co/python-windows-visual-c-14-required/).

This is a pre-release version and still contains many bugs. The vis module pre-loads a network from the `temp_data` directory, but you can change the path in line 61 to visualize any other network generated with PaIntDB.  

## Installation
1. Clone repo.
2. Navigate to the project root.
3. Create a new Python virtual environment:  `python -m venv env`
4. Activate virtual environment: `source env/bin/activate` on MacOS/Linux, `.\env\Scripts\activate` on Windows.
5. Download and install required libaries: `pip install -r requirements.txt`
6. To run the network generation app: `python -m dash_app.index`, and go to [http://127.0.0.1:8050/home](http://127.0.0.1:8050/home).
