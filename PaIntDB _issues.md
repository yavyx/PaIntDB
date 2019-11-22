# PaIntDB issues


### Back-end
- Make queries into dataframes instead of dictionaries
- Write tests! (check parameters in all main (user input) functions.
- Adding metabolites is too slow, need to check why.
- Use JSON to save useful dictionaries that only need to be created once (PAO1/PA14 mapping, gene/metabolite mapping, edge list dictionary, GO enrichment dictionary). Will probably improve perfomance dramatically.
- Check sources when building the network from a dataframe! Which source stays, which one leaves?
- Need to add an export to CSV function.
- New function to extract all db information from the selected interaction(s). Important for when the Dash application is running.
- Come up with first-order interaction algortihms. Need to think more about this algorithm.
- Double lists when extracting study items.
- Fix GO enrichment. Run GO enrichment on network or on full list (not all original genes are mapped) ? Both
- Use decorators? to modify network generating functions.
- Change gene_list to genes_df for all networks
- Fix mapped genes in combined networks

### Dash GUI
- How to handle exceptions
- State not working in update_download_link
- React errors

## Vis module
- Reset view button not working

### DONE!
- Make DB path a constant.
- Remove source info from network, only leave experimental info.
- Need to create a PPI network class. Have too many functions using the same parameters.
- Check SQL queries to ensure all info is getting retrieved properly. Remove metabolites during query if they're not needed (unnecessary memory usage).
- Add metabolite info to the networks.
- Removed metabolites from gene count.
- Metabolites are not added in computational networks. Not a bug! (There are no computational p-m interactions)
- Download button not working properly.


### Separate
Assigning all parent GO terms for Pseudomonas genes (to make subnetworks)




