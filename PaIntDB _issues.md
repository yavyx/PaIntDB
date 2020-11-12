# PaIntDB issues


### Back-end
- Make queries into dataframes instead of dictionaries
- Write tests! (check parameters in all main (user input) functions.
- Adding metabolites is too slow, need to check why.
- Come up with first-order interaction algorithms. Need to think more about this algorithm. Maybe
- Double lists when extracting study items.
- Change gene_list to genes_df for all networks

### Dash GUI
- Error messages when 0 genes are mapped.

## Vis module
- Padj decimal issues
- Change node selection to node ids, not node labels
- Fix PA14 GO term enrichment
- Fix subnetworks in flat lists ()
- Fix down regulated genes in first-order networks



## GO enrichment
- Change directory for PA14/PAO1 mapping in enrichment script.
- Create GO term sets for PA14 instead of mapping
- Clean it up, lots of unused code.



### DONE!
- Make DB path a constant.
- Remove source info from network, only leave experimental info.
- Need to create a PPI network class. Have too many functions using the same parameters.
- Check SQL queries to ensure all info is getting retrieved properly. Remove metabolites during query if they're not needed (unnecessary memory usage).
- Add metabolite info to the networks.
- Removed metabolites from gene count.
- Metabolites are not added in computational networks. Not a bug! (There are no computational p-m interactions)
- Download button not working properly.
- IMPORTANT: Inherited classes (DE network and combined network calling init more than once)
- Fix combined networks gene count
- Check TnSeq genes in combined networks
- Only keep enriched terms, not purified
- How to handle exceptions:
    - CSV input:
        - Check no extra columns
        - Check file format
- Save mapped dictionary to JSON file.
- Use os module for paths

### Separate
Assigning all parent GO terms for Pseudomonas genes (to make subnetworks)




