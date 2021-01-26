# CombiFF script - Collect reference data for force-field optimization.

* Search for identifier (name, InChIKey, CAS registry number) in [PubChem](https://pubchem.ncbi.nlm.nih.gov/) given smiles strings.
* Search for themodynamic, dielectric, transport and solvation data of organic molecules in an in-house database of experimental data of organic compounds.
* For every molecule, plot all available data and let the user select the reference data points.
* Write selected data to csv file.

![](/images/chap_1_vic.png)

[J. Chem. Theory Comput., 16, 7525-​7555 (2020).](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00683)

### The five main options of this script are:

1. Get identifiers

    * Given list of SMILES strings, search for identifiers in [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
    
      If a file with CID (id) and smiles from pubchem is given,
      then the CID is used directly to retrieve the molecule from pubchem.
      
      If an empty file is given,
      then the smiles is used to search for the molecule in the database.
      
      In both cases [RDKit](https://www.rdkit.org/) is used to canonicalize the smiles strings.

2. Write identifiers to csv file

    * Write matched identifiers (name, smiles, standardInChiKey, CAS registry number and molecular foruma)
    to csv file.

3. Search for all availabe reference data in database

    * Search for reference data in the tables.
    * The tables are saved in tmp/ directory.
    * The output matched data is saved in data/ directory.

4. Select data points

    * Plot all available data for each property and each molecule
      so that the user can select the data points that will be used as reference data
      for force-field optimization.

5. Plot data

    * For each molecule and each property
    save plot with all available data and simulated results
    if input file with simulated values is given.
