
# xbenv: Halogen Bond (XB) prediction based in protein Environments 


![figure1](https://github.com/lemyp-cadd/XBenv/blob/main/figure1.jpg)


This module, termed xbenv, is an implementation of the protocol for XB prediction based on halogen environment described in 
**"Targeting Protein Pockets with Halogen Bonds: The role of the halogen environment"** (*Bogado ML, Villafañe RN, Gomez Chavez JL, Angelina EL, Sosa GL, Peruchena NM (2022) to be published...*). 

-------------
INSTALLATION: 
-------------

Current version of xbenv works under python2. You can setup a python2 environment in conda with the required dependencies. 

Required dependencies:
```
* numpy  
* scipy 
* pandas 
* matplotlib
```

## install necessary packages:

OS: Linux

[MGLTools](https://ccsb.scripps.edu/mgltools/downloads/)

[Anaconda](https://www.anaconda.com/products/individual)


## Download the package and then enter the folder:  

`git clone https://github.com/emilioangelina/xbenv.git`

`cd xbenv`

## create a python 2 pearsonal conda environment

`conda create -n xbenv python=2.7`

`conda activate xbenv`

## dependencies

`conda install scipy numpy pandas matplotlib`

or 

`pip2 install scipy numpy pandas matplotlib`

Verify that the environment uses version 2.7 of the numpy library, if not, it is recommended to remove the default version and install the appropriate version again.
Or alternatively, install the packages through the environment file: 

`conda env create -f xbenv_env.yaml` 

----- 
USAGE: 
-----

Current version of the protocol accepts ligand-protein complexes in pdb format or as autodock .dlg files (docking log files). The charge density derived properties have to be provided, as well as docking poses. (Protocol implements only steps 8 to 11 of Figure 1) 

The protocol can be run by invoking the main script evaluate.py which takes 3 arguments:

`python evaluate.py arg1 arg2 arg3` 

arg1 can take any value. It indicates the folder where the input data is stored, i.e. "data"
By input data we mean: 


	a) Ligands and receptors from X-bonded complexes, have to be provided separatedly in two folders whitin the input data folder. Receptors should be stored in a folder named "receptors". Ligands are stored in a folder whose name is provided by arg4 (see below). Ligands and receptor file names should start with 4 letter code, usually a pdb id, i.e. 3qgt_ligand.pdb and 3qgt_receptor.pdb 

	b) csv files with charge density derived properties from previous survey on X-bonded complexes

	c) numpy binary file parameters.npy containing the adjusted parameters of the machine learning model that predict BCP formation

	d) txt file fit_param.txt with adjusted parameters for 1) BCP formation cuttoff probability (prob_cutoff) and 2) the number of radial and angular bins (nbins) for partitioning the charge density derived data.  


arg2 can take one of two values: "pdbs" or "poses", that indicate whether the ligand structures to be processed are single structures in pdb format  or docking poses in dlg format, respectively. 

arg3 indicates the name of the folder containing the input ligand structures within the input data folder, i.e. pdbs_ligands or dlgs_ligands 


 For instance for pdb structures: 

`python evaluate.py "data" "pdbs" "pdbs_ligands"`


or for docking poses: 

`python evaluate.py "data "poses" "dlgs_ligands"` 


Only the H atoms predicted to form a BCP with chlorine are printed. 


## Citing this work

If you use the code or data in this package, please cite:

```bibtex
@Article{XBenv,
  author  = {Bogado ML, Villafañe RN, Gomez Chavez JL, Angelina EL, Sosa GL, Peruchena NM},
  journal = {ChemRxiv},
  title   = {Targeting Protein Pockets with Halogen Bonds: The role of the halogen environment},
  year    = {2021},
  volume  = {---},
  number  = {---},
  pages   = {----},
  doi     = {----}
}
```
