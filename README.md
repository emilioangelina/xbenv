
# xbenv: Halogen Bond (XB) prediction based in protein Environments 


![figure1](https://github.com/lemyp-cadd/XBenv/blob/main/figure1.jpg)


This module, termed xbenv, is an implementation of the protocol for XB prediction based on halogen environment described in 
**"Targeting Protein Pockets with Halogen Bonds"** (*Bogado ML, Villafañe RN, Gomez Chavez JL, Angelina EL, Sosa GL, Peruchena NM (2022) to be published...*). 

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

`git clone https://github.com/lemyp-cadd/XBenv.git`

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

Current version of the protocol accepts ligand-protein complexes in pdb format or as autodock .dlg files (docking log files). The charge density derived properties have to be provided as well as docking poses. (Protocol implements only steps 8 to 11 of Figure 1) 

The protocol can be run by invoking the main script evaluate.py which takes five arguments:

`python evaluate.py arg1 arg2 arg3 arg4 arg5`

arg1 can take any value. It indicates the folder where the input data is stored, i.e. "data"
By input data we mean: 


	a) Ligands and receptors from X-bonded complexes, have to be provided separatedly in two folders whitin the input data folder. Receptors should be stored in a folder named "receptors". Ligands are stored in a folder whose name is provided by arg4 (see below). Ligands and receptor file names should start with 4 letter code, usually a pdb id, i.e. 3qgt_ligand.pdb and 3qgt_receptor.pdb 

	b) csv files with charge density derived properties from previous survey on X-bonded complexes

	c) numpy binary file parameters.npy containing the adjusted parameters of the machine learning model that predict BCP formation

	d) txt file fit_param.txt with adjusted parameters for 1) BCP formation cuttoff probability (prob_cutoff) and 2) the number of radial and angular bins (nbins) for partitioning the charge density derived data.  



arg2 can take one of two values: "pdbs" or "poses", that indicate whether the ligand structures to be processed are single structures in pdb format  or docking poses in dlg format, respectively. 

arg3 can take one of two values: "retrospective" or "prospective", the first option for a retrospective evaluation of the protocol on the same structures from which the charge density-derived XB environments were extracted. In a prospective setting the input structures are different than those used for obtain the charge density data. 

arg4 indicates the name of the folder containing the input ligand structures within the input data folder, i.e. pdbs_xbonds or dlgs_xbonds 

arg5 can take 4 arguments: "fit", "predict", "evaluate" or any pdb identifier for an individual structure, i.e. "3qgt". 

"fit" performs the ajustements of the parameters prob_cutoff and nbins by minimizing the relative error of the predictions on pdb structures or docking poses. It can be run only in a retrospective setting. 

For instance:

`python evaluate.py "data" "pdbs" "retrospective" "pdbs_xbonds" "fit" ` 

or

`python evaluate.py "data" "poses" "retrospective" "dlgs_xbonds" "fit"`
 
For adjustement of the machine learning parameters, run jupyer notebook predict_BCP.ipynb. Once adjusted, the numpy binary file parameters.npy in input data folder will be updated with the new ML parameters.  

"predict" perform the XB environment prediction with already adjuted parameters. It can be run in retrospective as well as prospective setting in both pdbs or docking poses. For instance: 

`python evaluate.py "data" "pdbs" "retrospective" "pdbs_xbonds" "predict"`  

or

`python evaluate.py "data_prospective" "pdbs" "prospective" "pdbs_xbonds_prospective" "predict"` 

or 

`python evaluate.py "data" "poses" "retrospective" "dlgs_xbonds" "predict"` 

or

`python evaluate.py "data_prospective" "poses" "prospective" "dlgs_xbonds_prospective" "predict"`

"evaluate" compare the computed XB environment with usual XB environments in known X-bonded complexes to get an estimate of how likely is to obtain the predicted value. For instance: 

`python evaluate.py "data" "pdbs" "retrospective" "pdbs_xbonds" "evaluate"` 

or

`python evaluate.py "data" "poses" "retrospective" "dlgs_xbonds" "evaluate"`

Finally, if the pdb id is provided, the analysis of the individual XB environment hydrogen atoms for the specified structure is performed. For instance: 

`python evaluate.py "data" "pdbs" "retrospective" "pdbs_xbonds" "3qgt"`   

or 

`python evaluate.py "data" "poses" "retrospective" "dlgs_xbonds" "3qgt"` 

or 

`python evaluate.py "data" "pdbs" "prospective" "pdbs_xbonds" "3qgt"` 

or 

`python evaluate.py "data" "poses" "prospective" "dlgs_xbonds" "3qgt"`  


An example output is: 


`ID=3qgt,best_pose=0,best_pred_prop=-0.31999108931,true_prop=-0.31999108931`

`('HB2', 'SER', 108, 0.97910425647902366, -0.17115947758737141)`

`('HD12', 'ILE', 112, 0.71028922382027027, -0.14883161172288742)`


For instance, in the last line:  HD12 is atom name, ILE the residue name, 112 the residue number, 0.71028922382027027 is the predicted probability of BCP formation and -0.14883161172288742 is the value of the predicted charge density-derived property

Only the H atoms predicted to form a BCP with chlorine are printed. 


## Citing this work

If you use the code or data in this package, please cite:

```bibtex
@Article{XBenv,
  author  = {Bogado ML, Villafañe RN, Gomez Chavez JL, Angelina EL, Sosa GL, Peruchena NM},
  journal = {ChemRxiv},
  title   = {Targeting Protein Pockets with Halogen Bonds},
  year    = {2021},
  volume  = {---},
  number  = {---},
  pages   = {----},
  doi     = {----}
}
```
