#!/bin/sh

mkdir $1
mkdir $1/$2
mkdir $1/$2/pdbqts
mkdir $1/$2/dlgs
cp $3/$2/*.dlg $1/$2/dlgs/.
cd $1/$2/dlgs
for f in *.dlg; do 
	rec=$(echo $f | cut -c 15-18)
	mkdir $rec
	~/MGLTools-*/bin/pythonsh ../../../write_conformations_from_dlg.py -d $f
	mv *.pdbqt $rec
	mv $rec ../pdbqts/.
	echo 'done with ' $rec
done
cd ..
cd ..
cd ..
