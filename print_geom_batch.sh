mkdir $1
mkdir $1/$2
mkdir $1/$2/geom_props
if [[ $1 = "pdbs" ]]; then 
	for rec in $3/receptors/*.pdb; do 
		filename=$(basename $rec)
		recname=$(echo $filename| cut -c 1-4)
		mkdir $1/$2/geom_props/$recname
		cat $3/$2/*${recname}*.pdb $rec > $1/$2/geom_props/$recname/${recname}_chlorobenzene_0.pdb
	done
	counter=0
	for f in $1/$2/geom_props/*; do 
		python hx_bond_geom.py 'pdb' $f
#		if [[ "$counter" -lt 1 ]]; then
#			cat $f/dist_angles_H.txt > $1/$2/geom_props/all_HXs.txt 
#		else
#			tail -n +2 $f/dist_angles_H.txt >> $1/$2/geom_props/all_HXs.txt
#		fi
#		counter=$((counter+1))
	done 
else 
	cp -r $1/$2/pdbqts/* $1/$2/geom_props/. 
	for f in $1/$2/geom_props/*; do 
		fname=$(basename $f)
		for pose in $f/*.pdbqt; do 
			posename=$(basename $pose .pdbqt)
			cat $pose $3/receptors/$fname*.pdb > $f/$fname'_'$posename.pdbqt
			rm $pose
		done 
	done
	
	for f in $1/$2/geom_props/*; do
                python hx_bond_geom_pose.py 'pdbqt' $f
	done
fi



