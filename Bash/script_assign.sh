Folder=('4Myr' 'Renzo' 'Halfres' 'b0.1' 'b1' 'b3' 'v10' 'v100' 'v500')

for F in "${Folder[@]}"
do
    echo "Creating the txt: " ${F}

    python ./Assignation.py ${F}

done

