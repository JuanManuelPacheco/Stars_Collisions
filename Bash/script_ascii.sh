Folder=('coll_Renzo_refosco' 'coll_b0.1' 'coll_b1' 'coll_noangm_v100' 'coll_noangm_v500')

for F in "${Folder[@]}"
do
    cd ${F}/
    for output in $(ls out*.sph | xargs -n 1 basename)
    do
        echo "Creating the ascii file: " ${output}

        splash -f starsmasher to ascii ${output}
    done
    cd ..
done
