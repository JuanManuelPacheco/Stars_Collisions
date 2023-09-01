Folder=('4Myr' 'Renzo' 'Halfres' 'b0.1' 'b1' 'b3' 'v10' 'v100' 'v500')

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
