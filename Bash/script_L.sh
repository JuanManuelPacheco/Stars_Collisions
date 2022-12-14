Folder=('coll_4Myr' 'coll_Renzo_refosco' 'coll_b0.1' 'coll_b1' 'coll_noangm_v100' 'coll_noangm_v500')

for F in "${Folder[@]}"
do
    cd ${F}/
    for output in $(ls out*.ascii | xargs -n 1 basename)
    do
	name=$(echo ${output} | cut -f 1 -d '.')

	if [ -f "./Vel_${name}.png" ]
	then
            echo "Existing png"
	else
            echo "Creating the png: " ${name}
            python ../L_evo.py ${F} ${output} ${name}
        fi
    done
    cd ..
done

