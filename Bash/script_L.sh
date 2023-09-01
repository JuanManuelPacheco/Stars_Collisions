Folder=('4Myr' 'Renzo' 'Halfres' 'b0.1' 'b1' 'b3' 'v10' 'v100' 'v500')

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

