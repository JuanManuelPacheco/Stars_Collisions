Folder='coll_b0.1'
for output in $(ls ${Folder}/*.sph | xargs -n 1 basename)
do
echo "Creating the ascii file: " ${output}

#name=$(echo ${output} | cut -f 1 -d '.')
#python Setup/velocity_curves.py ${simulation}/ ${name}

splash -f starsmasher to ascii ${output}
done
