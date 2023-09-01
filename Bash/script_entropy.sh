Folder=('4Myr' 'Renzo' 'Halfres' 'b0.1' 'b1' 'b3' 'v10' 'v100' 'v500')
fin=('1221' '1606' '1464' '1269' '1090' '1342' '1263' '0246' '1247')

sum=0
i=1

for F in "${Folder[@]}"
do 
    echo "Creating the txt: " ${F}
    echo "final out: " ${fin[sum]} 
    python ./Entropy.py ${F} ${fin[sum]}
    sum=$(($sum + $i))
done

