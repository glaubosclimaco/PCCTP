linha=7
for file in $(ls -r *); do
informacao=$(sed -n $linha'p' $file)
echo $informacao >> saida.txt
done;
