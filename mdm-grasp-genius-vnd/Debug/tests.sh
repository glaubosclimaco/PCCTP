#! /bin/bash

# Compilando FPMax
# cd FPMax
# ls
# rm *.o
# rm fim_closed
# rm fim_maximal
# rm fpmax_hnmp
# make
# ls
# cd .. 

ls
rm *.o
rm *.d
rm grasp-genius-vnd
make

# ./grasp-genius-vnd data/g_att48_VT16_T16_W16_75 3 2
# ./grasp-genius-vnd data/g_att48_VT16_T16_W16_75 1 4

# ./grasp-genius-vnd data/${INSTANCE[$i]} $seed $sup

# # Parametros
num_exec=10
num_instances=1
stable_percentage=5


INSTANCE[0]="g_att48_VT9_T10_W29_25"

# INSTANCE[0]="g_att48_VT9_T10_W29_25"
# INSTANCE[1]="g_tsp225_VT75_T75_W75_25"
# INSTANCE[2]="g_pr299_VT59_T60_W180_50"
# INSTANCE[3]="g_pr226_VT76_T75_W75_25"
# INSTANCE[4]="g_lin318_VT63_T64_W191_25"
# INSTANCE[5]="g_gr202_VT39_T122_W41_25"
# INSTANCE[6]="g_gil262_VT88_T87_W87_50"
# INSTANCE[7]="g_a280_VT94_T93_W93_25"
# INSTANCE[8]="g_si175_VT59_T58_W58_75"
# INSTANCE[9]="g_si175_VT35_T105_W35_50"
# INSTANCE[10]="g_rat195_VT39_T117_W39_75"
# INSTANCE[11]="g_tsp225_VT75_T75_W75_50"
# INSTANCE[12]="g_pr107_VT20_T22_W65_50"
# INSTANCE[13]="g_att48_VT16_T16_W16_75"


# echo "INICIO - tests.sh"

for ((sup=2;sup<=5;sup++));
do
	mkdir -p resultados/sup_$sup
	for ((i=0;i<$num_instances;i++));
	do
		for ((seed=1;seed<=$num_exec;seed++));
		do
			echo "executando data/${INSTANCE[$i]} "
			echo "./grasp-genius-vnd data/${INSTANCE[$i]} $seed $sup $stable_percentage"
		    ./grasp-genius-vnd data/${INSTANCE[$i]} $seed $sup $stable_percentage > resultados/sup_$sup/saida_${INSTANCE[$i]}_$seed.txt

		    echo "fim executacao data/${INSTANCE[$i]}"
		    echo "\n"
		done;
	done;
done;

echo "FIM - tests.sh"
