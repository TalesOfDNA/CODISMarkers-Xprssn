SAMPLE_SIZE=200
# 40 60 80 100 120 140 160 180 200)
NUM_REPS=1000

for i in $(seq $NUM_REPS); do
	sbatch RunSimis_BySTR.sh $SAMPLE_SIZE $i 
done

