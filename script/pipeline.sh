for r in 0.01 0.1 0.15 0.2 0.25 0.3
do
	../src/kbhashing "$1" "$2" "$3" "$4" "$5" < ../data/simdata_"$6"_"$r"_"$7"
done