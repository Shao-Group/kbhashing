for r in 0.01 0.1 0.15 0.2 0.25 0.3
do
	../src/seqsim "$1" "$r" "$2" > ../data/simdata_"$1"_"$r"_"$2"
done