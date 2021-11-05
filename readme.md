# Overview
kb-hashing is a locality sensitive hashing tool for genetic sequence comparison in bioinformatics.

# Compile kb-hashing

Use the following to compile kb-hashing:
```
cd src
g++ -o kbhashing kbhashing.cpp -std=c++11
```

# Simulation data generation

We provide a data simulation program and corresponding script to generate simulated data. Use the following to compile it and run the script:

```
cd src
g++ -o seqsim seqsim.cpp -std=c++11
cd ../src/script
./simulation.sh n m
```

Here n means the length of generated sequences and m represents the number of pairs of sequence in each data set. The generated data will be the folder data. There will be six datasets in that folder with different error rate from 0.05 to 0.3.

