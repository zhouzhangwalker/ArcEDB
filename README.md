# ArcEDB:  An Arbitrary-Precision Encrypted Database via (Amortized) Modular Homomorphic Encryption
This code is the implmentation of the paper ArcEDB:  An Arbitrary-Precision Encrypted Database via (Amortized) Modular Homomorphic Encryption.
## Requirements

```
git 
gcc >= 10
cmake >= 3.16
GMP 6.2.0
OpenMP
```

## Building ArcEDB
You can build the ArcEDB (out-of-source) for your machine by executing the following commands:
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DSEAL_THROW_ON_TRANSPARENT_CIPHERTEXT=OFF -DSEAL_USE_INTEL_HEXL=ON
make
```

Following the instructions above, output binaries will be generated in the `build/bin/` directory. You can run these binaries by:
```
$./bin/comparison_test
$./bin/simd_comparison_test
$./bin/relational_queries_test
$./bin/time_series_queries_test
```

If you have Docker on your system, this will do above on docker.
```
docker build -t arcedb .
docker run -i -t arcedb
```
After the build completes, the following examples can be found in 
the docker container.

## Examples
### Homomorphic Comparison
- codes `test/comparison_test.cpp`  
- output binary `build/bin/comparison_test` 
- This demo shows the homomorphic comparison operator `ArbHCMP` in ArcEDB. 

### Homomorphic SIMD Comparison
- codes `test/simd_comparison_test.cpp`  
- output binary `build/bin/simd_comparison_test` 
- This demo shows the homomorphic comparison operator `SIMDArbHCMP` in ArcEDB.

### Relational Query Evaluation
- codes `test/relational_queries_test.cpp`
- output binary `build/bin/relational_queries_test`
- This demo shows the evaluation of relational queries over a 32K rows of encrypted database.

### Time-Series Query Evaluation
- codes `test/time_series_queries_test.cpp`
- output binary `build/bin/time_series_queries_test`
- This demo shows the evaluation of time-series queries over a 32K rows of encrypted database.
