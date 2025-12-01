# ZCURVE 2026
A Prokaryotic Protein-coding Gene Recognition System
## Setup
```bash
g++ -DZLIB -o zcurve -fopenmp -mavx -mfma Main.cpp BioIO.cpp BioUtil.cpp Encoding.cpp Model.cpp svm.cpp -lz -O2
```
If you don't want to use zlib, just remove the -DZLIB option.