# ZCURVE 2026
Z-curve-based Prokaryotic Protein-Coding Gene Recognition System
## Setup
```bash
g++ -DZLIB -o zcurve2026.exe -fopenmp -mavx -static -mfma Main.cpp BioIO.cpp BioUtil.cpp Encoding.cpp Model.cpp svm.cpp -lz -O2
```