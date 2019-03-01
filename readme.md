## SFB 925 QMC workshop 
Hamburg, 08.05.2015


The example programs were based on the examples provided by the creators of the TRIQS based applications cthyb and SOM. More information \
https://github.com/TRIQS/triqs
https://triqs.ipht.cnrs.fr/1.4/applications/cthyb/ \
http://krivenko.github.io/som/ 

### Usage
0. Clone this repository on your machine `git clone https://github.com/patryk-kubiczek/qmc_workshop.git`
1. Addin pytriqs to PATH: edit load_triqs.sh and then `source load_triqs.sh`
2. Heavy calculations should be run in parallel `mpirun -np 4 pytriqs *.py`
3. At the end, plot the results with `pytriqs plot.py`
