# CLAS22_gev_ana

Quick test to show dst2root is working and output a few histograms to a root file.

Prerequisites:
* [cern root](https://root.cern.ch/)

To build:

```bash
git clone https://github.com/kneupane1/CLAS22_gev_ana
mkdir -p clas22_analysis/build
cd clas22_analysis/build
cmake ..
make
```
To run:
1. To produce csv file in order to plot in jupyter notebook or python
```bash
BEAM_E=22.0 NUM_THREADS=4 ./clas22_yields yields.csv /path/to/input/*.root
```
2. To produce root file
```bash
BEAM_E=22.0 NUM_THREADS=4 ./clas22_analysis clas22_output_file.root /path/to/input/*.root
```
If it breaks, reduce the number of threads for the number of files. In general each thread should have 2 or more files and the number of threads should be less than or equal to the number of cores you are using.

So for 16 files on a 4 core computer use NUM_THREADS=4 and each thread will process 4 files. For 4 files on a 4 core computer use NUM_THREADS=1 or NUM_THREADS=2 so each thread will have 4 or 2 files respecfully.
