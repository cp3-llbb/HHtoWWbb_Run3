![](https://img.shields.io/github/v/tag/cp3-llbb/HHtoWWbb_Run3)
[![Run test](https://github.com/cp3-llbb/HHtoWWbb_Run3/actions/workflows/python_test.yml/badge.svg)](https://github.com/cp3-llbb/HHtoWWbb_Run3/actions/workflows/python_test.yml)
![](https://img.shields.io/badge/CMS-Run3-blue)

# HH->WWbb Run-3 analysis
----WORK IN PROGRESS----

This repository uses the **bamboo analysis framework**, you can install it via the instructions here: https://bamboo-hep.readthedocs.io/en/latest/install.html#fresh-install

Then clone this repository in the parent directory containing the bamboo installation:

```bash
git clone https://github.com/cp3-llbb/HHtoWWbb_Run3.git && cd HHtoWWbb_Run3
```

Execute these each time you start from a clean shell on lxplus or any other machine with an cvmfs:
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
source (path to your bamboo installation)/bamboovenv/bin/activate
export PYTHONPATH="${PYTHONPATH}:${PWD}/python/"
```

and the followings before submitting jobs to the batch system (HTCondor, Slurm, Dask and Spark are supported):

```bash
voms-proxy-init --voms cms -rfc --valid 192:00 
export X509_USER_PROXY=$(voms-proxy-info -path)
```
if you encounter problems with accessing files when using batch, the following lines may solve your problem

```bash
voms-proxy-init --voms cms -rfc --valid 192:00  --out ~/private/gridproxy/x509
export X509_USER_PROXY=$HOME/private/gridproxy/x509
```

Then plot various control regions via the following command line using batch (you can pass `--maxFiles 1` to use only 1 file from each sample for a quick test):

```bash
bambooRun -m python/controlPlotter.py config/analysis_2022.yml -o ./outputDir/ --distributed driver --envConfig config/cern.ini --eras combined -c <DL or SL>
```
Instead of passing everytime `--envConfig config/cern.ini`, you can copy the content of that file to `~/.config/bamboorc`.

Pass `--mvaSkim` to produce skims for MVA. as well.

then to produce plots, just execute:

```sh
./scripts/plot_<SL or DL>.sh <path to output/plots>
```

using the `parquet` output file that contains skims and the DNN jupyter notebook, you can perform machine learning applications. Then passing `--mvaEval` option, you can apply DNN score cuts on your analysis.
