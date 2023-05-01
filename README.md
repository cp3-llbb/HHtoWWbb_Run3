![](https://img.shields.io/github/v/tag/Oguz-Guzel/HHWWbbRun3)
[![Run test](https://github.com/cp3-llbb/HHtoWWbb_Run3/actions/workflows/python_test.yml/badge.svg)](https://github.com/cp3-llbb/HHtoWWbb_Run3/actions/workflows/python_test.yml)
[![CodeQL](https://github.com/Oguz-Guzel/HHWWbbRun3/actions/workflows/github-code-scanning/codeql/badge.svg)](https://github.com/Oguz-Guzel/HHWWbbRun3/actions/workflows/github-code-scanning/codeql)
![](https://img.shields.io/badge/CMS-Run3-blue)

----WORK IN PROGRESS----

**HH->WWbb Run-3 analysis**

Install **bamboo analysis framework** with the instructions here: https://bamboo-hep.readthedocs.io/en/latest/install.html#fresh-install

Then clone this repository in the parent directory containing the bamboo installation:

```bash
git clone https://github.com/cp3-llbb/HHtoWWbb_Run3.git && cd HHtoWWbb_Run3
```

Execute these each time you start from a clean shell:
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
source (path to your bamboo installation)/bamboovenv/bin/activate
export PYTHONPATH="${PYTHONPATH}:${PWD}/python/"
```

and the followings before submitting to the batch system:

```bash
voms-proxy-init --voms cms -rfc --valid 192:00 
export X509_USER_PROXY=$(voms-proxy-info -path)
```
if you encounter problems with accessing files when using batch, the following lines may solve your problem

```bash
voms-proxy-init --voms cms -rfc --valid 192:00  --out ~/private/gridproxy/x509
export X509_USER_PROXY=$HOME/private/gridproxy/x509
```

Then plot various control regions via the following command line using batch (you can pass `--maxFiles 1` to use only 1 file from each sample for a test):

```bash
bambooRun -m python/controlPlotter.py config/analysis_2022.yml -o ./outputDir/ --distributed driver --envConfig config/cern.ini
```
Instead of passing everytime `--envConfig config/cern.ini`, you can copy the content of that file to `~/.config/bamboorc`.
