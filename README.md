[![Run test](https://github.com/Oguz-Guzel/HHWWbbRun3/actions/workflows/python_test.yml/badge.svg)](https://github.com/Oguz-Guzel/HHWWbbRun3/actions/workflows/python_test.yml)
[![CodeQL](https://github.com/Oguz-Guzel/HHWWbbRun3/actions/workflows/github-code-scanning/codeql/badge.svg)](https://github.com/Oguz-Guzel/HHWWbbRun3/actions/workflows/github-code-scanning/codeql)

----WORK IN PROGRESS----

**HH->WWbb Run-3 analysis**

Install **bamboo analysis framework** with the instructions here: https://bamboo-hep.readthedocs.io/en/latest/install.html

Then clone this repository in the parent directory containing the bamboo installation:

```
git clone https://github.com/Oguz-Guzel/HHWWbbRun3.git && cd HHWWbbRun3
```

Execute these each time you start from a clean shell:
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
source (path to your bamboo installation)/bamboovenv/bin/activate
export PYTHONPATH="${PYTHONPATH}:${PWD}/python/"
```

and the followings before submitting to the batch system:

```
voms-proxy-init --voms cms -rfc --valid 192:00 
export X509_USER_PROXY=$(voms-proxy-info -path)
```

Then plot various control regions via the following command line:

```
bambooRun -m python/controlPlotter.py:NanoBaseHHWWbb config/analysis_2022.yml -o output-Run_2022 --envConfig ../bamboo/examples/ingrid.ini --distributed=driver
```
