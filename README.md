----WORK IN PROGRESS----

**HH->WWbb Run-3 analysis**

Install **bamboo analysis framework** with the instructions here: https://bamboo-hep.readthedocs.io/en/latest/install.html

Then clone this repository inside the directory containing the bamboo installation:

```
git clone https://github.com/Oguz-Guzel/HHWWbb.git && cd HHWWbb
```

Execute followings before submitting to the batch system:

```
mkdir -p ~/private/gridproxy
voms-proxy-init --voms cms -rfc --valid 96:00  --out ~/private/gridproxy/x509
export X509_USER_PROXY=$HOME/private/gridproxy/x509
```
Then check various control plots via the following command line:

```
bambooRun -m analysis.py:NanoBaseHHWWbb config/analysis_2022.yml -o output-Run_2022 --envConfig ../bamboo/examples/ingrid.ini --distributed=driver
```

versioning

**v1**: first succesful run on all eras

**v2**: add triggers on MC

**v3**: add Wjets and VV MC

**v4**: add DY-M-10to50

**v5**: add SL and DL channels

