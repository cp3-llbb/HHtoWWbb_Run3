----WORK IN PROGRESS----

**HH->WWbb Run-3 analysis**

Install bamboo using the instructions here: https://bamboo-hep.readthedocs.io/en/latest/install.html

Then clone this repository inside the directory containing the bamboo installation:

`git clone https://github.com/Oguz-Guzel/HHWWbb.git`
`cd HHWWbb`

Execute followings before submitting to the batch system:

`$ mkdir -p ~/private/gridproxy`

`$ voms-proxy-init --voms cms -rfc --valid 96:00  --out ~/private/gridproxy/x509`

`$ export X509_USER_PROXY=$HOME/private/gridproxy/x509`

Then start the analysis using the following command line:

`$ bambooRun -m checkZpeak.py:NanoBaseHHWWbb config/analysis_2022.yml -o output-Run_2022 --envConfig ../examples/ingrid.ini --distributed=driver `

For data and background samples, see `config` folder

Check ecal endcap issue in the era E-F-G data only:

`bambooRun -m checkZpeak.py:checkEE config/data_2022_erasEFG.yml -o test --envConfig ../bamboo/examples/cern.ini --distributed=driver`
