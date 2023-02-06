----WORK IN PROGRESS----

**HH->WWbb Run-3 analysis**

Execute followings before submitting to the batch system:

`$ mkdir -p ~/private/gridproxy`

`$ voms-proxy-init --voms cms -rfc --valid 96:00  --out ~/private/gridproxy/x509`

`$ export X509_USER_PROXY=$HOME/private/gridproxy/x509`

Use the following command line (replace X):

`$ bambooRun -m checkZpeak.py:NanoBaseHHWWbb config/analysis_2022.yml -o output-Run_2022 --envConfig path-to-your-bamboo-build/bamboo/examples/ingrid.ini --distributed=driver `

For data and background samples, see `config` folder
