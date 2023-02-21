----WORK IN PROGRESS----

**HH->WWbb Run-3 analysis**

Install bamboo with the instructions here: https://bamboo-hep.readthedocs.io/en/latest/install.html

Then clone this repository inside the directory containing the bamboo installation:

`git clone https://github.com/Oguz-Guzel/HHWWbb.git && cd HHWWbb`

Execute followings before submitting to the batch system:

`$ mkdir -p ~/private/gridproxy`

`$ voms-proxy-init --voms cms -rfc --valid 96:00  --out ~/private/gridproxy/x509`

`$ export X509_USER_PROXY=$HOME/private/gridproxy/x509`

Then check the Z peak and various control plots via the following command line:

`$ bambooRun -m checkZpeak.py:NanoBaseHHWWbb config/analysis_2022.yml -o output-Run_2022 --envConfig ../examples/ingrid.ini --distributed=driver `

You can check also the Ecal Endcap issue (a pipe leak) happened during the era E-F-G data only or before (eras C and D):

`bambooRun -m checkZpeak.py:checkEE config/data_2022_splitErasCDandEFG.yml -o test --envConfig ../bamboo/examples/cern.ini --distributed=driver --era 2022<CD or EFG>`
