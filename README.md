----WORK IN PROGRESS----

**HH->WWbb Run-3 analysis**

Install **bamboo analysis framework** with the instructions here: https://bamboo-hep.readthedocs.io/en/latest/install.html

Then clone this repository inside the directory containing the bamboo installation:

`git clone https://github.com/Oguz-Guzel/HHWWbb.git && cd HHWWbb`

Execute followings before submitting to the batch system:

`$ mkdir -p ~/private/gridproxy`

`$ voms-proxy-init --voms cms -rfc --valid 96:00`

Then check the Z peak and various control plots via the following command line:

`$ bambooRun -m checkZpeak.py:NanoBaseHHWWbb config/analysis_2022.yml -o output-Run_2022 --envConfig ../bamboo/examples/ingrid.ini --distributed=driver `
