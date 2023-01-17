**HH->WWbb Run-3 analysis**

Start with:
` bambooRun -m checkZpeak.py:NanoBaseHHWWbb config/hh_wwbb_run3.yml -o output-test `

or use the following if running on slurm (batch) (replace <X>):
` bambooRun -m checkZpeak.py:NanoBaseHHWWbb config/datasets_2022<X>.yml -o output-Run_2022<X> --envConfig path-to-your-bamboo-build/bamboo/examples/ingrid.ini --distributed=driver --era 2022<X> `

For data and background samples, see `config` folder

Execute followings before submitting to a batch system

`$ mkdir -p ~/private/gridproxy`

`$ voms-proxy-init --voms cms -rfc --valid 96:00  --out ~/private/gridproxy/x509`

`$ export X509_USER_PROXY=$HOME/private/gridproxy/x509`
