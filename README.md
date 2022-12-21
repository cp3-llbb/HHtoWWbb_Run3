**HH->WWbb Run-3 analysis**

Start with:
` bambooRun -m checkZpeak.py:NanoBaseHHWWbb config/hh_wwbb_run3.yml -o output-test `

or use the following if running on slurm (batch):
` bambooRun -m checkZpeak.py:NanoBaseHHWWbb config/hh_wwbb_run3.yml -o output-test --envConfig path-to-your-bamboo-build/bamboo/examples/ingrid.ini --distributed=driver `

For data and background samples, see ` config/hh_wwbb_run3.yml `

Execute followings before submitting to a batch system

`mkdir -p ~/private/gridproxy`

`$ voms-proxy-init --voms cms -rfc --valid 96:00  --out ~/private/gridproxy/x509`

`$ export X509_USER_PROXY=$HOME/private/gridproxy/x509`
