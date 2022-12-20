**HH->WWbb Run-3 analysis**

Start with:
` bambooRun -m checkZpeak.py:NanoBaseHHWWbb config/hh_wwbb_run3.yml -o output-test `

or use the following if running on slurm (batch):
` bambooRun -m checkZpeak.py:NanoBaseHHWWbb config/hh_wwbb_run3.yml -o output-test --envConfig path-to-your-bamboo-build/bamboo/examples/ingrid.ini --maxFiles=1 --distributed=driver `

For data and background samples, see ` config/hh_wwbb_run3.yml `
