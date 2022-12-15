**Bamboo**
Start with:
<!-- bambooRun -m hh_wwbb_run3.py:NanoBaseHHWWbb config/hh_wwbb_run3.yml -o output-test -->

**DATA**
<!-- dasgoclient -query "dataset dataset=/MuonEG/\*2022\*/NANOAOD" -->
- /store/data/Run2022C/MuonEG/NANOAOD/PromptNanoAODv10-v1/40000/55dfef50-3397-4c6a-a3d8-8d3d6ec5feb5.root
- /store/data/Run2022C/SingleMuon/NANOAOD/PromptNanoAODv10_v1-v1/50000/488e74c8-d09f-4724-ae30-99ab7c516d93.root
- /store/data/Run2022C/EGamma/NANOAOD/PromptNanoAODv10-v2/50000/41d53e82-7df7-4ddd-b4e0-45282baefc1d.root

**MC**
**TTbar**
<!-- dataset dataset=/*TTbarToDilepton*/*mcRun3*2022*/NANOAOD* -->
/store/relval/CMSSW\_12\_4\_0/RelValTTbarToDilepton\_14TeV/NANOAODSIM/124X\_mcRun3\_2022\_realistic\_v5-v1/2580000/7c06dffd-839c-411d-9630-4ba84a197b42.root

**DY**
<!-- dataset dataset=/*DYJets*14TeV*/*mcRun3*/NANOAOD* -->
/store/mc/RunIIAutumn18NanoAODv7/DYJetsToLL\_M-50\_HT-1200to2500\_TuneCP5\_14TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020\_102X\_upgrade2018\_realistic\_v21-v1/100000/238D3F8E-234B-AC41-A221-F943653E40DC.root
changed to **/store/mc/Run3Summer19NanoAOD/DYJets\_incl\_MLL-50\_TuneCP5\_14TeV-madgraphMLM-pythia8/NANOAODSIM/2021Scenario\_106X\_mcRun3\_2021\_realistic_v3-v1/50000/C4465203-C5DB-B848-A601-CB396449051A.root**

<!-- xrdcp root://cms-xrd-global.cern.ch////store....... -->
