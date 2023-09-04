from bamboo.analysismodules import NanoAODModule, HistogramsModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection, configureJets, configureType1MET
from bamboo import treefunctions as op
from bamboo import treedecorators as td

from itertools import chain
import os

class NanoBaseHHWWbb(NanoAODModule, HistogramsModule):
    def __init__(self, args):
        super(NanoBaseHHWWbb, self).__init__(args)

    def addArgs(self, parser):
        super(NanoBaseHHWWbb, self).addArgs(parser)
        parser.add_argument("-c", "--channel",
                            dest="channel",
                            type=str,
                            required=True,
                            help='Channel to be selected between SL and DL')
        parser.add_argument("--mvaSkim",
                            dest="mvaSkim",
                            action="store_true",
                            help="Produce skims for MVA")
        parser.add_argument("--mva",
                            dest="mva",
                            action="store_true",
                            help="Run MVA")

    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        def isMC():
            if sampleCfg['type'] == 'data':
                return False
            elif sampleCfg['type'] in ['mc', 'signal']:
                return True
            else:
                raise RuntimeError(
                    f"The type '{sampleCfg['type']}' of {sample} dataset not understood.")

        era = sampleCfg['era']
        self.is_MC = isMC()
        self.triggersPerPrimaryDataset = {}

        def addHLTPath(PD, HLT):
            if PD not in self.triggersPerPrimaryDataset.keys():
                self.triggersPerPrimaryDataset[PD] = []
            try:
                self.triggersPerPrimaryDataset[PD].append(
                    getattr(tree.HLT, HLT))
            except AttributeError:
                print("Couldn't find branch tree.HLT.%s, will omit it!" % HLT)

        def getNanoAODDescription():
            groups = ["HLT_", "MET_", "RawMET_"]
            collections = ["nElectron", "nJet",
                           "nMuon", "nFatJet", "nSubJet", "nTau"]
            mcCollections = ["nGenDressedLepton", "nGenJet", "nGenPart"]
            varReaders = []
            if isMC:
                varReaders.append(td.CalcCollectionsGroups(Jet=("pt", "mass")))
                varReaders.append(td.CalcCollectionsGroups(GenJet=("pt", "mass")))
                varReaders.append(td.CalcCollectionsGroups(MET=("pt", "phi")))
                return td.NanoAODDescription(groups=groups, collections=collections + mcCollections, systVariations=varReaders)
            else:
                return td.NanoAODDescription(groups=groups, collections=collections, systVariations=varReaders)

        tree, noSel, backend, lumiArgs = super(NanoBaseHHWWbb, self).prepareTree(tree=tree,
                                                                                 sample=sample,
                                                                                 sampleCfg=sampleCfg,
                                                                                 description=getNanoAODDescription(),
                                                                                 backend=backend)
        ### Triggers ###
        # Muon
        addHLTPath('Muon', 'IsoMu24')
        addHLTPath('Muon', 'IsoMu27')
        addHLTPath('Muon', 'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8')
        # EGamma
        addHLTPath('EGamma', 'Ele32_WPTight_Gsf')
        addHLTPath('EGamma', 'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL')
        # MuonEG
        addHLTPath('MuonEG', 'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ')
        # SingleMuon
        addHLTPath('SingleMuon', 'IsoMu24')
        addHLTPath('SingleMuon', 'IsoMu27')
        # DoubleMuon
        addHLTPath('DoubleMuon', 'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8')

        # Gen Weight and Triggers
        if self.is_MC:
            noSel = noSel.refine('genWeight', weight=tree.genWeight, cut=(
                op.OR(*chain.from_iterable(self.triggersPerPrimaryDataset.values()))))
        else:
            noSel = noSel.refine('trigger', cut=[makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggersPerPrimaryDataset)])

        # if sampleCfg['type'] == 'mc':
        #     JECTagDatabase = {"2022": "Winter22Run3_V2_MC",
        #                       "2022EE": "Summer22EEPrompt22_V1_MC"}
        #     if era in JECTagDatabase.keys():
        #         print(f"Configure jet corrections for era {era}.")
        #         configureJets(
        #             variProxy               = tree._Jet,
        #             jetType                 = "AK4PFPuppi",
        #             jec                     = JECTagDatabase[era],
        #             #  smear                   = JERTagDatabase['2022EE'],
        #             jecLevels               = "default",
        #             jesUncertaintySources   = "All",
        #             mayWriteCache           = self.args.distributed != "worker",
        #             isMC                    = self.is_MC,
        #             backend                 = backend,
        #             uName                   = sample
        #             )
        #         print(f"Finished configuring jet corrections for era {era}!")
        # if sampleCfg['type'] == 'data':
        #     JECTagDatabase = {"2022C": "Winter22Run3_RunC_V2_DATA",
        #                       "2022D": "Winter22Run3_RunD_V2_DATA",
        #                       "2022F": "Summer22EEPrompt22_RunF_V1_DATA",
        #                       "2022G": "Summer22EEPrompt22_RunG_V1_DATA"}
        #     if era in JECTagDatabase.keys():
        #         print(f"Configure jet corrections...")
        #         configureJets(
        #             variProxy               = tree._Jet,
        #             jetType                 = "AK4PFPuppi",
        #             jec                     = JECTagDatabase[era],
        #             jecLevels               = "default",
        #             jesUncertaintySources   = "All",
        #             mayWriteCache           = self.args.distributed != "worker",
        #             isMC                    = self.is_MC,
        #             backend                 = backend,
        #             uName                   = sample
        #             )
        #         print(f"Finished configuring jet corrections for era {era}!")
        return tree, noSel, backend, lumiArgs
        
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Postprocess: run plotIt

        The list of plots is created if needed (from a representative file,
        this enables rerunning the postprocessing step on the results files),
        and then plotIt is executed
        """
        if not self.plotList:
            self.plotList = self.getPlotList(resultsdir=resultsdir, config=config)
        from bamboo.plots import Plot, DerivedPlot, CutFlowReport
        plotList_cutflowreport = [ap for ap in self.plotList if isinstance(ap, CutFlowReport)]
        plotList_plotIt = [ap for ap in self.plotList
                           if (isinstance(ap, Plot) or isinstance(ap, DerivedPlot))
                           and len(ap.binnings) == 1]
        eraMode, eras = self.args.eras
        if eras is None:
            eras = list(config["eras"].keys())
        if plotList_cutflowreport:
            from bamboo.analysisutils import printCutFlowReports
            printCutFlowReports(
                config, plotList_cutflowreport, workdir=workdir, resultsdir=resultsdir,
                readCounters=self.readCounters, eras=(eraMode, eras), verbose=self.args.verbose)
        if plotList_plotIt:
            from bamboo.analysisutils import writePlotIt, runPlotIt
            import os
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(
                config, plotList_plotIt, cfgName, eras=eras, workdir=workdir, resultsdir=resultsdir,
                readCounters=self.readCounters, plotDefaults=self.plotDefaults,
                vetoFileAttributes=self.__class__.CustomSampleAttributes)
            runPlotIt(
                cfgName, workdir=workdir, plotIt=self.args.plotIt, eras=(eraMode, eras),
                verbose=self.args.verbose)
        
        from bamboo.plots import Skim
        skims = [ap for ap in self.plotList if isinstance(ap, Skim)]

        from bamboo.analysisutils import loadPlotIt
        p_config, samples, _, systematics, legend = loadPlotIt(config, [], eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes)
        
        if self.args.mvaSkim and skims:
            from bamboo.analysisutils import loadPlotIt
            from bamboo.root import gbl
            import pandas as pd
            import os.path
            
            for skim in skims:
                frames = []
                for smp in samples:
                    for cb in (smp.files if hasattr(smp, "files") else [smp]):
                        tree = cb.tFile.Get(skim.treeName)
                        if not tree:
                            print("WARNING: tree %s not found in file %s" % (skim.treeName, cb.tFile.GetName()))
                            print("         skipping...")
                        else:
                            N = tree.GetEntries()
                            cols = gbl.ROOT.RDataFrame(tree).AsNumpy()
                            cols["weight"] *= cb.scale
                            cols["process"] = [smp.name]*len(cols["weight"])
                            frames.append(pd.DataFrame(cols))
                df = pd.concat(frames)
                df["process"] = pd.Categorical(df["process"], categories=pd.unique(df["process"]))
                pqoutname = os.path.join(resultsdir, f"{skim.name}.parquet")
                df.to_parquet(pqoutname)
                print(f"Saved dataframe for skim {skim.name} to {pqoutname}")
