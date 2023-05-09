from bamboo.analysismodules import NanoAODModule, HistogramsModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection, configureJets, configureType1MET
from bamboo import treefunctions as op
from bamboo import treedecorators as td

from itertools import chain


class NanoBaseHHWWbb(NanoAODModule, HistogramsModule):
    def __init__(self, args):
        super(NanoBaseHHWWbb, self).__init__(args)

    def addArgs(self, parser):
        super(NanoBaseHHWWbb, self).addArgs(parser)
        parser.add_argument("--era",
                            action='store',
                            type=str,
                            default=None,
                            help='This has no use right now!')

    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        def isMC():
            if sampleCfg['type'] == 'data':
                return False
            elif sampleCfg['type'] == 'mc':
                return True
            else:
                raise RuntimeError(
                    f"The type '{sampleCfg['type']}' of {sample} dataset not understood.")

        era = sampleCfg['era']  # reserved for future use
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
            collections = ["nElectron", "nJet", "nMuon", "nFatJet", "nSubJet"]
            mcCollections = ["nGenDressedLepton",
                             "nGenJet", "nGenPart", "nCorrT1METJet"]
            varReaders = []
            if isMC:
                varReaders.append(td.CalcCollectionsGroups(Jet=("pt", "mass")))
                varReaders.append(
                    td.CalcCollectionsGroups(GenJet=("pt", "mass")))
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

        if self.is_MC:
            configureJets(
                variProxy=tree._Jet,
                jetType="AK4PFchs",
                jsonFile="/home/ucl/cp3/aguzel/Bamboo/CMSJMECalculators/tests/data/2018_UL_jet_jerc.json.gz",
                jec="Summer19UL18_V5_MC",
                #   smear="Autumn18_V7b_MC",
                jesUncertaintySources=["Total"],
                enableSystematics=lambda v: not "jesTotal" in v,
                isMC=self.is_MC,
                backend=backend,
                uName=sample
            )
            configureType1MET(
                variProxy=tree._MET,
                jsonFile="/home/ucl/cp3/aguzel/Bamboo/CMSJMECalculators/tests/data/2018_UL_jet_jerc.json",
                smear="Summer19UL18_JRV2_MC",
                jec="Summer19UL18_V5_MC",
                # jesUncertainties=["Total"],
                splitJER=True,
                addHEM2018Issue=False,
                isMC=self.is_MC,
                backend=backend,
                uName=sample,
            )

        return tree, noSel, backend, lumiArgs
