
from bamboo.analysismodules import NanoAODModule, HistogramsModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.treedecorators import NanoAODDescription
from bamboo.plots import Plot, CutFlowReport
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

from itertools import chain


class NanoBaseHHWWbb(NanoAODModule, HistogramsModule):
    def __init__(self, args):
        super(NanoBaseHHWWbb, self).__init__(args)
        self.plotDefaults = {"show-ratio": True,
                             "y-axis-show-zero": True,
                             "normalized": False,
                             "y-axis": "Events",
                             "log-y": "both",
                             "ratio-y-axis-range": [0.8, 1.2],
                             "ratio-y-axis": '#frac{Data}{MC}',
                             "sort-by-yields": True}

    def addArgs(self, parser):
        super(NanoBaseHHWWbb, self).addArgs(parser)
        parser.add_argument("--era",
                            action='store',
                            type=str,
                            default=None,
                            help='This has no use right now!')

    def prepareTree(self, tree, sample=None, sampleCfg=None):

        def isMC():
            if sampleCfg['type'] == 'data':
                return False
            elif sampleCfg['type'] == 'mc':
                return True
            else:
                print(
                    f"Please specify the type of {sample} dataset in the configuration file (data or mc) and re-run.")
                exit()

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

        def getNanoAODDescription():  # implemented from Sebastien's analysis (mentioned on issue #101 on bamboo gitlab page)
            groups = ["HLT_", "MET_", "Pileup_"]
            collections = ["nElectron", "nJet", "nMuon", "nFatJet"]
            varReaders = []
            return NanoAODDescription(groups=groups, collections=collections, systVariations=varReaders)

        tree, noSel, backend, lumiArgs = super(NanoBaseHHWWbb, self).prepareTree(tree=tree,
                                                                                 sample=sample,
                                                                                 sampleCfg=sampleCfg,
                                                                                 description=getNanoAODDescription(),
                                                                                 backend="lazy")
        ### Triggers ###
        # EGamma
        addHLTPath('EGamma', 'Ele32_WPTight_Gsf')
        addHLTPath('EGamma', 'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL')


        return tree, noSel, backend, lumiArgs

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []
        yields = CutFlowReport("yields", printInLog=True, recursive=True)
        plots.append(yields)
        yields.add(noSel, 'No Selection')
        #############################################################################
        #                                 Electrons                                 #
        #############################################################################
        electrons = op.sort(op.select(tree.Electron, lambda el: op.AND(
            el.pt >= 7.,
            op.abs(el.eta) <= 2.5,
            op.abs(el.dxy) <= 0.05,
            op.abs(el.dz) <= 1.,
            el.miniPFRelIso_all <= 0.4,
            el.sip3d <= 8,
            # el.mvaNoIso_WPL,
            el.lostHits <= 1
        )), lambda el: -el.pt)

        #############################################################################
        #                          Gen Weight and Triggers                          #
        #############################################################################
        if self.is_MC:
            noSel = noSel.refine('genWeight', weight=tree.genWeight, cut=(
                op.OR(*chain.from_iterable(self.triggersPerPrimaryDataset.values()))))
        else:
            noSel = noSel.refine('trigger', cut=[makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggersPerPrimaryDataset)])
        #############################################################################
        #                               Selections                                  #
        #############################################################################

        # has at least one electron pair
        hasElEl = noSel.refine("hasOSElEl", cut=[op.rng_len(electrons) >= 2,
                                                 electrons[0].charge != electrons[1].charge, electrons[0].pt > 20., electrons[1].pt > 10.])

        #############################################################################
        #                                 Plots                                     #
        #############################################################################
        plots.extend([
            Plot.make2D("el_eta_vs_phi", [electrons[0].eta, electrons[0].phi], hasElEl, [
                        EqBin(50, -2.5, 2.5), EqBin(60, -6.29, 6.29)], title="electron eta vs phi")
        ])

        return plots
