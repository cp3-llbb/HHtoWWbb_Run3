
from bamboo.analysismodules import NanoAODModule, HistogramsModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection, printCutFlowReports

from bamboo.treedecorators import NanoAODDescription

from bamboo.plots import Plot, CutFlowReport
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op


class NanoBaseHHWWbb(NanoAODModule, HistogramsModule):
    def __init__(self, args):
        super(NanoBaseHHWWbb, self).__init__(args)
        self.plotDefaults = {"show-ratio": True,
                             "y-axis-show-zero": True,
                             # "normalized": True,
                             "y-axis": "Events",
                             "log-y": "both",
                             "ratio-y-axis-range": [0.8, 1.2],
                             "ratio-y-axis": '#frac{Data}{MC}',
                             "sort-by-yields": False}

    def addArgs(self, parser):
        super(NanoBaseHHWWbb, self).addArgs(parser)
        parser.add_argument("--era",
                            action='store',
                            type=int,
                            default=2022,
                            help='It has no use right now')

    def prepareTree(self, tree, sample=None, sampleCfg=None):
        era = sampleCfg['era']
        self.is_MC = self.isMC(sample)
        self.triggersPerPrimaryDataset = {}

        def addHLTPath(PD, HLT):
            if PD not in self.triggersPerPrimaryDataset.keys():
                self.triggersPerPrimaryDataset[PD] = []
            try:
                self.triggersPerPrimaryDataset[PD].append(
                    getattr(tree.HLT, HLT))
            except AttributeError:
                print("Couldn't find branch tree.HLT.%s, will omit it!" % HLT)
        tree, noSel, backend, lumiArgs = super(NanoBaseHHWWbb, self).prepareTree(tree=tree,
                                                                                 sample=sample,
                                                                                 sampleCfg=sampleCfg,
                                                                                 description=NanoAODDescription.get(
                                                                                     tag="v1",
                                                                                     year=era,
                                                                                     isMC=self.is_MC,
                                                                                     systVariations=None),
                                                                                 backend="lazy")
        if era == "2022":
            # MuonEG
            addHLTPath('MuonEG', 'Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL')
            # EGamma
            addHLTPath('EGamma', 'Ele32_WPTight_Gsf')
            addHLTPath('EGamma', 'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL')
            # SingleMuon
            addHLTPath('SingleMuon', 'IsoMu24')
            addHLTPath('SingleMuon', 'IsoMu27')
            # Tau
            # addHLTPath('Tau', '')

        return tree, noSel, backend, lumiArgs

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []
        yields = CutFlowReport("yields")
        plots.append(yields)
        yields.add(noSel, 'No Selection')
        #############################################################################
        #                                 Muons                                     #
        #############################################################################
        muons = op.sort(op.select(tree.Muon, lambda mu: op.AND(
            mu.pt >= 5.,
            op.abs(mu.eta) <= 2.4,
            op.abs(mu.dxy) <= 0.05,
            op.abs(mu.dz) <= 0.1,
            mu.miniPFRelIso_all <= 0.4,
            mu.sip3d <= 8,
            mu.tightId
        )), lambda mu: -mu.pt)
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
        #                                 AK8 Jets                                  #
        #############################################################################
        ak8Jets = op.sort(op.select(tree.FatJet, lambda j: op.AND(
            j.jetId & 2,  # tight
            j.pt > 200.,
            op.abs(j.eta) <= 2.4
        )), lambda jet: -jet.pt)
        ak8BJets = op.select(
            ak8Jets, lambda fatjet: fatjet.btagDeepB > 0.4184)  # 2018 WP

        #############################################################################
        #                                 AK4 Jets                                  #
        #############################################################################
        ak4Jets = op.sort(op.select(tree.Jet, lambda j: op.AND(
            j.jetId & 2,  # tight
            j.pt >= 25.,
            op.abs(j.eta) < 2.4
        )), lambda jet: -jet.pt)
        ak4BJets = op.select(
            ak4Jets, lambda jet: jet.btagDeepB > 0.2770)  # 2018 WP
        #############################################################################
        #                          Gen Weight and Triggers                          #
        #############################################################################
        if self.is_MC:
            noSel = noSel.refine('genWeight', weight=tree.genWeight)
            genWeightSel = noSel.refine('genWeightSel')
        else:
            noSel = noSel.refine('Triggers', cut=[makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggersPerPrimaryDataset)])
            triggerSel = noSel.refine('triggerSel')
        #############################################################################
        #                               Selections                                  #
        #############################################################################

        hasElEl = noSel.refine("hasOSElEl", cut=[op.rng_len(electrons) >= 2,
                                                 electrons[0].charge != electrons[1].charge, electrons[0].pt > 20., electrons[1].pt > 10.])

        hasTwoJetsElEl = hasElEl.refine(
            "hasTwoJetsElEl", cut=[op.rng_len(ak4Jets) >= 2])

        hasMuMu = noSel.refine("hasOSMuMu", cut=[op.rng_len(muons) >= 2,
                                                 muons[0].charge != muons[1].charge, muons[0].pt > 20., muons[1].pt > 10.])

        hasTwoJetsMuMu = hasMuMu.refine(
            'hasTwoJetsMuMu', cut=[op.rng_len(ak4Jets) >= 2])

        #############################################################################
        #                                 Plots                                     #
        #############################################################################

        # plots.extend([
        #     Plot.make1D("nEl_NoSel", op.rng_len(electrons), noSel, EqBin(
        #         10, 0., 10.), xTitle="Number of electrons"),
        #     Plot.make1D("nEl_HasElEl", op.rng_len(electrons), hasElEl, EqBin(
        #         10, 0., 10.), xTitle="Number of electrons"),
        #     Plot.make1D("nEl_hasTwoJetsElEl", op.rng_len(electrons), hasTwoJetsElEl, EqBin(
        #         10, 0., 10.), xTitle="Number of electrons"),
        #     Plot.make1D("nMu_NoSel", op.rng_len(muons), noSel, EqBin(
        #         10, 0., 10.), xTitle="Number of muons"),
        #     Plot.make1D("nMu_HasMuMu", op.rng_len(muons), hasMuMu, EqBin(
        #         10, 0., 10.), xTitle="Number of muons"),
        #     Plot.make1D("nMu_hasTwoJetsMuMu", op.rng_len(muons), hasTwoJetsMuMu, EqBin(
        #         10, 0., 10.), xTitle="Number of muons"),
        #     Plot.make1D("nJet_NoSel", op.rng_len(ak4Jets), noSel, EqBin(
        #         10, 0., 10.), xTitle="Number of jets"),
        #     Plot.make1D("nJet_HasElEl", op.rng_len(ak4Jets), hasElEl, EqBin(
        #         10, 0., 10.), xTitle="Number of jets"),
        #     Plot.make1D("nJet_HasMuMu", op.rng_len(ak4Jets), hasMuMu, EqBin(
        #         10, 0., 10.), xTitle="Number of jets"),
        #     Plot.make1D("nJet_hasTwoJetsElEl", op.rng_len(ak4Jets), hasTwoJetsElEl, EqBin(
        #         10, 0., 10.), xTitle="Number of jets"),
        #     Plot.make1D("nJet_hasTwoJetsMuMu", op.rng_len(ak4Jets), hasTwoJetsMuMu, EqBin(
        #         10, 0., 10.), xTitle="Number of jets"),
        #     Plot.make1D("massZto2e", op.invariant_mass(electrons[0].p4, electrons[1].p4),
        #                 hasElEl, EqBin(120, 40., 120.), title="mass of Z to 2e",
        #                 xTitle="Invariant Mass of Nelectrons=2 (in GeV/c^2)"),
        #     Plot.make1D("massZto2e_hasTwoJets", op.invariant_mass(electrons[0].p4, electrons[1].p4),
        #                 hasTwoJetsElEl, EqBin(120, 40., 120.), title="mass of Z to 2e",
        #                 xTitle="Invariant Mass of Nelectrons=2 (in GeV/c^2)"),
        #     Plot.make1D("massZto2mu", op.invariant_mass(muons[0].p4, muons[1].p4),
        #                 hasMuMu, EqBin(120, 40., 120.), title="mass of Z to 2mu",
        #                 xTitle="Invariant Mass of Nmuons=2 (in GeV/c^2)"),
        #     Plot.make1D("massZto2mu_hasTwoJets", op.invariant_mass(muons[0].p4, muons[1].p4),
        #                 hasTwoJetsMuMu, EqBin(120, 40., 120.), title="mass of Z to 2mu",
        #                 xTitle="Invariant Mass of Nmuons=2 (in GeV/c^2)")
        # ])

        # hasOneJet = noSel.refine('hasOneJet', cut=[op.rng_len(ak4Jets) >= 1])
        # hasTwoJets = noSel.refine('hasTwoJets', cut=[op.rng_len(ak4Jets) >= 2])

        # plots.append(Plot.make1D("leadingJetPt_hasOneJet", ak4Jets[0].pt,
        #                          hasOneJet, EqBin(250, 0., 250.), title="leading jet p_T",
        #                          xTitle="Leading Jet p_T (GeV/c^2)"))
        # plots.append(Plot.make1D("leadingJetPt_hasTwoJets", ak4Jets[0].pt,
        #                          hasTwoJets, EqBin(250, 0., 250.), title="leading jet p_T",
        #                          xTitle="Leading Jet p_T (GeV/c^2)"))
        # plots.append(Plot.make1D("leadingJetPt_hasTwoJetsElEl", ak4Jets[0].pt,
        #                          hasTwoJetsElEl, EqBin(250, 0., 250.), title="leading jet p_T",
        #                          xTitle="Leading Jet p_T (GeV/c^2)"))
        # plots.append(Plot.make1D("leadingJetPt_hasTwoJetsMuMu", ak4Jets[0].pt,
        #                          hasTwoJetsMuMu, EqBin(250, 0., 250.), title="leading jet p_T",
        #                          xTitle="Leading Jet p_T (GeV/c^2)"))

        # plots.append(Plot.make1D("subleadingJetPt_NoSel", ak4Jets[1].pt,
        #                          hasTwoJets, EqBin(250, 0., 250.), title="subleading jet p_T",
        #                          xTitle="Leading Jet p_T (GeV/c^2)"))
        # plots.append(Plot.make1D("subleadingJetPt_hasElEl", ak4Jets[1].pt,
        #                          hasTwoJetsElEl, EqBin(250, 0., 250.), title="subleading jet p_T",
        #                          xTitle="Leading Jet p_T (GeV/c^2)"))
        # plots.append(Plot.make1D("subleadingJetPt_hasMuMu", ak4Jets[1].pt,
        #                          hasTwoJetsMuMu, EqBin(250, 0., 250.), title="subleading jet p_T",
        #                          xTitle="Leading Jet p_T (GeV/c^2)"))
        if not self.is_MC:
            yields.add(triggerSel, 'Trigger sel')
        else:
            yields.add(genWeightSel, 'gen weight sel')
        yields.add(hasElEl, 'two electrons')
        yields.add(hasTwoJetsElEl, 'two el. two jets')
        yields.add(hasMuMu, 'two muons')
        yields.add(hasTwoJetsMuMu, 'two muons two jets')

        return plots
