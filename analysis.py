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
        # Muon
        addHLTPath('Muon', 'IsoMu24')
        addHLTPath('Muon', 'IsoMu27')
        addHLTPath('Muon', 'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8')
        # EGamma
        addHLTPath('EGamma', 'Ele32_WPTight_Gsf')
        addHLTPath('EGamma', 'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL')
        # MuonEG
        addHLTPath('MuonEG', 'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ')
        # SingleMuon
        addHLTPath('SingleMuon', 'IsoMu24')
        addHLTPath('SingleMuon', 'IsoMu27')
        # DoubleMuon
        addHLTPath('DoubleMuon', 'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8')

        return tree, noSel, backend, lumiArgs

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []
        yields = CutFlowReport("yields", printInLog=True, recursive=True)
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
        # and at least two ak4 jets
        hasTwoJetsElEl = hasElEl.refine(
            "hasTwoJetsElEl", cut=[op.rng_len(ak4Jets) >= 2])
        # and two b jets
        hasTwoBJetsElEl = hasTwoJetsElEl.refine(
            "hasTwoBJetsElEl", cut=[op.rng_len(ak4BJets) >= 2])
        # has at least one muon pair
        hasMuMu = noSel.refine("hasOSMuMu", cut=[op.rng_len(muons) >= 2,
                                                 muons[0].charge != muons[1].charge, muons[0].pt > 20., muons[1].pt > 10.])
        # and at least two ak4 jets
        hasTwoJetsMuMu = hasMuMu.refine(
            'hasTwoJetsMuMu', cut=[op.rng_len(ak4Jets) >= 2])
        # and two b jets
        hasTwoBJetsMuMu = hasTwoJetsMuMu.refine(
            'hasTwoBJetsMuMu', cut=[op.rng_len(ak4BJets) >= 2])
        # has at least one ak4 jet
        hasOneJet = noSel.refine('hasOneJet', cut=[op.rng_len(ak4Jets) >= 1])
        # has at least two ak4 jets
        hasTwoJets = noSel.refine('hasTwoJets', cut=[op.rng_len(ak4Jets) >= 2])

        ### Di-leptonic channel ###

        # has exactly two leptons
        hasTwoL = noSel.refine('hasTwoL', cut=(
            op.OR(
                op.AND(op.rng_len(electrons) == 2,
                       electrons[0].charge != electrons[1].charge, electrons[0].pt > 25., electrons[1].pt > 15.),
                op.AND(op.rng_len(muons) == 2,
                       muons[0].charge != muons[1].charge, muons[0].pt > 25., muons[1].pt > 15.),
                op.AND(op.rng_len(electrons) == 1, op.rng_len(muons)
                       == 1, electrons[0].charge != muons[0].charge, op.OR(op.AND(electrons[0].pt > 25., muons[0].pt > 15.), op.AND(electrons[0].pt > 15., muons[0].pt > 25.)))
            )
        ))

        emuPair = op.combine((electrons, muons), N=2,
                             pred=lambda l1, l2: l1.charge != l2.charge)
        eePair = op.combine(electrons, N=2, pred=lambda l1,
                            l2: l1.charge != l2.charge)
        mumuPair = op.combine(muons, N=2, pred=lambda l1,
                              l2: l1.charge != l2.charge)

        firstEMUpair = emuPair[0]
        firstEEpair = eePair[0]
        firstMUMUpair = mumuPair[0]
        # boosted -> and at least one b-tagged ak8 jet
        DL_boosted = hasTwoL.refine(
            'DL_boosted', cut=(op.rng_len(ak8BJets) >= 1))

        # resolved -> and at least two ak4 jets with at least one b-tagged and no ak8 jets
        DL_resolved = hasTwoL.refine('DL_resolved', cut=(op.AND(op.rng_len(
            ak4Jets) >= 2, op.rng_len(ak4BJets) >= 1, op.rng_len(ak8Jets) == 0)))

        ### Single-leptonic channel ###
        # has exactly one lepton
        hasOneL = noSel.refine('hasOneL', cut=(op.OR(op.AND(op.rng_len(
            electrons) == 1, electrons[0].pt > 32.), op.AND(op.rng_len(muons) == 1, muons[0].pt > 25.))))

        ak4jetPair = op.combine((ak4Jets, ak4BJets), N=2, pred=lambda j1, j2: op.AND(
            op.deltaR(j1.p4, j2.p4) > 0.8, op.AND(j1.pt > 25., j2.pt > 25.)))
        firstJetPair = ak4jetPair[0]

        ak4ak8Bpair = op.combine((ak4Jets, ak8BJets), N=2, pred=lambda j1, j2: op.AND(
            op.deltaR(j1.p4, j2.p4) > 1.2, op.AND(j1.pt > 25., j2.pt > 200.)))
        firstAK4AK8Bpair = ak4ak8Bpair[0]

        # boosted -> and at least one b-tagged ak8 jet and at least one ak4 jet outside the b-tagged ak8 jet
        SL_boosted = hasOneL.refine('SL_boosted', cut=(op.AND(op.rng_len(ak8BJets) >= 1, op.AND(
            op.rng_len(ak4Jets) >= 1, op.deltaR(firstAK4AK8Bpair[0].p4, firstAK4AK8Bpair[1].p4) > 1.2))))
        # resolved -> and at least three ak4 jets with at least one b-tagged and no ak8 jets
        SL_resolved = hasOneL.refine('SL_resolved', cut=(op.AND(op.rng_len(
            ak4Jets) >= 3, op.rng_len(ak4BJets) >= 1, op.rng_len(ak8Jets) == 0)))

        #############################################################################
        #                                 Plots                                     #
        #############################################################################
        plots.extend([
            # Plot.make1D("nEl_NoSel", op.rng_len(electrons), noSel, EqBin(
            #     10, 0., 10.), xTitle="Number of electrons"),
            # Plot.make1D("nEl_HasElEl", op.rng_len(electrons), hasElEl, EqBin(
            #     10, 0., 10.), xTitle="Number of electrons"),
            # Plot.make1D("nEl_hasTwoJetsElEl", op.rng_len(electrons), hasTwoJetsElEl, EqBin(
            #     10, 0., 10.), xTitle="Number of electrons"),
            # Plot.make1D("nMu_NoSel", op.rng_len(muons), noSel, EqBin(
            #     10, 0., 10.), xTitle="Number of muons"),
            # Plot.make1D("nMu_HasMuMu", op.rng_len(muons), hasMuMu, EqBin(
            #     10, 0., 10.), xTitle="Number of muons"),
            # Plot.make1D("nMu_hasTwoJetsMuMu", op.rng_len(muons), hasTwoJetsMuMu, EqBin(
            #     10, 0., 10.), xTitle="Number of muons"),
            # Plot.make1D("nJet_NoSel", op.rng_len(ak4Jets), noSel, EqBin(
            #     10, 0., 10.), xTitle="Number of jets"),
            # Plot.make1D("nJet_HasElEl", op.rng_len(ak4Jets), hasElEl, EqBin(
            #     10, 0., 10.), xTitle="Number of jets"),
            # Plot.make1D("nJet_HasMuMu", op.rng_len(ak4Jets), hasMuMu, EqBin(
            #     10, 0., 10.), xTitle="Number of jets"),
            # Plot.make1D("nJet_hasTwoJetsElEl", op.rng_len(ak4Jets), hasTwoJetsElEl, EqBin(
            #     10, 0., 10.), xTitle="Number of jets"),
            # Plot.make1D("nJet_hasTwoJetsMuMu", op.rng_len(ak4Jets), hasTwoJetsMuMu, EqBin(
            #     10, 0., 10.), xTitle="Number of jets"),
            # Plot.make1D("massZto2e", op.invariant_mass(electrons[0].p4, electrons[1].p4),
            #             hasElEl, EqBin(120, 40., 120.), title="mass of Z to 2e",
            #             xTitle="Invariant Mass of electrons (GeV/c^2)"),
            # Plot.make1D("massZto2e_hasTwoJets", op.invariant_mass(electrons[0].p4, electrons[1].p4),
            #             hasTwoJetsElEl, EqBin(120, 40., 120.), title="mass of Z to 2e",
            #             xTitle="Invariant Mass of electrons (GeV/c^2)"),
            # Plot.make1D("massZto2e_hasTwoBJets", op.invariant_mass(electrons[0].p4, electrons[1].p4),
            #             hasTwoBJetsElEl, EqBin(120, 40., 120.), title="mass of Z to 2e",
            #             xTitle="Invariant Mass of electrons (GeV/c^2)"),
            # Plot.make1D("massZto2mu", op.invariant_mass(muons[0].p4, muons[1].p4),
            #             hasMuMu, EqBin(120, 40., 120.), title="mass of Z to 2mu",
            #             xTitle="Invariant Mass of muons (GeV/c^2)"),
            # Plot.make1D("massZto2mu_hasTwoJets", op.invariant_mass(muons[0].p4, muons[1].p4),
            #             hasTwoJetsMuMu, EqBin(120, 40., 120.), title="mass of Z to 2mu",
            #             xTitle="Invariant Mass of muons (GeV/c^2)"),
            # Plot.make1D("massZto2mu_hasTwoBJets", op.invariant_mass(muons[0].p4, muons[1].p4),
            #             hasTwoBJetsMuMu, EqBin(120, 40., 120.), title="mass of Z to 2mu",
            #             xTitle="Invariant Mass of muons (GeV/c^2)"),
            # Plot.make1D("leadingJetPt_hasOneJet", ak4Jets[0].pt,
            #             hasOneJet, EqBin(250, 0., 250.), title="leading jet p_T",
            #             xTitle="Leading Jet p_T (GeV/c^2)"),
            # Plot.make1D("leadingJetPt_hasTwoJets", ak4Jets[0].pt,
            #             hasTwoJets, EqBin(250, 0., 250.), title="leading jet p_T",
            #             xTitle="Leading Jet p_T (GeV/c^2)"),
            # Plot.make1D("leadingJetPt_hasTwoJetsElEl", ak4Jets[0].pt,
            #             hasTwoJetsElEl, EqBin(250, 0., 250.), title="leading jet p_T",
            #             xTitle="Leading Jet p_T (GeV/c^2)"),
            # Plot.make1D("leadingJetPt_hasTwoJetsMuMu", ak4Jets[0].pt,
            #             hasTwoJetsMuMu, EqBin(250, 0., 250.), title="leading jet p_T",
            #             xTitle="Leading Jet p_T (GeV/c^2)"),
            # Plot.make1D("subleadingJetPt_NoSel", ak4Jets[1].pt,
            #             hasTwoJets, EqBin(250, 0., 250.), title="subleading jet p_T",
            #             xTitle="Leading Jet p_T (GeV/c^2)"),
            # Plot.make1D("subleadingJetPt_hasElEl", ak4Jets[1].pt,
            #             hasTwoJetsElEl, EqBin(250, 0., 250.), title="subleading jet p_T",
            #             xTitle="Leading Jet p_T (GeV/c^2)"),
            # Plot.make1D("subleadingJetPt_hasMuMu", ak4Jets[1].pt,
            #             hasTwoJetsMuMu, EqBin(250, 0., 250.), title="subleading jet p_T",
            #             xTitle="Leading Jet p_T (GeV/c^2)"),
            # Plot.make1D("nPU", tree.Pileup_nPU, hasElEl, EqBin(100, 0, 100), title="number of PU"),
            # Plot.make1D("nPU", tree.Pileup_nPU, hasMuMu, EqBin(100, 0, 100), title="number of PU")
            Plot.make1D("DL_InvM_emu_boosted", op.invariant_mass(firstEMUpair[0].p4, firstEMUpair[1].p4), DL_boosted, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (boosted) (GeV/c^2)"),
            Plot.make1D("DL_InvM_ee_boosted", op.invariant_mass(firstEEpair[0].p4, firstEEpair[1].p4), DL_boosted, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electrons (boosted) (GeV/c^2)"),
            Plot.make1D("DL_InvM_mumu_boosted", op.invariant_mass(firstMUMUpair[0].p4, firstMUMUpair[1].p4), DL_boosted, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muons (boosted) (GeV/c^2)"),

            Plot.make1D("DL_InvM_emu_resolved", op.invariant_mass(firstEMUpair[0].p4, firstEMUpair[1].p4), DL_resolved, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (boosted) (GeV/c^2)"),
            Plot.make1D("DL_InvM_ee_resolved", op.invariant_mass(firstEEpair[0].p4, firstEEpair[1].p4), DL_resolved, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electrons (boosted) (GeV/c^2)"),
            Plot.make1D("DL_InvM_mumu_resolved", op.invariant_mass(firstMUMUpair[0].p4, firstMUMUpair[1].p4), DL_resolved, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muons (boosted) (GeV/c^2)"),

            Plot.make1D("SL_InvM_jj_resolved", op.invariant_mass(firstJetPair[0].p4, firstJetPair[1].p4), SL_resolved, EqBin(
                160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),
            Plot.make1D("SL_InvM_jj_boosted", op.invariant_mass(firstAK4AK8Bpair[0].p4, firstAK4AK8Bpair[1].p4), SL_boosted, EqBin(
                160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)")
        ])

        yields.add(hasElEl, 'two electrons')
        yields.add(hasTwoJetsElEl, 'two el. two jets')
        yields.add(hasTwoBJetsElEl, 'two el. two Bjets')
        yields.add(hasMuMu, 'two muons')
        yields.add(hasTwoJetsMuMu, 'two muons two jets')
        yields.add(hasTwoBJetsMuMu, 'two muons two Bjets')
        yields.add(hasTwoL, 'two leptons')
        yields.add(DL_boosted, 'DL boosted')
        yields.add(DL_resolved, 'DL resolved')
        yields.add(SL_boosted, 'SL boosted')
        yields.add(SL_resolved, 'SL resolved')

        return plots
