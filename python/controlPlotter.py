
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.plots import Plot, CutFlowReport
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

import definitions as defs
from itertools import chain

from basePlotter import NanoBaseHHWWbb


class controlPlotter(NanoBaseHHWWbb):
    """"""

    def __init__(self, args):
        super(controlPlotter, self).__init__(args)

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []
        yields = CutFlowReport("yields", printInLog=True, recursive=True)
        plots.append(yields)
        yields.add(noSel, 'No Selection')

        # Muons
        muons = op.sort(
            op.select(tree.Muon, lambda mu: defs.muonDef(mu)),
            lambda mu: -defs.muonConePt(tree.Muon)[mu.idx]
        )
        # Electrons
        electrons = op.sort(
            op.select(tree.Electron, lambda el: defs.elDef(el)),
            lambda el: -defs.elConePt(tree.Electron)[el.idx]
        )

        # Cleaned Electrons
        clElectrons = defs.cleanElectrons(electrons, muons)
        # AK8 Jets
        ak8Jets = op.sort(
            op.select(tree.FatJet, lambda jet: defs.ak8jetDef(jet)), lambda jet: -jet.pt)

        ak8bJets = op.select(
            ak8Jets, lambda fatjet: fatjet.btagDeepB > 0.4184)  # 2018 WP

        # AK4 Jets(jet)
        ak4Jets = op.sort(
            op.select(tree.Jet, lambda jet: defs.ak4jetDef(jet)), lambda jet: -jet.pt)

        ak4bJets = op.select(
            ak4Jets, lambda jet: jet.btagDeepB > 0.2770)  # 2018 WP

        # Gen Weight and Triggers
        if self.is_MC:
            noSel = noSel.refine('genWeight', weight=tree.genWeight, cut=(
                op.OR(*chain.from_iterable(self.triggersPerPrimaryDataset.values()))))
        else:
            noSel = noSel.refine('trigger', cut=[makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggersPerPrimaryDataset)])

        # Selections

        # has at least one electron pair
        hasElEl = noSel.refine("hasOSElEl", cut=[op.rng_len(clElectrons) >= 2,
                                                 clElectrons[0].charge != clElectrons[1].charge, clElectrons[0].pt > 20., clElectrons[1].pt > 10.])
        # and at least two ak4 jets
        hasTwoJetsElEl = hasElEl.refine(
            "hasTwoJetsElEl", cut=[op.rng_len(ak4Jets) >= 2])
        # and two b jets
        hasTwoBJetsElEl = hasTwoJetsElEl.refine(
            "hasTwoBJetsElEl", cut=[op.rng_len(ak4bJets) >= 2])
        # has at least one muon pair
        hasMuMu = noSel.refine("hasOSMuMu", cut=[op.rng_len(muons) >= 2,
                                                 muons[0].charge != muons[1].charge, muons[0].pt > 20., muons[1].pt > 10.])
        # and at least two ak4 jets
        hasTwoJetsMuMu = hasMuMu.refine(
            'hasTwoJetsMuMu', cut=[op.rng_len(ak4Jets) >= 2])
        # and two b jets
        hasTwoBJetsMuMu = hasTwoJetsMuMu.refine(
            'hasTwoBJetsMuMu', cut=[op.rng_len(ak4bJets) >= 2])
        # has at least one ak4 jet
        hasOneJet = noSel.refine('hasOneJet', cut=[op.rng_len(ak4Jets) >= 1])
        # has at least two ak4 jets
        hasTwoJets = noSel.refine('hasTwoJets', cut=[op.rng_len(ak4Jets) >= 2])

        ### Di-leptonic channel ###

        # has exactly two leptons
        hasTwoL = noSel.refine('hasTwoL', cut=(
            op.OR(
                op.AND(op.rng_len(clElectrons) == 2, op.rng_len(muons) == 0,
                       clElectrons[0].charge != clElectrons[1].charge, clElectrons[0].pt > 25., clElectrons[1].pt > 15.),
                op.AND(op.rng_len(muons) == 2, op.rng_len(clElectrons) == 0,
                       muons[0].charge != muons[1].charge, muons[0].pt > 25., muons[1].pt > 15.),
                op.AND(op.rng_len(clElectrons) == 1, op.rng_len(muons) == 1,
                       clElectrons[0].charge != muons[0].charge, op.OR(op.AND(clElectrons[0].pt > 25., muons[0].pt > 15.), op.AND(clElectrons[0].pt > 15., muons[0].pt > 25.)))
            )
        ))

        emuPair = op.combine((clElectrons, muons), N=2,
                             pred=lambda el, mu: el.charge != mu.charge)
        eePair = op.combine(clElectrons, N=2, pred=lambda el1,
                            el2: el1.charge != el2.charge)
        mumuPair = op.combine(muons, N=2, pred=lambda mu1,
                              mu2: mu1.charge != mu2.charge)

        firstEMUpair = emuPair[0]
        firstEEpair = eePair[0]
        firstMUMUpair = mumuPair[0]
        # boosted -> and at least one b-tagged ak8 jet
        DL_boosted = hasTwoL.refine(
            'DL_boosted', cut=(op.rng_len(ak8bJets) >= 1))

        # resolved -> and at least two ak4 jets with at least one b-tagged and no ak8 jets
        DL_resolved = hasTwoL.refine('DL_resolved', cut=(op.AND(op.rng_len(
            ak4Jets) >= 2, op.rng_len(ak4bJets) >= 1, op.rng_len(ak8Jets) == 0)))

        ### Semi-leptonic channel ###
        # has exactly one lepton
        hasOneL = noSel.refine('hasOneL', cut=(op.OR(
            op.AND(
                op.rng_len(clElectrons) == 1,
                op.rng_len(muons) == 0,
                clElectrons[0].pt > 32.),
            op.AND(
                op.rng_len(muons) == 1,
                op.rng_len(clElectrons) == 0,
                muons[0].pt > 25.)
        )))

        ak4ak4bJetPair = op.combine((ak4Jets, ak4bJets), N=2, pred=lambda j1, j2:
                                    op.deltaR(j1.p4, j2.p4) > 0.8)
        firstJetPair = ak4ak4bJetPair[0]

        ak4ak8bPair = op.combine((ak4Jets, ak8bJets), N=2, pred=lambda ak4, ak8b: op.AND(
            op.deltaR(ak4.p4, ak8b.p4) >= 1.2))
        firstAK4AK8bPair = ak4ak8bPair[0]

        # boosted -> and at least one b-tagged ak8 jet and at least one ak4 jet outside the b-tagged ak8 jet
        SL_boosted = hasOneL.refine('SL_boosted', cut=(op.AND(
            op.rng_len(ak8bJets) >= 1,
            op.rng_len(ak4Jets) >= 1,
            op.deltaR(ak4Jets[0].p4, ak8bJets[0].p4) >= 1.2)
        ))
        # resolved -> and at least three ak4 jets with at least one b-tagged and no ak8 jets
        SL_resolved = hasOneL.refine('SL_resolved', cut=(op.AND(op.rng_len(
            ak4Jets) >= 3, op.rng_len(ak4bJets) >= 1, op.rng_len(ak8Jets) == 0)
        ))

        #############################################################################
        #                                 Plots                                     #
        #############################################################################
        plots.extend([
            # Plot.make1D("nEl_NoSel", op.rng_len(clElectrons), noSel, EqBin(
            #     10, 0., 10.), xTitle="Number of Electrons"),
            # Plot.make1D("nEl_HasElEl", op.rng_len(clElectrons), hasElEl, EqBin(
            #     10, 0., 10.), xTitle="Number of Electrons"),
            # Plot.make1D("nEl_hasTwoJetsElEl", op.rng_len(clElectrons), hasTwoJetsElEl, EqBin(
            #     10, 0., 10.), xTitle="Number of Electrons"),
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
            # Plot.make1D("massZto2e", op.invariant_mass(clElectrons[0].p4, clElectrons[1].p4),
            #             hasElEl, EqBin(120, 40., 120.), title="mass of Z to 2e",
            #             xTitle="Invariant Mass of electrons (GeV/c^2)"),
            # Plot.make1D("massZto2e_hasTwoJets", op.invariant_mass(clElectrons[0].p4, clElectrons[1].p4),
            #             hasTwoJetsElEl, EqBin(120, 40., 120.), title="mass of Z to 2e",
            #             xTitle="Invariant Mass of electrons (GeV/c^2)"),
            # Plot.make1D("massZto2e_hasTwoBJets", op.invariant_mass(clElectrons[0].p4, clElectrons[1].p4),
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

            Plot.make1D("DL_InvM_emu_boosted", op.invariant_mass(firstEMUpair[0].p4, firstEMUpair[1].p4), DL_boosted, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (boosted) (GeV/c^2)"),
            Plot.make1D("DL_InvM_ee_boosted", op.invariant_mass(firstEEpair[0].p4, firstEEpair[1].p4), DL_boosted, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electrons (boosted) (GeV/c^2)"),
            Plot.make1D("DL_InvM_mumu_boosted", op.invariant_mass(firstMUMUpair[0].p4, firstMUMUpair[1].p4), DL_boosted, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muons (boosted) (GeV/c^2)"),
            Plot.make1D("DL_InvM_jj_boosted", op.invariant_mass(ak8Jets[0].subJet1.p4, ak8Jets[0].subJet2.p4), DL_boosted, EqBin(
                160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),

            Plot.make1D("DL_InvM_emu_resolved", op.invariant_mass(firstEMUpair[0].p4, firstEMUpair[1].p4), DL_resolved, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (resolved) (GeV/c^2)"),
            Plot.make1D("DL_InvM_ee_resolved", op.invariant_mass(firstEEpair[0].p4, firstEEpair[1].p4), DL_resolved, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electrons (resolved) (GeV/c^2)"),
            Plot.make1D("DL_InvM_mumu_resolved", op.invariant_mass(firstMUMUpair[0].p4, firstMUMUpair[1].p4), DL_resolved, EqBin(
                160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muons (resolved) (GeV/c^2)"),
            Plot.make1D("DL_InvM_jj_resolved", op.invariant_mass(firstJetPair[0].p4, firstJetPair[1].p4), SL_resolved, EqBin(
                160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),

            Plot.make1D("SL_InvM_jj_resolved", op.invariant_mass(firstJetPair[0].p4, firstJetPair[1].p4), SL_resolved, EqBin(
                160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),
            Plot.make1D("SL_InvM_jj_boosted", op.invariant_mass(ak8bJets[0].subJet1.p4, ak8bJets[0].subJet2.p4), SL_boosted, EqBin(
                160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)")
        ])

        # Cutflow report
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
