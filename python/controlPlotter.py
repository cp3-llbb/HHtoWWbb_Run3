
from bamboo.plots import Plot, CutFlowReport
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

import definitions as defs

from basePlotter import NanoBaseHHWWbb


class controlPlotter(NanoBaseHHWWbb):
    """"""

    def __init__(self, args):
        super(controlPlotter, self).__init__(args)
        self.channel = self.args.channel

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []
        yields = CutFlowReport("yields", printInLog=True, recursive=True)
        plots.append(yields)
        yields.add(noSel, 'No Selection')

        # Muons
        muon_conept = defs.muonConePt(tree.Muon)

        muons = op.sort(
            op.select(tree.Muon, lambda mu: defs.muonDef(mu)),
            lambda mu: -muon_conept[mu.idx]
        )

        fakeMuons = defs.muonFakeSel(muons)

        tightMuons = op.select(fakeMuons, lambda mu: defs.muonTightSel(mu))

        # Electrons
        electron_conept = defs.elConePt(tree.Electron)

        electrons = op.sort(
            op.select(tree.Electron, lambda el: op.AND(
                defs.elDef(el), electron_conept[el.idx] > 7)),
            lambda el: -electron_conept[el.idx]
        )
        # Cleaned Electrons
        clElectrons = defs.cleanElectrons(electrons, muons)

        # Fake Electrons
        fakeElectrons = defs.elFakeSel(clElectrons)

        tightElectrons = op.select(
            fakeElectrons, lambda el: defs.elTightSel(el))

        # AK8 Jets
        ak8Jets = op.sort(
            op.select(tree.FatJet, lambda jet: defs.ak8jetDef(jet)), lambda jet: -jet.pt)

        cleanedAK8Jets = op.select(ak8Jets, lambda jet: op.AND(op.rng_any(fakeElectrons, lambda el: op.deltaR(
            el.p4, jet.p4) > 0.8), op.rng_any(fakeMuons, lambda mu: op.deltaR(mu.p4, jet.p4) > 0.8)))

        ak8bJets = op.select(
            ak8Jets, lambda fatjet: fatjet.btagDeepB > 0.4184)  # 2018 WP

        # AK4 Jets
        ak4Jets = op.sort(
            op.select(tree.Jet, lambda jet: defs.ak4jetDef(jet)), lambda jet: -jet.pt)

        # remove jets within cone of DR<0.4 of leading leptons at each channel

        if self.channel == 'SL':
            def cleaningWithRespectToLeadingLepton(DR):
                return lambda jet: op.multiSwitch(
                    (op.AND(op.rng_len(fakeElectrons) >= 1, op.rng_len(
                        fakeMuons) == 0), op.deltaR(jet.p4, fakeElectrons[0].p4) >= DR),
                    (op.AND(op.rng_len(fakeElectrons) == 0, op.rng_len(
                        fakeMuons) >= 1), op.deltaR(jet.p4, fakeMuons[0].p4) >= DR),
                    (op.AND(op.rng_len(fakeMuons) >= 1, op.rng_len(fakeElectrons) >= 1), op.switch(
                        electron_conept[0] >= muon_conept[0],
                        op.deltaR(jet.p4, fakeElectrons[0].p4) >= DR,
                        op.deltaR(jet.p4, fakeMuons[0].p4) >= DR)),
                    op.c_bool(True)
                )
            cleanAk4Jets = cleaningWithRespectToLeadingLepton(0.4)

        if self.channel == 'DL':
            def cleaningWithRespectToLeadingLeptons(DR):
                return lambda j: op.multiSwitch(
                    # Only electrons
                    (op.AND(op.rng_len(fakeElectrons) >= 2, op.rng_len(fakeMuons) == 0),
                     op.AND(op.deltaR(j.p4, fakeElectrons[0].p4) >= DR, op.deltaR(j.p4, fakeElectrons[1].p4) >= DR)),
                    # Only muons
                    (op.AND(op.rng_len(fakeElectrons) == 0, op.rng_len(fakeMuons) >= 2),
                     op.AND(op.deltaR(j.p4, fakeMuons[0].p4) >= DR, op.deltaR(j.p4, fakeMuons[1].p4) >= DR)),
                    # One electron + one muon
                    (op.AND(op.rng_len(fakeElectrons) == 1, op.rng_len(fakeMuons) == 1),
                     op.AND(op.deltaR(j.p4, fakeElectrons[0].p4) >= DR, op.deltaR(j.p4, fakeMuons[0].p4) >= DR)),
                    # At least one electron + at least one muon
                    (op.AND(op.rng_len(fakeElectrons) >= 1, op.rng_len(fakeMuons) >= 1),
                     op.switch(
                        # Electron is the leading lepton
                        electron_conept[0] > muon_conept[0],
                        op.switch(op.rng_len(fakeElectrons) == 1,
                                  op.AND(op.deltaR(j.p4, fakeElectrons[0].p4) >= DR, op.deltaR(
                                      j.p4, fakeMuons[0].p4) >= DR),
                                  op.switch(electron_conept[1] > muon_conept[0],
                                            op.AND(op.deltaR(j.p4, fakeElectrons[0].p4) >= DR, op.deltaR(
                                                j.p4, fakeElectrons[1].p4) >= DR),
                                            op.AND(op.deltaR(j.p4, fakeElectrons[0].p4) >= DR, op.deltaR(j.p4, fakeMuons[0].p4) >= DR))),
                        # Muon is the leading lepton
                        op.switch(op.rng_len(fakeMuons) == 1,
                                  op.AND(op.deltaR(j.p4, fakeMuons[0].p4) >= DR, op.deltaR(
                                      j.p4, fakeElectrons[0].p4) >= DR),
                                  op.switch(muon_conept[1] > electron_conept[0],
                                            op.AND(op.deltaR(j.p4, fakeMuons[0].p4) >= DR, op.deltaR(
                                                j.p4, fakeMuons[1].p4) >= DR),
                                            op.AND(op.deltaR(j.p4, fakeMuons[0].p4) >= DR, op.deltaR(j.p4, fakeElectrons[0].p4) >= DR))))),
                    op.c_bool(True)
                )
            cleanAk4Jets = cleaningWithRespectToLeadingLeptons(0.4)
        
        clAK4Jets = op.select(ak4Jets, cleanAk4Jets)

        # Taus

        taus = defs.tauDef(tree.Tau)

        cleanedTaus = op.select(taus, lambda tau: op.AND(op.rng_any(fakeElectrons, lambda el: op.deltaR(
            el.p4, tau.p4) > 0.3), op.rng_any(fakeMuons, lambda mu: op.deltaR(mu.p4, tau.p4) > 0.3)))

        ### Di-leptonic channel ###

        # has exactly two leptons
        hasTwoL = noSel.refine('hasTwoL', cut=(
            op.OR(
                op.AND(op.rng_len(clElectrons) == 2, op.rng_len(muons) == 0,
                       clElectrons[0].charge != clElectrons[1].charge, clElectrons[0].pt > 25., clElectrons[1].pt > 15.),
                op.AND(op.rng_len(muons) == 2, op.rng_len(clElectrons) == 0,
                       muons[0].charge != muons[1].charge, muons[0].pt > 25., muons[1].pt > 15.),
                op.AND(op.rng_len(clElectrons) == 1, op.rng_len(muons) == 1,
                       clElectrons[0].charge != muons[0].charge,
                       op.switch(
                    clElectrons[0].pt >= muons[0].pt,
                    op.AND(clElectrons[0].pt > 25., muons[0].pt > 15.),
                    op.AND(clElectrons[0].pt > 15., muons[0].pt > 25.))
                ))
        ))

        # lepton channels
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
        # DL_resolved = hasTwoL.refine('DL_resolved', cut=(op.AND(op.rng_len(
        #     ak4Jets) >= 2, op.rng_len(ak4bJets) >= 1, op.rng_len(ak8Jets) == 0)))

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

        # ak4ak4bJetPair = op.combine((ak4Jets, ak4bJets), N=2, pred=lambda j1, j2:
        #                             op.deltaR(j1.p4, j2.p4) > 0.8)
        # firstJetPair = ak4ak4bJetPair[0]

        # boosted -> and at least one b-tagged ak8 jet and at least one ak4 jet outside the b-tagged ak8 jet
        SL_boosted = hasOneL.refine('SL_boosted', cut=(op.AND(
            op.rng_len(ak8bJets) >= 1,
            op.rng_len(ak4Jets) >= 1,
            op.deltaR(ak4Jets[0].p4, ak8bJets[0].p4) >= 1.2)
        ))
        # resolved -> and at least three ak4 jets with at least one b-tagged and no ak8 jets
        # SL_resolved = hasOneL.refine('SL_resolved', cut=(op.AND(op.rng_len(
        #     ak4Jets) >= 3, op.rng_len(ak4bJets) >= 1, op.rng_len(ak8Jets) == 0)
        # ))

        #############################################################################
        #                                 Plots                                     #
        #############################################################################
        plots.extend([
            Plot.make1D("nFakeElectrons", op.rng_len(fakeElectrons), noSel, EqBin(
                15, 0., 15.), xTitle="Number of fake electrons"),
            Plot.make1D("nFakeMuons", op.rng_len(fakeMuons), noSel, EqBin(
                15, 0., 15.), xTitle="Number of fake muons"),
            Plot.make1D("ncleanedAK8Jets", op.rng_len(cleanedAK8Jets), noSel, EqBin(
                15, 0., 15.), xTitle="Number of cleaned AK8 jets"),
            Plot.make1D("nTaus", op.rng_len(taus), noSel, EqBin(
                15, 0., 15.), xTitle="Number of taus"),
            Plot.make1D("nCleanedTaus", op.rng_len(cleanedTaus), noSel, EqBin(
                15, 0., 15.), xTitle="Number of cleaned taus"),
            Plot.make1D("nCleanedAK4Jets", op.rng_len(clAK4Jets), noSel, EqBin(
                15, 0., 15.), xTitle="Number of cleaned AK4 jets"),
        ])

        # Cutflow report
        yields.add(hasOneL, 'one lepton')
        yields.add(hasTwoL, 'two leptons')
        yields.add(DL_boosted, 'DL boosted')
        # yields.add(DL_resolved, 'DL resolved')
        yields.add(SL_boosted, 'SL boosted')
        # yields.add(SL_resolved, 'SL resolved')

        return plots
