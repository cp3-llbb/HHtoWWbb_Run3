
from bamboo.plots import Plot, CutFlowReport
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

import definitions as defs

from basePlotter import NanoBaseHHWWbb
from selections import makeDLSelection, makeSLSelection


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

        # lepton cone-pt definitions
        self.muon_conept = defs.muonConePt(tree.Muon)
        self.electron_conept = defs.elConePt(tree.Electron)

        # lepton definitions sorted by their cone-pt
        self.muons = op.sort(defs.muonDef(tree.Muon), lambda mu: -self.muon_conept[mu.idx])
        self.electrons = op.sort(defs.elDef(tree.Electron), lambda el: -self.electron_conept[el.idx])

        # cleaning electrons wrt muons
        self.clElectrons = defs.cleanElectrons(self.electrons, self.muons)

        # Fakeable leptons
        self.fakeMuons = defs.muonFakeSel(self.muons)
        self.fakeElectrons = defs.elFakeSel(self.clElectrons)

        # tight leptons
        self.tightMuons = defs.muonTightSel(self.fakeMuons)
        self.tightElectrons = defs.elTightSel(self.fakeElectrons)

        # Taus
        taus = defs.tauDef(tree.Tau)
        self.cleanedTaus = defs.cleanTaus(taus, self.fakeElectrons, self.fakeMuons)

        # AK4 Jets sorted by their pt
        ak4JetsPreSel = op.sort(defs.ak4jetDef(tree.Jet), lambda jet: -jet.pt)

        # remove jets within cone of DR<0.4 of leading leptons at each channel
        if self.channel == 'SL':
            def cleaningWithRespectToLeadingLepton(DR):
                return lambda jet: op.multiSwitch(
                    (op.AND(op.rng_len(self.fakeElectrons) >= 1, op.rng_len(
                        self.fakeMuons) == 0), op.deltaR(jet.p4, self.fakeElectrons[0].p4) >= DR),
                    (op.AND(op.rng_len(self.fakeElectrons) == 0, op.rng_len(
                        self.fakeMuons) >= 1), op.deltaR(jet.p4, self.fakeMuons[0].p4) >= DR),
                    (op.AND(op.rng_len(self.fakeMuons) >= 1, op.rng_len(self.fakeElectrons) >= 1), op.switch(
                        self.electron_conept[0] >= self.muon_conept[0],
                        op.deltaR(jet.p4, self.fakeElectrons[0].p4) >= DR,
                        op.deltaR(jet.p4, self.fakeMuons[0].p4) >= DR)),
                    op.c_bool(True)
                )
            cleanAk4Jets = cleaningWithRespectToLeadingLepton(0.4)

        if self.channel == 'DL':
            def cleaningWithRespectToLeadingLeptons(DR):
                return lambda j: op.multiSwitch(
                    # Only electrons
                    (op.AND(op.rng_len(self.fakeElectrons) >= 2, op.rng_len(self.fakeMuons) == 0),
                     op.AND(op.deltaR(j.p4, self.fakeElectrons[0].p4) >= DR, op.deltaR(j.p4, self.fakeElectrons[1].p4) >= DR)),
                    # Only muons
                    (op.AND(op.rng_len(self.fakeElectrons) == 0, op.rng_len(self.fakeMuons) >= 2),
                     op.AND(op.deltaR(j.p4, self.fakeMuons[0].p4) >= DR, op.deltaR(j.p4, self.fakeMuons[1].p4) >= DR)),
                    # One electron + one muon
                    (op.AND(op.rng_len(self.fakeElectrons) == 1, op.rng_len(self.fakeMuons) == 1),
                     op.AND(op.deltaR(j.p4, self.fakeElectrons[0].p4) >= DR, op.deltaR(j.p4, self.fakeMuons[0].p4) >= DR)),
                    # At least one electron + at least one muon
                    (op.AND(op.rng_len(self.fakeElectrons) >= 1, op.rng_len(self.fakeMuons) >= 1),
                     op.switch(
                        # Electron is the leading lepton
                        self.electron_conept[0] > self.muon_conept[0],
                        op.switch(op.rng_len(self.fakeElectrons) == 1,
                                  op.AND(op.deltaR(j.p4, self.fakeElectrons[0].p4) >= DR, op.deltaR(
                                      j.p4, self.fakeMuons[0].p4) >= DR),
                                  op.switch(self.electron_conept[1] > self.muon_conept[0],
                                            op.AND(op.deltaR(j.p4, self.fakeElectrons[0].p4) >= DR, op.deltaR(
                                                j.p4, self.fakeElectrons[1].p4) >= DR),
                                            op.AND(op.deltaR(j.p4, self.fakeElectrons[0].p4) >= DR, op.deltaR(j.p4, self.fakeMuons[0].p4) >= DR))),
                        # Muon is the leading lepton
                        op.switch(op.rng_len(self.fakeMuons) == 1,
                                  op.AND(op.deltaR(j.p4, self.fakeMuons[0].p4) >= DR, op.deltaR(
                                      j.p4, self.fakeElectrons[0].p4) >= DR),
                                  op.switch(self.muon_conept[1] > self.electron_conept[0],
                                            op.AND(op.deltaR(j.p4, self.fakeMuons[0].p4) >= DR, op.deltaR(
                                                j.p4, self.fakeMuons[1].p4) >= DR),
                                            op.AND(op.deltaR(j.p4, self.fakeMuons[0].p4) >= DR, op.deltaR(j.p4, self.fakeElectrons[0].p4) >= DR))))),
                    op.c_bool(True)
                )
            cleanAk4Jets = cleaningWithRespectToLeadingLeptons(0.4)

        self.ak4Jets = op.select(ak4JetsPreSel, cleanAk4Jets)
        self.ak4JetsByBtagScore = op.sort(self.ak4Jets, lambda j: -j.btagDeepFlavB)

        # bTagging for ak4 jets
        def ak4BtagLooseSel(jet): return jet.btagDeepFlavB > 0.0494
        def ak4BtagSel(jet): return jet.btagDeepFlavB > 0.2770
        def ak4NoBtagSel(jet): return jet.btagDeepFlavB <= 0.2770

        self.ak4BJets = op.select(self.ak4Jets, ak4BtagSel)
        self.ak4BJetsLoose = op.select(self.ak4Jets, ak4BtagLooseSel)
        self.ak4LightJetsByPt = op.select(self.ak4Jets, ak4NoBtagSel)
        self.ak4LightJetsByBtagScore = op.sort(
            self.ak4LightJetsByPt, lambda jet: -jet.btagDeepFlavB)
        self.remainingJets = op.select(
            self.ak4LightJetsByPt, lambda jet: jet.idx != self.ak4LightJetsByBtagScore[0].idx)

        def makeJetPairs(jets): return op.combine(
            jets, N=2, pred=lambda j1, j2: j1.pt > j2.pt, samePred=lambda j1, j2: j1.idx != j2.idx)
        # --------------------------------------------- #
        bJetsByScore = self.ak4JetsByBtagScore[:op.min(op.rng_len(
            self.ak4JetsByBtagScore), op.static_cast("std::size_t", op.c_int(2)))]
        probableWJets = op.select(self.ak4Jets, lambda jet: op.NOT(
            op.rng_any(bJetsByScore, lambda bjet: jet.idx == bjet.idx)))
        wJetsByPt = probableWJets[:op.min(op.rng_len(
            probableWJets), op.static_cast("std::size_t", op.c_int(2)))]

        def passWMassCutSel(wjets): return op.switch(op.rng_len(wjets) == 2, op.abs(
            op.invariant_mass(wjets[0].p4, wjets[1].p4)-80.4) < op.c_float(15.0), op.c_bool(False))

        # AK8 Jets
        self.ak8JetsDef = defs.ak8jetDef(tree.FatJet)

        if self.channel == 'SL': # sorted by btag score
            ak8JetsPreSel = op.sort(self.ak8JetsDef, lambda j: -j.btagDeepB)
        if self.channel == 'DL': # sorted by pt
            ak8JetsPreSel = op.sort(self.ak8JetsDef, lambda j: -j.pt)

        # cleaning ak8 jets wrt to leptons
        if self.channel == 'SL':
            cleanAk8Jets = cleaningWithRespectToLeadingLepton(0.8)
        if self.channel == 'DL':
            cleanAk8Jets = cleaningWithRespectToLeadingLeptons(0.8)

        self.ak8Jets = op.select(ak8JetsPreSel, cleanAk8Jets)

        # 2018 DeepJet WP
        def subjetBtag(subjet): return subjet.btagDeepB > 0.4184

        def ak8Btag(fatjet): return op.OR(op.AND(fatjet.subJet1.pt >= 30, subjetBtag(fatjet.subJet1)),
                                          op.AND(fatjet.subJet2.pt >= 30, subjetBtag(fatjet.subJet2)))

        def ak8noBtag(fatjet): return op.NOT(op.OR(op.AND(fatjet.subJet1.pt >= 30, subjetBtag(fatjet.subJet1)),
                                                   op.AND(fatjet.subJet2.pt >= 30, subjetBtag(fatjet.subJet2))))

        def ak8Btag_bothSubJets(fatjet): return op.AND(op.AND(fatjet.subJet1.pt >= 30, subjetBtag(fatjet.subJet1)),
                                                       op.AND(fatjet.subJet2.pt >= 30, subjetBtag(fatjet.subJet2)))

        self.ak8BJets = op.select(self.ak8Jets, ak8Btag)
        self.ak8nonBJets = op.select(self.ak8Jets, ak8noBtag)
        # Ak4 Jet Collection cleaned from Ak8b #

        def cleanAk4FromAk8b(ak4j): return op.AND(op.rng_len(
            self.ak8BJets) > 0, op.deltaR(ak4j.p4, self.ak8BJets[0].p4) > 1.2)
        self.ak4JetsCleanedFromAk8b = op.select(self.ak4Jets, cleanAk4FromAk8b)

        # used as a BDT input for SemiBoosted category
        def btaggedSubJets(fjet): return op.switch(
            ak8Btag_bothSubJets(fjet), op.c_float(2.0), op.c_float(1.0))
        nMediumBTaggedSubJets = op.rng_sum(self.ak8BJets, btaggedSubJets)

        if self.channel == 'DL':

            DL_boosted_ee, DL_boosted_mumu,\
            DL_boosted_emu, DL_resolved_1b_ee,\
            DL_resolved_1b_mumu, DL_resolved_1b_emu,\
            DL_resolved_2b_ee, DL_resolved_2b_mumu,\
            DL_resolved_2b_emu = makeDLSelection(self, noSel)

            yields.add(DL_boosted_ee, 'DL boosted ee')
            yields.add(DL_boosted_mumu, 'DL boosted mumu')
            yields.add(DL_boosted_emu, 'DL boosted emu')
            yields.add(DL_resolved_1b_ee, 'DL resolved_1b_ee')
            yields.add(DL_resolved_1b_mumu, 'DL resolved_1b_mumu')
            yields.add(DL_resolved_1b_emu, 'DL resolved_1b_emu')
            yields.add(DL_resolved_2b_ee, 'DL resolved_2b_ee')
            yields.add(DL_resolved_2b_mumu, 'DL resolved_2b_mumu')
            yields.add(DL_resolved_2b_emu, 'DL resolved_2b_emu')

        if self.channel == 'SL':
            SL_resolved, SL_resolved_e,\
            SL_resolved_mu, SL_boosted,\
            SL_boosted_e, SL_boosted_mu = makeSLSelection(self, noSel)

            yields.add(SL_boosted, 'SL boosted')
            yields.add(SL_boosted_e, 'SL boosted_e')
            yields.add(SL_boosted_mu, 'SL boosted_mu')
            yields.add(SL_resolved, 'SL resolved')
            yields.add(SL_resolved_e, 'SL resolved_e')
            yields.add(SL_resolved_mu, 'SL resolved_mu')

        #############################################################################
        #                                 Plots                                     #
        #############################################################################
        if self.channel == 'DL':
            plots.extend([
                # DL boosted plots
                Plot.make1D("DL_boosted_fatJet_pt_ee", self.ak8Jets[0].pt, DL_boosted_ee, EqBin(
                    400, 200, 1000), title="pT(j1)", xTitle="pT(ak8jet) (GeV/c)"),
                Plot.make1D("DL_boosted_fatJet_pt_mumu", self.ak8Jets[0].pt, DL_boosted_mumu, EqBin(
                    400, 200, 1000), title="pT(j1)", xTitle="pT(ak8jet) (GeV/c)"),
                Plot.make1D("DL_boosted_fatJet_pt_emu", self.ak8Jets[0].pt, DL_boosted_emu, EqBin(
                    400, 200, 1000), title="pT(j1)", xTitle="pT(ak8jet) (GeV/c)"),

                Plot.make1D("DL_boosted_subjet1_pt_ee", self.ak8Jets[0].subJet1.pt, DL_boosted_ee, EqBin(
                    250, 0, 500), title=" pT(j1 subjet1)", xTitle="pT(j1 subjet1) (GeV/c)"),
                Plot.make1D("DL_boosted_subjet1_pt_mumu", self.ak8Jets[0].subJet1.pt, DL_boosted_mumu, EqBin(
                    250, 0, 500), title=" pT(j1 subjet1)", xTitle="pT(j1 subjet1) (GeV/c)"),
                Plot.make1D("DL_boosted_subjet1_pt_emu", self.ak8Jets[0].subJet1.pt, DL_boosted_emu, EqBin(
                    250, 0, 500), title=" pT(j1 subjet1)", xTitle="pT(j1 subjet1) (GeV/c)"),
                Plot.make1D("DL_boosted_subjet2_pt_ee", self.ak8Jets[0].subJet2.pt, DL_boosted_ee, EqBin(
                    250, 0, 500), title=" pT(j1 subjet2)", xTitle="pT(j1 subjet2) (GeV/c)"),
                Plot.make1D("DL_boosted_subjet2_pt_mumu", self.ak8Jets[0].subJet2.pt, DL_boosted_mumu, EqBin(
                    250, 0, 500), title=" pT(j1 subjet2)", xTitle="pT(j1 subjet2) (GeV/c)"),
                Plot.make1D("DL_boosted_subjet2_pt_emu", self.ak8Jets[0].subJet2.pt, DL_boosted_emu, EqBin(
                    250, 0, 500), title=" pT(j1 subjet2)", xTitle="pT(j1 subjet2) (GeV/c)"),
                Plot.make1D("DL_boosted_fatJet_eta_ee", self.ak8Jets[0].eta, DL_boosted_ee, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(ak8jet)"),
                Plot.make1D("DL_boosted_fatJet_eta_mumu", self.ak8Jets[0].eta, DL_boosted_mumu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(ak8jet)"),
                Plot.make1D("DL_boosted_fatJet_eta_emu", self.ak8Jets[0].eta, DL_boosted_emu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(ak8jet)"),
                Plot.make1D("DL_boosted_InvM_ee", op.invariant_mass(self.ElElDileptonPreSel[0][0].p4, self.ElElDileptonPreSel[0][1].p4), DL_boosted_ee, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electrons (boosted) (GeV/c^2)"),
                Plot.make1D("DL_boosted_InvM_mumu", op.invariant_mass(self.MuMuDileptonPreSel[0][0].p4, self.MuMuDileptonPreSel[0][1].p4), DL_boosted_mumu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muons (boosted) (GeV/c^2)"),
                Plot.make1D("DL_boosted_InvM_emu", op.invariant_mass(self.ElMuDileptonPreSel[0][0].p4, self.ElMuDileptonPreSel[0][1].p4), DL_boosted_emu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (boosted) (GeV/c^2)"),
                Plot.make1D("DL_boosted_InvM_jj_ee", op.invariant_mass(self.ak8Jets[0].subJet1.p4, self.ak8Jets[0].subJet2.p4), DL_boosted_ee, EqBin(
                    160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),
                Plot.make1D("DL_boosted_InvM_jj_mumu", op.invariant_mass(self.ak8Jets[0].subJet1.p4, self.ak8Jets[0].subJet2.p4), DL_boosted_mumu, EqBin(
                    160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),
                Plot.make1D("DL_boosted_InvM_jj_emu", op.invariant_mass(self.ak8Jets[0].subJet1.p4, self.ak8Jets[0].subJet2.p4), DL_boosted_emu, EqBin(
                    160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),

                # DL resolved 1b plots
                Plot.make1D("DL_resolved_1b_ee_nJets", op.rng_len(self.ak4Jets), DL_resolved_1b_ee, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_1b_ee_InvM_leadingJet_pt", self.ak4Jets[0].pt, DL_resolved_1b_ee, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_ee_InvM_subleadingJet_pt", self.ak4Jets[1].pt, DL_resolved_1b_ee, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_ee_InvM_leadingJet_eta", self.ak4Jets[0].eta, DL_resolved_1b_ee, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_1b_ee_InvM_subleadingJet_eta", self.ak4Jets[1].eta, DL_resolved_1b_ee, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_1b_ee_DR_jets", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_1b_ee, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_1b_ee_InvM_ee", op.invariant_mass(self.ElElDileptonPreSel[0][0].p4, self.ElElDileptonPreSel[0][1].p4), DL_resolved_1b_ee, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electrons (resolved) (GeV/c^2)"),

                Plot.make1D("DL_resolved_1b_mumu_nJets", op.rng_len(self.ak4Jets), DL_resolved_1b_mumu, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_1b_mumu_InvM_leadingJet_pt", self.ak4Jets[0].pt, DL_resolved_1b_mumu, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_mumu_InvM_subleadingJet_pt", self.ak4Jets[1].pt, DL_resolved_1b_mumu, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_mumu_InvM_leadingJet_eta", self.ak4Jets[0].eta, DL_resolved_1b_mumu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_1b_mumu_InvM_subleadingJet_eta", self.ak4Jets[1].eta, DL_resolved_1b_mumu, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_1b_mumu_DR_jets", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_1b_mumu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_1b_mumu_InvM_mumu", op.invariant_mass(self.MuMuDileptonPreSel[0][0].p4, self.MuMuDileptonPreSel[0][1].p4), DL_resolved_1b_mumu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muon (resolved) (GeV/c^2)"),

                Plot.make1D("DL_resolved_1b_emu_nJets", op.rng_len(self.ak4Jets), DL_resolved_1b_emu, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_1b_emu_InvM_leadingJet_pt", self.ak4Jets[0].pt, DL_resolved_1b_emu, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_emu_InvM_subleadingJet_pt", self.ak4Jets[1].pt, DL_resolved_1b_emu, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_emu_InvM_leadingJet_eta", self.ak4Jets[0].eta, DL_resolved_1b_emu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_1b_emu_InvM_subleadingJet_eta", self.ak4Jets[1].eta, DL_resolved_1b_emu, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_1b_emu_DR_jets", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_1b_emu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_1b_emu_InvM_emu", op.invariant_mass(self.ElMuDileptonPreSel[0][0].p4, self.ElMuDileptonPreSel[0][1].p4), DL_resolved_1b_emu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muon (resolved) (GeV/c^2)"),


                # DL resolved 2b plots
                Plot.make1D("DL_resolved_2b_ee_nJets", op.rng_len(self.ak4Jets), DL_resolved_2b_ee, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_2b_ee_InvM_leadingJet_pt", self.ak4Jets[0].pt, DL_resolved_2b_ee, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_ee_InvM_subleadingJet_pt", self.ak4Jets[1].pt, DL_resolved_2b_ee, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_ee_InvM_leadingJet_eta", self.ak4Jets[0].eta, DL_resolved_2b_ee, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_2b_ee_InvM_subleadingJet_eta", self.ak4Jets[1].eta, DL_resolved_2b_ee, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_2b_ee_DR_jets", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_2b_ee, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_2b_ee_InvM_ee", op.invariant_mass(self.ElElDileptonPreSel[0][0].p4, self.ElElDileptonPreSel[0][1].p4), DL_resolved_2b_ee, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electrons (resolved) (GeV/c^2)"),

                Plot.make1D("DL_resolved_2b_mumu_nJets", op.rng_len(self.ak4Jets), DL_resolved_2b_mumu, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_2b_mumu_InvM_leadingJet_pt", self.ak4Jets[0].pt, DL_resolved_2b_mumu, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_mumu_InvM_subleadingJet_pt", self.ak4Jets[1].pt, DL_resolved_2b_mumu, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_mumu_InvM_leadingJet_eta", self.ak4Jets[0].eta, DL_resolved_2b_mumu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_2b_mumu_InvM_subleadingJet_eta", self.ak4Jets[1].eta, DL_resolved_2b_mumu, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_2b_mumu_DR_jets", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_2b_mumu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_2b_mumu_InvM_mumu", op.invariant_mass(self.MuMuDileptonPreSel[0][0].p4, self.MuMuDileptonPreSel[0][1].p4), DL_resolved_2b_mumu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muons (resolved) (GeV/c^2)"),

                Plot.make1D("DL_resolved_2b_emu_nJets", op.rng_len(self.ak4Jets), DL_resolved_2b_emu, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_2b_emu_InvM_leadingJet_pt", self.ak4Jets[0].pt, DL_resolved_2b_emu, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_emu_InvM_subleadingJet_pt", self.ak4Jets[1].pt, DL_resolved_2b_emu, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_emu_InvM_leadingJet_eta", self.ak4Jets[0].eta, DL_resolved_2b_emu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_2b_emu_InvM_subleadingJet_eta", self.ak4Jets[1].eta, DL_resolved_2b_emu, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_2b_emu_DR_jets", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_2b_emu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_2b_emu_InvM_emu", op.invariant_mass(self.ElMuDileptonPreSel[0][0].p4, self.ElMuDileptonPreSel[0][1].p4), DL_resolved_2b_emu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pairs (resolved) (GeV/c^2)")

            ])
        if self.channel == "SL":
            plots.extend([
                # SL boosted plots
                Plot.make1D("SL_boosted_fatJet_pt", self.ak8BJets[0].pt, SL_boosted, EqBin(
                    400, 200, 1000), title="pT(j)", xTitle="pT(j) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet1_pt", self.ak8BJets[0].subJet1.pt, SL_boosted, EqBin(
                    250, 0, 500), title=" pT(subjet1)", xTitle="pT(subjet1) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet2_pt", self.ak8BJets[0].subJet2.pt, SL_boosted, EqBin(
                    250, 0, 500), title=" pT(subjet2)", xTitle="pT(subjet2) (GeV/c)"),
                Plot.make1D("SL_boosted_fatJet_eta", self.ak8BJets[0].eta, SL_boosted, EqBin(
                    80, -3, 3), title="eta(j)", xTitle="eta(j)"),
                Plot.make1D("SL_boosted_subjet1_eta", self.ak8BJets[0].subJet1.eta, SL_boosted, EqBin(
                    80, -3, 3), title="eta(subjet1)", xTitle="eta(subjet1)"),
                Plot.make1D("SL_boosted_subjet2_eta", self.ak8BJets[0].subJet2.eta, SL_boosted, EqBin(
                    80, -3, 3), title="eta(subjet2)", xTitle="eta(subjet2)"),
                Plot.make1D("SL_boosted_InvM_jj", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), SL_boosted, EqBin(
                    160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^2)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet1_eta", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].eta], SL_boosted, [
                    EqBin(160, 40., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet1 eta", xTitle="Invariant Mass of jets (GeV/c^2)", yTitle="eta(j1)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet2_eta", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].subJet2.eta], SL_boosted, [
                    EqBin(160, 40., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet2 eta", xTitle="Invariant Mass of jets (GeV/c^2)", yTitle="eta(j2)"),

                # SL boosted electron final state plots
                Plot.make1D("SL_boosted_fatJet_pt_e", self.ak8BJets[0].pt, SL_boosted_e, EqBin(
                    400, 200, 1000), title="pT(j)", xTitle="pT(j) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet1_pt_e", self.ak8BJets[0].subJet1.pt, SL_boosted_e, EqBin(
                    250, 0, 500), title=" pT(subjet1)", xTitle="pT(subjet1) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet2_pt_e", self.ak8BJets[0].subJet2.pt, SL_boosted_e, EqBin(
                    250, 0, 500), title=" pT(subjet2)", xTitle="pT(subjet2) (GeV/c)"),
                Plot.make1D("SL_boosted_fatJet_eta_e", self.ak8BJets[0].eta, SL_boosted_e, EqBin(
                    80, -3, 3), title="eta(j)", xTitle="eta(j)"),
                Plot.make1D("SL_boosted_subjet1_eta_e", self.ak8BJets[0].subJet1.eta, SL_boosted_e, EqBin(
                    80, -3, 3), title="eta(subjet1)", xTitle="eta(subjet1)"),
                Plot.make1D("SL_boosted_subjet2_eta_e", self.ak8BJets[0].subJet2.eta, SL_boosted_e, EqBin(
                    80, -3, 3), title="eta(subjet2)", xTitle="eta(subjet2)"),
                Plot.make1D("SL_boosted_InvM_jj_e", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), SL_boosted_e, EqBin(
                    160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^2)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet1_eta_e", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].eta], SL_boosted_e, [
                    EqBin(160, 40., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet1 eta", xTitle="Invariant Mass of jets (GeV/c^2)", yTitle="eta(j1)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet2_eta_e", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].subJet2.eta], SL_boosted_e, [
                    EqBin(160, 40., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet2 eta", xTitle="Invariant Mass of jets (GeV/c^2)", yTitle="eta(j2)"),

                # SL boosted muon final state plots
                Plot.make1D("SL_boosted_fatJet_pt_mu", self.ak8BJets[0].pt, SL_boosted_mu, EqBin(
                    400, 200, 1000), title="pT(j)", xTitle="pT(j) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet1_pt_mu", self.ak8BJets[0].subJet1.pt, SL_boosted_mu, EqBin(
                    250, 0, 500), title=" pT(subjet1)", xTitle="pT(subjet1) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet2_pt_mu", self.ak8BJets[0].subJet2.pt, SL_boosted_mu, EqBin(
                    250, 0, 500), title=" pT(subjet2)", xTitle="pT(subjet2) (GeV/c)"),
                Plot.make1D("SL_boosted_fatJet_eta_mu", self.ak8BJets[0].eta, SL_boosted_mu, EqBin(
                    80, -3, 3), title="eta(j)", xTitle="eta(j)"),
                Plot.make1D("SL_boosted_subjet1_eta_mu", self.ak8BJets[0].subJet1.eta, SL_boosted_mu, EqBin(
                    80, -3, 3), title="eta(subjet1)", xTitle="eta(subjet1)"),
                Plot.make1D("SL_boosted_subjet2_eta_mu", self.ak8BJets[0].subJet2.eta, SL_boosted_mu, EqBin(
                    80, -3, 3), title="eta(subjet2)", xTitle="eta(subjet2)"),
                Plot.make1D("SL_boosted_InvM_jj_mu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), SL_boosted_mu, EqBin(
                    160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^2)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet1_eta_mu", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].eta], SL_boosted_mu, [
                    EqBin(160, 40., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet1 eta", xTitle="Invariant Mass of jets (GeV/c^2)", yTitle="eta(j1)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet2_eta_mu", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].subJet2.eta], SL_boosted_mu, [
                    EqBin(160, 40., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet2 eta", xTitle="Invariant Mass of jets (GeV/c^2)", yTitle="eta(j2)"),

                # SL resolved plots
                Plot.make1D("SL_resolved_nJets", op.rng_len(self.ak4BJets), SL_resolved, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("SL_resolved_InvM_leadingJet_pt", self.ak4BJets[0].pt, SL_resolved, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_leadingJet_eta", self.ak4BJets[0].eta, SL_resolved, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_pt", self.ak4BJets[1].pt, SL_resolved, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_eta", self.ak4BJets[1].eta, SL_resolved, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("SL_resolved_DR_jets", op.deltaR(self.ak4BJets[0].p4, self.ak4BJets[1].p4), SL_resolved, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                # SL resolved electron final state plots
                Plot.make1D("SL_resolved_nJets_e", op.rng_len(self.ak4BJets), SL_boosted_e, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("SL_resolved_InvM_leadingJet_pt_e", self.ak4BJets[0].pt, SL_boosted_e, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_leadingJet_eta_e", self.ak4BJets[0].eta, SL_boosted_e, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_pt_e", self.ak4BJets[1].pt, SL_boosted_e, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_eta_e", self.ak4BJets[1].eta, SL_boosted_e, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("SL_resolved_DR_jets_e", op.deltaR(self.ak4BJets[0].p4, self.ak4BJets[1].p4), SL_boosted_e, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                # SL resolved electron final state plots
                Plot.make1D("SL_resolved_nJets_mu", op.rng_len(self.ak4BJets), SL_boosted_mu, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("SL_resolved_InvM_leadingJet_pt_mu", self.ak4BJets[0].pt, SL_boosted_mu, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_leadingJet_eta_mu", self.ak4BJets[0].eta, SL_boosted_mu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_pt_mu", self.ak4BJets[1].pt, SL_boosted_mu, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_eta_mu", self.ak4BJets[1].eta, SL_boosted_mu, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("SL_resolved_DR_jets_mu", op.deltaR(self.ak4BJets[0].p4, self.ak4BJets[1].p4), SL_boosted_mu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
            ])

        return plots
