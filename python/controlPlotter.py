
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
        self.tightMuons = defs.muonTightSel(self.muons)
        self.tightElectrons = defs.elTightSel(self.clElectrons)

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

        self.ak8BJets = op.select(self.ak8Jets, ak8Btag)

        # Ak4 Jet Collection cleaned from Ak8b #
        def cleanAk4FromAk8b(ak4j): return op.AND(op.rng_len(
            self.ak8BJets) > 0, op.deltaR(ak4j.p4, self.ak8BJets[0].p4) > 1.2)
        self.ak4JetsCleanedFromAk8b = op.select(self.ak4Jets, cleanAk4FromAk8b)
        
        def labeler(label):
            return {'labels': [{'text': label, 'position': [0.23, 0.87], 'size': 25}]}

        if self.channel == 'DL':

            DL_boosted_ee, DL_boosted_mumu,\
            DL_boosted_emu, DL_resolved_ee,\
            DL_resolved_mumu, DL_resolved_emu = makeDLSelection(self, noSel)

            # cutflow report
            yields.add(DL_boosted_ee, 'DL boosted ee')
            yields.add(DL_boosted_mumu, 'DL boosted mumu')
            yields.add(DL_boosted_emu, 'DL boosted emu')
            yields.add(DL_resolved_ee, 'DL resolved ee')
            yields.add(DL_resolved_mumu, 'DL resolved mumu')
            yields.add(DL_resolved_emu, 'DL resolved emu')
            
            # labels on plots
            DLboostedEE_label = labeler('DL boosted EE')
            DLboostedMuMu_label = labeler('DL boosted MuMu')
            DLboostedEMU_label = labeler('DL boosted EMu')
            
            DLresolvedEE_label = labeler('DL resolved EE')
            DLresolvedMuMu_label = labeler('DL resolved MuMu')
            DLresolvedEMu_label = labeler('DL resolved EMu')

        if self.channel == 'SL':
            SL_resolved, SL_resolved_e,\
            SL_resolved_mu, SL_boosted,\
            SL_boosted_e, SL_boosted_mu = makeSLSelection(self, noSel)

            yields.add(SL_boosted, 'SL boosted')
            yields.add(SL_boosted_e, 'SL boosted e')
            yields.add(SL_boosted_mu, 'SL boosted mu')
            yields.add(SL_resolved, 'SL resolved')
            yields.add(SL_resolved_e, 'SL resolved e')
            yields.add(SL_resolved_mu, 'SL resolved mu')

        #############################################################################
        #                                 Plots                                     #
        #############################################################################
        
        if self.channel == 'DL':
            plots.extend([
                
                #########################################
                ######                             ######
                ######       DL boosted plots      ######
                ######                             ######
                #########################################

                # number of ak8 b-jets
                Plot.make1D("DL_boosted_nfatJet_ee", op.rng_len(self.ak8Jets), DL_boosted_ee, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_nfatJet_mumu", op.rng_len(self.ak8Jets), DL_boosted_mumu, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_nfatJet_emu", op.rng_len(self.ak8Jets), DL_boosted_emu, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet", plotopts=DLboostedEMU_label),

                # fatjet pt
                Plot.make1D("DL_boosted_fatJet_pt_ee", self.ak8Jets[0].pt, DL_boosted_ee, EqBin(
                    100, 200, 800), title="pT(ak8jet)", xTitle="Fatjet p_{T} (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_fatJet_pt_mumu", self.ak8Jets[0].pt, DL_boosted_mumu, EqBin(
                    100, 200, 800), title="pT(ak8jet)", xTitle="Fatjet p_{T} (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_fatJet_pt_emu", self.ak8Jets[0].pt, DL_boosted_emu, EqBin(
                    100, 200, 800), title="pT(ak8jet)", xTitle="Fatjet p_{T} (GeV/c)", plotopts=DLboostedEMU_label),

                # subjet1 pt
                Plot.make1D("DL_boosted_subjet1_pt_ee", self.ak8Jets[0].subJet1.pt, DL_boosted_ee, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="First sub-jet p_{T} (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subjet1_pt_mumu", self.ak8Jets[0].subJet1.pt, DL_boosted_mumu, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="First sub-jet p_{T} (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subjet1_pt_emu", self.ak8Jets[0].subJet1.pt, DL_boosted_emu, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="First sub-jet p_{T} (GeV/c)", plotopts=DLboostedEMU_label),

                # subjet2 pt
                Plot.make1D("DL_boosted_subjet2_pt_ee", self.ak8Jets[0].subJet2.pt, DL_boosted_ee, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Second sub-jet p_{T} (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subjet2_pt_mumu", self.ak8Jets[0].subJet2.pt, DL_boosted_mumu, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Second sub-jet p_{T} (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subjet2_pt_emu", self.ak8Jets[0].subJet2.pt, DL_boosted_emu, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Second sub-jet p_{T} (GeV/c)", plotopts=DLboostedEMU_label),

                # fatjet eta
                Plot.make1D("DL_boosted_fatJet_eta_ee", self.ak8Jets[0].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_fatJet_eta_mumu", self.ak8Jets[0].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_fatJet_eta_emu", self.ak8Jets[0].eta, DL_boosted_emu, EqBin(
                    30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedEMU_label),

                # Invariant mass of leptons
                Plot.make1D("DL_boosted_InvM_ee", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_boosted_ee, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_InvM_mumu", op.invariant_mass(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_boosted_mumu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons  (GeV/c^{2})", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_InvM_emu", op.invariant_mass(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_boosted_emu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair  (GeV/c^{2})", plotopts=DLboostedEMU_label),

                # total charge of leptons
                Plot.make1D("DL_boosted_totalCharge_ee", op.sum(self.firstOSElEl[0].charge, self.firstOSElEl[1].charge), DL_boosted_ee, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_totalCharge_mumu", op.sum(self.firstOSMuMu[0].charge, self.firstOSMuMu[1].charge), DL_boosted_mumu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons ", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_totalCharge_emu", op.sum(self.firstOSElMu[0].charge, self.firstOSElMu[1].charge), DL_boosted_emu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair ", plotopts=DLboostedEMU_label),

                # invariant mass of subjets
                Plot.make1D("DL_boosted_InvM_jj_ee", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_ee, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_InvM_jj_mumu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_mumu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_InvM_jj_emu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_emu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedEMU_label),
                
                # leading lepton pt
                Plot.make1D("DL_boosted_leadingLepton_pt_ee", self.firstOSElEl[0].pt, DL_boosted_ee, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_leadingLepton_pt_mumu", self.firstOSMuMu[0].pt, DL_boosted_mumu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_leadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].pt, self.firstOSElMu[1].pt), DL_boosted_emu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedEMU_label),
                
                # sub-leading lepton pt
                Plot.make1D("DL_boosted_subleadingLepton_pt_ee", self.firstOSElEl[1].pt, DL_boosted_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subleadingLepton_pt_mumu", self.firstOSMuMu[1].pt, DL_boosted_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subleadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].pt, self.firstOSElMu[0].pt), DL_boosted_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedEMU_label),
                
                # leading lepton eta
                Plot.make1D("DL_boosted_leadingLepton_eta_ee", self.firstOSElEl[0].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_leadingLepton_eta_mumu", self.firstOSMuMu[0].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_leadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].eta, self.firstOSElMu[1].eta), DL_boosted_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedEMU_label),
                
                # sub-leading lepton eta
                Plot.make1D("DL_boosted_subleadingLepton_eta_ee", self.firstOSElEl[1].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subleadingLepton_eta_mumu", self.firstOSMuMu[1].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subleadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].eta, self.firstOSElMu[0].eta), DL_boosted_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedEMU_label),
                
                # DR between leading and sub-leading lepton
                Plot.make1D("DL_boosted_DR_leptons_ee", op.deltaR(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_boosted_ee, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_leptons_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_boosted_mumu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_leptons_emu", op.deltaR(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_boosted_emu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedEMU_label),

                # DR between leading lepton and ak8 jet
                Plot.make1D("DL_boosted_DR_leptonANDak8bjet_ee", op.deltaR(self.firstOSElEl[0].p4, self.ak8Jets[0].p4), DL_boosted_ee, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_leptonANDak8bjet_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.ak8Jets[0].p4), DL_boosted_mumu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_leptonANDak8bjet_emu", op.deltaR(op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), self.ak8Jets[0].p4), DL_boosted_emu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedEMU_label),
                
                # number of electrons
                Plot.make1D("DL_boosted_nElectrons_ee", op.rng_len(self.tightElectrons), DL_boosted_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_nElectrons_mumu", op.rng_len(self.tightElectrons), DL_boosted_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_nElectrons_emu", op.rng_len(self.tightElectrons), DL_boosted_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEMU_label),
                
                # number of muons
                Plot.make1D("DL_boosted_nMuons_ee", op.rng_len(self.tightMuons), DL_boosted_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_nMuons_mumu", op.rng_len(self.tightMuons), DL_boosted_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_nMuons_emu", op.rng_len(self.tightMuons), DL_boosted_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEMU_label),
                
                #########################################
                ######                             ######
                ######      DL resolved plots      ######
                ######                             ######
                #########################################
                
                # number of ak4 bjets
                Plot.make1D("DL_resolved_nAK4bJets_ee", op.rng_len(self.ak4BJets), DL_resolved_ee, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nAK4bJets_mumu", op.rng_len(self.ak4BJets), DL_resolved_mumu, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nAK4bJets_emu", op.rng_len(self.ak4BJets), DL_resolved_emu, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolvedEMu_label),
                
                # ak4 bjet pt
                Plot.make1D("DL_resolved_ak4BJet_pt_ee", self.ak4BJets[0].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_ak4BJet_pt_mumu", self.ak4BJets[0].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_ak4BJet_pt_emu", self.ak4BJets[0].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolvedEMu_label),

                # ak4 bjet eta
                Plot.make1D("DL_resolved_ak4BJet_eta_ee", self.ak4BJets[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_ak4BJet_eta_mumu", self.ak4BJets[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_ak4BJet_eta_emu", self.ak4BJets[0].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolvedEMu_label),

                # number of ak4 jets
                Plot.make1D("DL_resolved_nak4Jets_ee", op.rng_len(self.ak4Jets), DL_resolved_ee, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nak4Jets_mumu", op.rng_len(self.ak4Jets), DL_resolved_mumu, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nak4Jets_emu", op.rng_len(self.ak4Jets), DL_resolved_emu, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolvedEMu_label),
                
                # leading jet pt
                Plot.make1D("DL_resolved_leadingJet_pt_ee", self.ak4Jets[0].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingJet_pt_mumu", self.ak4Jets[0].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingJet_pt_emu", self.ak4Jets[0].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolvedEMu_label),
                
                # leading jet eta
                Plot.make1D("DL_resolved_leadingJet_eta_ee", self.ak4Jets[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingJet_eta_mumu", self.ak4Jets[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingJet_eta_emu", self.ak4Jets[0].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolvedEMu_label),
                
                # sub-leading jet pt
                Plot.make1D("DL_resolved_subleadingJet_pt_ee", self.ak4Jets[1].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingJet_pt_mumu", self.ak4Jets[1].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingJet_pt_emu", self.ak4Jets[1].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolvedEMu_label),

                # sub-leading jet eta
                Plot.make1D("DL_resolved_subleadingJet_eta_ee", self.ak4Jets[1].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingJet_eta_mumu", self.ak4Jets[1].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingJet_eta_emu", self.ak4Jets[1].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolvedEMu_label),
                
                # DR between leading and sub-leading jet
                Plot.make1D("DL_resolved_DR_jets_ee", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_jets_mumu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_jets_emu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolvedEMu_label),
                
                # Invariant mass of leptons
                Plot.make1D("DL_resolved_InvM_ee", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_resolved_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_InvM_mumu", op.invariant_mass(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_resolved_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_InvM_emu", op.invariant_mass(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_resolved_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (GeV/c^{2})", plotopts=DLresolvedEMu_label),

                # total charge of leptons
                Plot.make1D("DL_resolved_totalCharge_ee", op.sum(self.firstOSElEl[0].charge, self.firstOSElEl[1].charge), DL_resolved_ee, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_totalCharge_mumu", op.sum(self.firstOSMuMu[0].charge, self.firstOSMuMu[1].charge), DL_resolved_mumu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_totalCharge_emu", op.sum(self.firstOSElMu[0].charge, self.firstOSElMu[1].charge), DL_resolved_emu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair", plotopts=DLresolvedEMu_label),

                # leading lepton pt
                Plot.make1D("DL_resolved_leadingLepton_pt_ee", self.firstOSElEl[0].pt, DL_resolved_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingLepton_pt_mumu", self.firstOSMuMu[0].pt, DL_resolved_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].pt, self.firstOSElMu[1].pt), DL_resolved_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolvedEMu_label),
                
                # sub-leading lepton pt
                Plot.make1D("DL_resolved_subleadingLepton_pt_ee", self.firstOSElEl[1].pt, DL_resolved_ee, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingLepton_pt_mumu", self.firstOSMuMu[1].pt, DL_resolved_mumu, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].pt, self.firstOSElMu[0].pt), DL_resolved_emu, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolvedEMu_label),
                
                # leading lepton eta
                Plot.make1D("DL_resolved_leadingLepton_eta_ee", self.firstOSElEl[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingLepton_eta_mumu", self.firstOSMuMu[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].eta, self.firstOSElMu[1].eta), DL_resolved_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolvedEMu_label),
                
                # sub-leading lepton eta
                Plot.make1D("DL_resolved_subleadingLepton_eta_ee", self.firstOSElEl[1].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingLepton_eta_mumu", self.firstOSMuMu[1].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].eta, self.firstOSElMu[0].eta), DL_resolved_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolvedEMu_label),
                
                # DR between leading and sub-leading lepton
                Plot.make1D("DL_resolved_DR_leptons_ee", op.deltaR(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_leptons_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_leptons_emu", op.deltaR(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolvedEMu_label),

                # DR between leading lepton and ak4 b jet
                Plot.make1D("DL_resolved_DR_leptonANDak8bjet_ee", op.deltaR(self.firstOSElEl[0].p4, self.ak4BJets[0].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_leptonANDak8bjet_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.ak4BJets[0].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_leptonANDak8bjet_emu", op.deltaR(op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), self.ak4BJets[0].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolvedEMu_label),
                
                # number of electrons
                Plot.make1D("DL_resolved_nElectrons_ee", op.rng_len(self.tightElectrons), DL_resolved_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nElectrons_mumu", op.rng_len(self.tightElectrons), DL_resolved_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nElectrons_emu", op.rng_len(self.tightElectrons), DL_resolved_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEMu_label),
                
                # number of muons
                Plot.make1D("DL_resolved_nMuons_ee", op.rng_len(self.tightMuons), DL_resolved_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nMuons_mumu", op.rng_len(self.tightMuons), DL_resolved_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nMuons_emu", op.rng_len(self.tightMuons), DL_resolved_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEMu_label),
            ])
        if self.channel == "SL":
            plots.extend([
                # SL boosted plots
                Plot.make1D("SL_boosted_fatJet_pt", self.ak8BJets[0].pt, SL_boosted, EqBin(
                    400, 200, 1000), title="pT(j)", xTitle="pT(j) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet1_pt", self.ak8BJets[0].subJet1.pt, SL_boosted, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="pT(subjet1) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet2_pt", self.ak8BJets[0].subJet2.pt, SL_boosted, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="pT(subjet2) (GeV/c)"),
                Plot.make1D("SL_boosted_fatJet_eta", self.ak8BJets[0].eta, SL_boosted, EqBin(
                    30, -3, 3), title="eta(j)", xTitle="eta(j)"),
                Plot.make1D("SL_boosted_subjet1_eta", self.ak8BJets[0].subJet1.eta, SL_boosted, EqBin(
                    30, -3, 3), title="eta(subjet1)", xTitle="eta(subjet1)"),
                Plot.make1D("SL_boosted_subjet2_eta", self.ak8BJets[0].subJet2.eta, SL_boosted, EqBin(
                    30, -3, 3), title="eta(subjet2)", xTitle="eta(subjet2)"),
                Plot.make1D("SL_boosted_InvM_jj", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), SL_boosted, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet1_eta", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].eta], SL_boosted, [
                    EqBin(100, 0., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet1 eta", xTitle="Invariant Mass of jets (GeV/c^{2})", yTitle="eta(j1)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet2_eta", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].subJet2.eta], SL_boosted, [
                    EqBin(100, 0., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet2 eta", xTitle="Invariant Mass of jets (GeV/c^{2})", yTitle="eta(j2)"),

                # SL boosted electron final state plots
                Plot.make1D("SL_boosted_fatJet_pt_e", self.ak8BJets[0].pt, SL_boosted_e, EqBin(
                    400, 200, 1000), title="pT(j)", xTitle="pT(j) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet1_pt_e", self.ak8BJets[0].subJet1.pt, SL_boosted_e, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="pT(subjet1) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet2_pt_e", self.ak8BJets[0].subJet2.pt, SL_boosted_e, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="pT(subjet2) (GeV/c)"),
                Plot.make1D("SL_boosted_fatJet_eta_e", self.ak8BJets[0].eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(j)", xTitle="eta(j)"),
                Plot.make1D("SL_boosted_subjet1_eta_e", self.ak8BJets[0].subJet1.eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(subjet1)", xTitle="eta(subjet1)"),
                Plot.make1D("SL_boosted_subjet2_eta_e", self.ak8BJets[0].subJet2.eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(subjet2)", xTitle="eta(subjet2)"),
                Plot.make1D("SL_boosted_InvM_jj_e", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), SL_boosted_e, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet1_eta_e", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].eta], SL_boosted_e, [
                    EqBin(100, 0., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet1 eta", xTitle="Invariant Mass of jets (GeV/c^{2})", yTitle="eta(j1)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet2_eta_e", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].subJet2.eta], SL_boosted_e, [
                    EqBin(100, 0., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet2 eta", xTitle="Invariant Mass of jets (GeV/c^{2})", yTitle="eta(j2)"),

                # SL boosted muon final state plots
                Plot.make1D("SL_boosted_fatJet_pt_mu", self.ak8BJets[0].pt, SL_boosted_mu, EqBin(
                    400, 200, 1000), title="pT(j)", xTitle="pT(j) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet1_pt_mu", self.ak8BJets[0].subJet1.pt, SL_boosted_mu, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="pT(subjet1) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet2_pt_mu", self.ak8BJets[0].subJet2.pt, SL_boosted_mu, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="pT(subjet2) (GeV/c)"),
                Plot.make1D("SL_boosted_fatJet_eta_mu", self.ak8BJets[0].eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(j)", xTitle="eta(j)"),
                Plot.make1D("SL_boosted_subjet1_eta_mu", self.ak8BJets[0].subJet1.eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(subjet1)", xTitle="eta(subjet1)"),
                Plot.make1D("SL_boosted_subjet2_eta_mu", self.ak8BJets[0].subJet2.eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(subjet2)", xTitle="eta(subjet2)"),
                Plot.make1D("SL_boosted_InvM_jj_mu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), SL_boosted_mu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet1_eta_mu", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].eta], SL_boosted_mu, [
                    EqBin(100, 0., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet1 eta", xTitle="Invariant Mass of jets (GeV/c^{2})", yTitle="eta(j1)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet2_eta_mu", [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), self.ak8BJets[0].subJet2.eta], SL_boosted_mu, [
                    EqBin(100, 0., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet2 eta", xTitle="Invariant Mass of jets (GeV/c^{2})", yTitle="eta(j2)"),

                # SL resolved plots
                Plot.make1D("SL_resolved_nJets", op.rng_len(self.ak4BJets), SL_resolved, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("SL_resolved_InvM_leadingJet_pt", self.ak4BJets[0].pt, SL_resolved, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_leadingJet_eta", self.ak4BJets[0].eta, SL_resolved, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_pt", self.ak4BJets[1].pt, SL_resolved, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_eta", self.ak4BJets[1].eta, SL_resolved, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("SL_resolved_DR_jets", op.deltaR(self.ak4BJets[0].p4, self.ak4BJets[1].p4), SL_resolved, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                # SL resolved electron final state plots
                Plot.make1D("SL_resolved_nJets_e", op.rng_len(self.ak4BJets), SL_boosted_e, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("SL_resolved_InvM_leadingJet_pt_e", self.ak4BJets[0].pt, SL_boosted_e, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_leadingJet_eta_e", self.ak4BJets[0].eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_pt_e", self.ak4BJets[1].pt, SL_boosted_e, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_eta_e", self.ak4BJets[1].eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("SL_resolved_DR_jets_e", op.deltaR(self.ak4BJets[0].p4, self.ak4BJets[1].p4), SL_boosted_e, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                # SL resolved electron final state plots
                Plot.make1D("SL_resolved_nJets_mu", op.rng_len(self.ak4BJets), SL_boosted_mu, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("SL_resolved_InvM_leadingJet_pt_mu", self.ak4BJets[0].pt, SL_boosted_mu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_leadingJet_eta_mu", self.ak4BJets[0].eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_pt_mu", self.ak4BJets[1].pt, SL_boosted_mu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_subleadingJet_eta_mu", self.ak4BJets[1].eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("SL_resolved_DR_jets_mu", op.deltaR(self.ak4BJets[0].p4, self.ak4BJets[1].p4), SL_boosted_mu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
            ])

        return plots
