
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

        # lepton cone-pt definitions
        muon_conept = defs.muonConePt(tree.Muon)
        electron_conept = defs.elConePt(tree.Electron)

        # lepton definitions sorted by their cone-pt
        muons = op.sort(defs.muonDef(tree.Muon), lambda mu: -muon_conept[mu.idx])
        electrons = op.sort(defs.elDef(tree.Electron), lambda el: -electron_conept[el.idx])

        # cleaning electrons wrt muons
        clElectrons = defs.cleanElectrons(electrons, muons)

        # Fakeable leptons
        fakeMuons = defs.muonFakeSel(muons)
        fakeElectrons = defs.elFakeSel(clElectrons)

        # tight leptons
        tightMuons = defs.muonTightSel(fakeMuons)
        tightElectrons = defs.elTightSel(fakeElectrons)

        # Dilepton selections
        if self.channel == "DL":
            def leptonOS(l1, l2): return l1.charge != l2.charge

            ElElDileptonPreSel = op.combine(electrons, N=2)
            MuMuDileptonPreSel = op.combine(muons, N=2)
            ElMuDileptonPreSel = op.combine((electrons, muons))

            OSElElDileptonPreSel = op.combine(electrons, N=2, pred=leptonOS)
            OSMuMuDileptonPreSel = op.combine(muons, N=2, pred=leptonOS)
            OSElMuDileptonPreSel = op.combine((electrons, muons), pred=leptonOS)

        # Dilepton for selection #
            ElElFakeSel = op.combine(fakeElectrons, N=2)
            MuMuFakeSel = op.combine(fakeMuons, N=2)
            ElMuFakeSel = op.combine((fakeElectrons, fakeMuons))

            ElElTightSel = op.combine(tightElectrons, N=2)
            MuMuTightSel = op.combine(tightMuons, N=2)
            ElMuTightSel = op.combine((tightElectrons, tightMuons))

            if self.is_MC:
                def is_matched(lep): return op.OR(lep.genPartFlav == 1,  # Prompt muon or electron
                                                  lep.genPartFlav == 15)  # From tau decay
                # lep.genPartFlav==22) # From photon conversion (only available for electrons)
            else:
                def is_matched(lep): return op.c_bool(True)

            def dilepton_matched(dilep): return op.AND(
                is_matched(dilep[0]), is_matched(dilep[1]))

            def electronTightSel(ele): return ele.mvaTTH >= 0.30

            def muonTightSel(mu): return op.AND(mu.mvaTTH >= 0.50, mu.mediumId)

            def tightpair_ElEl(dilep): return op.AND(
                dilepton_matched(dilep),
                electronTightSel(dilep[0]),
                electronTightSel(dilep[1])
            )

            def tightpair_MuMu(dilep): return op.AND(
                dilepton_matched(dilep),
                muonTightSel(dilep[0]),
                muonTightSel(dilep[1])
            )

            def tightpair_ElMu(dilep): return op.AND(
                dilepton_matched(dilep),
                electronTightSel(dilep[0]),
                muonTightSel(dilep[1])
            )
        # Taus
        taus = defs.tauDef(tree.Tau)
        cleanedTaus = defs.cleanTaus(taus, fakeElectrons, fakeMuons)

        # AK4 Jets sorted by their pt
        ak4JetsPreSel = op.sort(defs.ak4jetDef(tree.Jet), lambda jet: -jet.pt)

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

        ak4Jets = op.select(ak4JetsPreSel, cleanAk4Jets)
        ak4JetsByBtagScore = op.sort(ak4Jets, lambda j: -j.btagDeepFlavB)

        # bTagging for ak4 jets
        def ak4BtagLooseSel(jet): return jet.btagDeepFlavB > 0.0494
        def ak4BtagSel(jet): return jet.btagDeepFlavB > 0.2770
        def ak4NoBtagSel(jet): return jet.btagDeepFlavB <= 0.2770

        ak4BJets = op.select(ak4Jets, ak4BtagSel)
        ak4BJetsLoose = op.select(ak4Jets, ak4BtagLooseSel)
        ak4LightJetsByPt = op.select(ak4Jets, ak4NoBtagSel)
        ak4LightJetsByBtagScore = op.sort(
            ak4LightJetsByPt, lambda jet: -jet.btagDeepFlavB)
        remainingJets = op.select(
            ak4LightJetsByPt, lambda jet: jet.idx != ak4LightJetsByBtagScore[0].idx)

        def makeJetPairs(jets): return op.combine(
            jets, N=2, pred=lambda j1, j2: j1.pt > j2.pt, samePred=lambda j1, j2: j1.idx != j2.idx)
        # --------------------------------------------- #
        bJetsByScore = ak4JetsByBtagScore[:op.min(op.rng_len(
            ak4JetsByBtagScore), op.static_cast("std::size_t", op.c_int(2)))]
        probableWJets = op.select(ak4Jets, lambda jet: op.NOT(
            op.rng_any(bJetsByScore, lambda bjet: jet.idx == bjet.idx)))
        wJetsByPt = probableWJets[:op.min(op.rng_len(
            probableWJets), op.static_cast("std::size_t", op.c_int(2)))]

        def passWMassCutSel(wjets): return op.switch(op.rng_len(wjets) == 2, op.abs(
            op.invariant_mass(wjets[0].p4, wjets[1].p4)-80.4) < op.c_float(15.0), op.c_bool(False))

        # AK8 Jets
        ak8Jets = defs.ak8jetDef(tree.FatJet)

        if self.channel == 'SL': # sorted by btag score
            ak8JetsPreSel = op.sort(ak8Jets, lambda j: -j.btagDeepB)
        if self.channel == 'DL': # sorted by pt
            ak8JetsPreSel = op.sort(ak8Jets, lambda j: -j.pt)

        # cleaning ak8 jets wrt to leptons
        if self.channel == 'SL':
            cleanAk8Jets = cleaningWithRespectToLeadingLepton(0.8)
        if self.channel == 'DL':
            cleanAk8Jets = cleaningWithRespectToLeadingLeptons(0.8)

        ak8Jets = op.select(ak8JetsPreSel, cleanAk8Jets)

        # 2018 DeepJet WP
        def subjetBtag(subjet): return subjet.btagDeepB > 0.4184

        def ak8Btag(fatjet): return op.OR(op.AND(fatjet.subJet1.pt >= 30, subjetBtag(fatjet.subJet1)),
                                          op.AND(fatjet.subJet2.pt >= 30, subjetBtag(fatjet.subJet2)))

        def ak8noBtag(fatjet): return op.NOT(op.OR(op.AND(fatjet.subJet1.pt >= 30, subjetBtag(fatjet.subJet1)),
                                                   op.AND(fatjet.subJet2.pt >= 30, subjetBtag(fatjet.subJet2))))

        def ak8Btag_bothSubJets(fatjet): return op.AND(op.AND(fatjet.subJet1.pt >= 30, subjetBtag(fatjet.subJet1)),
                                                       op.AND(fatjet.subJet2.pt >= 30, subjetBtag(fatjet.subJet2)))

        ak8BJets = op.select(ak8Jets, ak8Btag)
        ak8nonBJets = op.select(ak8Jets, ak8noBtag)
        # Ak4 Jet Collection cleaned from Ak8b #

        def cleanAk4FromAk8b(ak4j): return op.AND(op.rng_len(
            ak8BJets) > 0, op.deltaR(ak4j.p4, ak8BJets[0].p4) > 1.2)
        ak4JetsCleanedFromAk8b = op.select(ak4Jets, cleanAk4FromAk8b)

        # used as a BDT input for SemiBoosted category
        def btaggedSubJets(fjet): return op.switch(
            ak8Btag_bothSubJets(fjet), op.c_float(2.0), op.c_float(1.0))
        nMediumBTaggedSubJets = op.rng_sum(ak8BJets, btaggedSubJets)

        # common variables for DL and SL channels
        Zmass = 91.1876

        def lowMllCut(dileptons): return op.NOT(op.rng_any(
            dileptons, lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4) < 12.))

        def outZ(dileptons): return op.NOT(op.rng_any(
            dileptons, lambda dilep: op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4) - Zmass) < 10.))

        # end of common variables

        ### Di-leptonic channel ###
        if self.channel == 'DL':
            def OSDilepton(dilep): return dilep[0].charge != dilep[1].charge

            # Pt cuts #
            def lowPtCutElEl(dilep): return op.AND(
                electron_conept[dilep[0].idx] > 15, electron_conept[dilep[1].idx] > 15)  # subleading above 15 GeV

            def lowPtCutMuMu(dilep): return op.AND(
                muon_conept[dilep[0].idx] > 15, muon_conept[dilep[1].idx] > 15)  # subleading above 15 GeV

            def lowPtCutElMu(dilep): return op.AND(
                electron_conept[dilep[0].idx] > 15, muon_conept[dilep[1].idx] > 15)  # subleading above 15 GeV

            def leadingPtCutElEl(dilep): return op.OR(
                electron_conept[dilep[0].idx] > 25, electron_conept[dilep[1].idx] > 25)  # leading above 25 GeV

            def leadingPtCutMuMu(dilep): return op.OR(
                muon_conept[dilep[0].idx] > 25, muon_conept[dilep[1].idx] > 25)  # leading above 25 GeV

            def leadingPtCutElMu(dilep): return op.OR(
                electron_conept[dilep[0].idx] > 25, muon_conept[dilep[1].idx] > 25)  # leading above 25 GeV

            # Fake lepton selection #
            # elelSel = noSel.refine('elelSel', cut=[op.rng_len(ElElFakeSel) >= 1,
            #                                        op.OR(op.rng_len(fakeMuons) == 0,
            #                                              op.AND(op.rng_len(fakeMuons) == 1,
            #                                                     electron_conept[ElElFakeSel[0]
            #                                                                     [0].idx] > muon_conept[fakeMuons[0].idx],
            #                                                     electron_conept[ElElFakeSel[0][1].idx] > muon_conept[fakeMuons[0].idx]),
            #                                              op.AND(op.rng_len(fakeMuons) >= 2,
            #                                                     electron_conept[ElElFakeSel[0]
            #                                                                     [0].idx] > muon_conept[fakeMuons[0].idx],
            #                                                     electron_conept[ElElFakeSel[0]
            #                                                                     [1].idx] > muon_conept[fakeMuons[0].idx],
            #                                                     electron_conept[ElElFakeSel[0]
            #                                                                     [0].idx] > muon_conept[fakeMuons[1].idx],
            #                                                     electron_conept[ElElFakeSel[0][1].idx] > muon_conept[fakeMuons[1].idx]))])

            # mumuSel = noSel.refine('mumuSel', cut=[op.rng_len(MuMuFakeSel) >= 1,
            #                                        op.OR(op.rng_len(fakeElectrons) == 0,
            #                                              op.AND(op.rng_len(fakeElectrons) == 1,
            #                                                     muon_conept[MuMuFakeSel[0]
            #                                                                 [0].idx] > electron_conept[fakeElectrons[0].idx],
            #                                                     muon_conept[MuMuFakeSel[0][1].idx] > electron_conept[fakeElectrons[0].idx]),
            #                                              op.AND(op.rng_len(fakeElectrons) >= 2,
            #                                                     muon_conept[MuMuFakeSel[0]
            #                                                                 [0].idx] > electron_conept[fakeElectrons[0].idx],
            #                                                     muon_conept[MuMuFakeSel[0]
            #                                                                 [1].idx] > electron_conept[fakeElectrons[0].idx],
            #                                                     muon_conept[MuMuFakeSel[0]
            #                                                                 [0].idx] > electron_conept[fakeElectrons[1].idx],
            #                                                     muon_conept[MuMuFakeSel[0][1].idx] > electron_conept[fakeElectrons[1].idx]))])
            # elmuSel = noSel.refine('elmuSel', cut=[op.rng_len(ElMuFakeSel) >= 1,
            #                                        op.OR(op.AND(op.rng_len(fakeElectrons) == 1,
            #                                                     op.rng_len(fakeMuons) == 1),
            #                                              op.AND(op.rng_len(fakeElectrons) >= 2,
            #                                                     op.rng_len(
            #                                                  fakeMuons) == 1,
            #                                            muon_conept[ElMuFakeSel[0][1].idx] > electron_conept[fakeElectrons[1].idx]),
            #     op.AND(op.rng_len(fakeMuons) >= 2,
            #                                            op.rng_len(
            #         fakeElectrons) == 1,
            #                                            electron_conept[ElMuFakeSel[0][0].idx] > muon_conept[fakeMuons[1].idx]),
            #     op.AND(op.rng_len(fakeElectrons) >= 2,
            #                                            op.rng_len(
            #         fakeMuons) >= 2,
            #                                            muon_conept[ElMuFakeSel[0]
            #                                                        [1].idx] > electron_conept[fakeElectrons[1].idx],
            #                                            electron_conept[ElMuFakeSel[0][0].idx] > muon_conept[fakeMuons[1].idx]))])

            # OS cut
            elelSel = noSel.refine('OSelelSel', cut=[OSDilepton(ElElFakeSel[0])])
            mumuSel = noSel.refine('OSmumuSel', cut=[OSDilepton(MuMuFakeSel[0])])
            elmuSel = noSel.refine('OSelmuSel', cut=[OSDilepton(ElMuFakeSel[0])])

            # Pt cuts
            elelSel.refine('elelptSel', cut=[lowPtCutElEl(ElElFakeSel[0]),
                           leadingPtCutElEl(ElElFakeSel[0])])
            mumuSel.refine('mumuptSel', cut=[lowPtCutMuMu(MuMuFakeSel[0]),
                           leadingPtCutMuMu(MuMuFakeSel[0])])
            elmuSel.refine('elmuptSel', cut=[lowPtCutElMu(ElMuFakeSel[0]),
                           leadingPtCutElMu(ElMuFakeSel[0])])

            # Mll cut
            mllCut = [lowMllCut(ElElDileptonPreSel), lowMllCut(
                MuMuDileptonPreSel), lowMllCut(ElMuDileptonPreSel)]
            elelSel.refine('elelMllSel', cut=mllCut)
            mumuSel.refine('mumuMllSel', cut=mllCut)
            elmuSel.refine('elmuMllSel', cut=mllCut)

            # Z-veto
            outZCut = [outZ(OSElElDileptonPreSel), outZ(OSMuMuDileptonPreSel)]
            elelSel.refine('OSoutZelelSel', cut=outZCut)
            mumuSel.refine('OSoutZmumuSel', cut=outZCut)
            elmuSel.refine('OSoutZelmuSel', cut=outZCut)

            # di-lepton cut
            elelSel.refine('dileptonCut_ee', cut=[
                tightpair_ElEl(ElElFakeSel[0]),
                op.rng_len(tightElectrons) == 2,
                op.rng_len(tightMuons) == 0,
                ElElTightSel[0][0].idx == ElElFakeSel[0][0].idx,
                ElElTightSel[0][1].idx == ElElFakeSel[0][1].idx]
            )
            mumuSel.refine('dileptonCut_mumu', cut=[
                tightpair_MuMu(MuMuFakeSel[0]),
                op.rng_len(tightMuons) == 2,
                op.rng_len(tightElectrons) == 0,
                MuMuTightSel[0][0].idx == MuMuFakeSel[0][0].idx,
                MuMuTightSel[0][1].idx == MuMuFakeSel[0][1].idx]
            )
            elmuSel.refine('dileptonCut_emu', cut=[
                tightpair_ElMu(ElMuFakeSel[0]),
                op.rng_len(tightElectrons) == 1,
                op.rng_len(tightMuons) == 1,
                ElMuTightSel[0][0].idx == ElMuFakeSel[0][0].idx,
                ElMuTightSel[0][1].idx == ElMuFakeSel[0][1].idx]
            )

            # boosted -> and at least one b-tagged ak8 jet
            DL_boosted_ee = elelSel.refine(
                'DL_boosted_ee', cut=(op.rng_len(ak8BJets) >= 1))
            DL_boosted_mumu = mumuSel.refine(
                'DL_boosted_mumu', cut=(op.rng_len(ak8BJets) >= 1))
            DL_boosted_emu = elmuSel.refine(
                'DL_boosted_emu', cut=(op.rng_len(ak8BJets) >= 1))

            # resolved -> and at least two ak4 jets with at least one b-tagged and no ak8 jets
            DL_resolved_1b_ee = elelSel.refine('DL_resolved_1b_ee', cut=(op.AND(op.rng_len(ak4BJets) >= 1, op.rng_len(ak8Jets) == 0)))
            DL_resolved_1b_mumu = mumuSel.refine('DL_resolved_1b_mumu', cut=(op.AND(op.rng_len(ak4BJets) >= 1, op.rng_len(ak8Jets) == 0)))
            DL_resolved_1b_emu = elmuSel.refine('DL_resolved_1b_emu', cut=(op.AND(op.rng_len(ak4BJets) >= 1, op.rng_len(ak8Jets) == 0)))

            DL_resolved_2b_ee = DL_resolved_1b_ee.refine('DL_resolved_ee', cut=(op.rng_len(ak4BJets) >= 2))
            DL_resolved_2b_mumu = DL_resolved_1b_mumu.refine('DL_resolved_mumu', cut=(op.rng_len(ak4BJets) >= 2))
            DL_resolved_2b_emu = DL_resolved_1b_emu.refine('DL_resolved_emu', cut=(op.rng_len(ak4BJets) >= 2))

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
            def elPtCut(lep): return electron_conept[lep[0].idx] > 32.0

            def muPtCut(lep): return muon_conept[lep[0].idx] > 25.0

            def tau_h_veto(taus): return op.rng_len(taus) == 0
            
            def leptonOS(l1, l2): return l1.charge != l2.charge

            ElElDileptonPreSel = op.combine(clElectrons, N=2)
            MuMuDileptonPreSel = op.combine(muons, N=2)
            ElMuDileptonPreSel = op.combine((clElectrons, muons))

            OSElElDileptonPreSel = op.combine(clElectrons, N=2, pred=leptonOS)
            OSMuMuDileptonPreSel = op.combine(muons, N=2, pred=leptonOS)
            
            outZcut = [outZ(OSElElDileptonPreSel), outZ(OSMuMuDileptonPreSel)]

            mllCut = [lowMllCut(ElElDileptonPreSel), lowMllCut(
                MuMuDileptonPreSel), lowMllCut(ElMuDileptonPreSel)]

            SL_resolved = noSel.refine('SL_resolved', cut=[
                mllCut, outZcut, tau_h_veto(cleanedTaus),
                op.rng_len(ak4Jets) >= 3,
                op.rng_len(ak4BJets) >= 1,
                op.rng_len(ak8BJets) == 0])

            SL_resolved_e = SL_resolved.refine('SL_resolved_e', cut=[
                elPtCut(tightElectrons),
                op.rng_len(tightElectrons) == 1,
                op.rng_len(tightMuons) == 0])

            SL_resolved_mu = SL_resolved.refine('SL_resolved_mu', cut=[
                muPtCut(tightMuons),
                op.rng_len(tightElectrons) == 0,
                op.rng_len(tightMuons) == 1])

            SL_boosted = noSel.refine('SL_boosted', cut=[
                op.OR(elPtCut, muPtCut), lowMllCut, outZ, tau_h_veto,
                op.rng_len(ak8BJets) >= 1,
                op.rng_len(ak4JetsCleanedFromAk8b) >= 1])

            yields.add(SL_boosted, 'SL boosted')
            yields.add(SL_resolved, 'SL resolved')

        #############################################################################
        #                                 Plots                                     #
        #############################################################################
        if self.channel == 'DL':
            plots.extend([
                # DL boosted plots
                Plot.make1D("DL_boosted_fatJet_pt_ee", ak8Jets[0].pt, DL_boosted_ee, EqBin(
                    400, 200, 1000), title="pT(j1)", xTitle="pT(ak8jet) (GeV/c)"),
                Plot.make1D("DL_boosted_fatJet_pt_mumu", ak8Jets[0].pt, DL_boosted_mumu, EqBin(
                    400, 200, 1000), title="pT(j1)", xTitle="pT(ak8jet) (GeV/c)"),
                Plot.make1D("DL_boosted_fatJet_pt_emu", ak8Jets[0].pt, DL_boosted_emu, EqBin(
                    400, 200, 1000), title="pT(j1)", xTitle="pT(ak8jet) (GeV/c)"),

                Plot.make1D("DL_boosted_subjet1_pt_ee", ak8Jets[0].subJet1.pt, DL_boosted_ee, EqBin(
                    250, 0, 500), title=" pT(j1 subjet1)", xTitle="pT(j1 subjet1) (GeV/c)"),
                Plot.make1D("DL_boosted_subjet1_pt_mumu", ak8Jets[0].subJet1.pt, DL_boosted_mumu, EqBin(
                    250, 0, 500), title=" pT(j1 subjet1)", xTitle="pT(j1 subjet1) (GeV/c)"),
                Plot.make1D("DL_boosted_subjet1_pt_emu", ak8Jets[0].subJet1.pt, DL_boosted_emu, EqBin(
                    250, 0, 500), title=" pT(j1 subjet1)", xTitle="pT(j1 subjet1) (GeV/c)"),
                Plot.make1D("DL_boosted_subjet2_pt_ee", ak8Jets[0].subJet2.pt, DL_boosted_ee, EqBin(
                    250, 0, 500), title=" pT(j1 subjet2)", xTitle="pT(j1 subjet2) (GeV/c)"),
                Plot.make1D("DL_boosted_subjet2_pt_mumu", ak8Jets[0].subJet2.pt, DL_boosted_mumu, EqBin(
                    250, 0, 500), title=" pT(j1 subjet2)", xTitle="pT(j1 subjet2) (GeV/c)"),
                Plot.make1D("DL_boosted_subjet2_pt_emu", ak8Jets[0].subJet2.pt, DL_boosted_emu, EqBin(
                    250, 0, 500), title=" pT(j1 subjet2)", xTitle="pT(j1 subjet2) (GeV/c)"),
                Plot.make1D("DL_boosted_fatJet_eta_ee", ak8Jets[0].eta, DL_boosted_ee, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(ak8jet)"),
                Plot.make1D("DL_boosted_fatJet_eta_mumu", ak8Jets[0].eta, DL_boosted_mumu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(ak8jet)"),
                Plot.make1D("DL_boosted_fatJet_eta_emu", ak8Jets[0].eta, DL_boosted_emu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(ak8jet)"),
                Plot.make1D("DL_boosted_InvM_ee", op.invariant_mass(ElElDileptonPreSel[0][0].p4, ElElDileptonPreSel[0][1].p4), DL_boosted_ee, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electrons (boosted) (GeV/c^2)"),
                Plot.make1D("DL_boosted_InvM_mumu", op.invariant_mass(MuMuDileptonPreSel[0][0].p4, MuMuDileptonPreSel[0][1].p4), DL_boosted_mumu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muons (boosted) (GeV/c^2)"),
                Plot.make1D("DL_boosted_InvM_emu", op.invariant_mass(ElMuDileptonPreSel[0][0].p4, ElMuDileptonPreSel[0][1].p4), DL_boosted_emu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (boosted) (GeV/c^2)"),
                Plot.make1D("DL_boosted_InvM_jj_ee", op.invariant_mass(ak8Jets[0].subJet1.p4, ak8Jets[0].subJet2.p4), DL_boosted_ee, EqBin(
                    160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),
                Plot.make1D("DL_boosted_InvM_jj_mumu", op.invariant_mass(ak8Jets[0].subJet1.p4, ak8Jets[0].subJet2.p4), DL_boosted_mumu, EqBin(
                    160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),
                Plot.make1D("DL_boosted_InvM_jj_emu", op.invariant_mass(ak8Jets[0].subJet1.p4, ak8Jets[0].subJet2.p4), DL_boosted_emu, EqBin(
                    160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),

                # DL resolved 1b plots
                Plot.make1D("DL_resolved_1b_ee_nJets", op.rng_len(ak4Jets), DL_resolved_1b_ee, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_1b_ee_InvM_leadingJet_pt", ak4Jets[0].pt, DL_resolved_1b_ee, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_ee_InvM_subleadingJet_pt", ak4Jets[1].pt, DL_resolved_1b_ee, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_ee_InvM_leadingJet_eta", ak4Jets[0].eta, DL_resolved_1b_ee, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_1b_ee_InvM_subleadingJet_eta", ak4Jets[1].eta, DL_resolved_1b_ee, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_1b_ee_DR_jets", op.deltaR(ak4Jets[0].p4, ak4Jets[1].p4), DL_resolved_1b_ee, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_1b_ee_InvM_ee", op.invariant_mass(ElElDileptonPreSel[0][0].p4, ElElDileptonPreSel[0][1].p4), DL_resolved_1b_ee, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electrons (resolved) (GeV/c^2)"),

                Plot.make1D("DL_resolved_1b_mumu_nJets", op.rng_len(ak4Jets), DL_resolved_1b_mumu, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_1b_mumu_InvM_leadingJet_pt", ak4Jets[0].pt, DL_resolved_1b_mumu, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_mumu_InvM_subleadingJet_pt", ak4Jets[1].pt, DL_resolved_1b_mumu, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_mumu_InvM_leadingJet_eta", ak4Jets[0].eta, DL_resolved_1b_mumu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_1b_mumu_InvM_subleadingJet_eta", ak4Jets[1].eta, DL_resolved_1b_mumu, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_1b_mumu_DR_jets", op.deltaR(ak4Jets[0].p4, ak4Jets[1].p4), DL_resolved_1b_mumu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_1b_mumu_InvM_mumu", op.invariant_mass(MuMuDileptonPreSel[0][0].p4, MuMuDileptonPreSel[0][1].p4), DL_resolved_1b_mumu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muon (resolved) (GeV/c^2)"),

                Plot.make1D("DL_resolved_1b_emu_nJets", op.rng_len(ak4Jets), DL_resolved_1b_emu, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_1b_emu_InvM_leadingJet_pt", ak4Jets[0].pt, DL_resolved_1b_emu, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_emu_InvM_subleadingJet_pt", ak4Jets[1].pt, DL_resolved_1b_emu, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_1b_emu_InvM_leadingJet_eta", ak4Jets[0].eta, DL_resolved_1b_emu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_1b_emu_InvM_subleadingJet_eta", ak4Jets[1].eta, DL_resolved_1b_emu, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_1b_emu_DR_jets", op.deltaR(ak4Jets[0].p4, ak4Jets[1].p4), DL_resolved_1b_emu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_1b_emu_InvM_emu", op.invariant_mass(ElMuDileptonPreSel[0][0].p4, ElMuDileptonPreSel[0][1].p4), DL_resolved_1b_emu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muon (resolved) (GeV/c^2)"),


                # DL resolved 2b plots
                Plot.make1D("DL_resolved_2b_ee_nJets", op.rng_len(ak4Jets), DL_resolved_2b_ee, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_2b_ee_InvM_leadingJet_pt", ak4Jets[0].pt, DL_resolved_2b_ee, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_ee_InvM_subleadingJet_pt", ak4Jets[1].pt, DL_resolved_2b_ee, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_ee_InvM_leadingJet_eta", ak4Jets[0].eta, DL_resolved_2b_ee, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_2b_ee_InvM_subleadingJet_eta", ak4Jets[1].eta, DL_resolved_2b_ee, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_2b_ee_DR_jets", op.deltaR(ak4Jets[0].p4, ak4Jets[1].p4), DL_resolved_2b_ee, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_2b_ee_InvM_ee", op.invariant_mass(ElElDileptonPreSel[0][0].p4, ElElDileptonPreSel[0][1].p4), DL_resolved_2b_ee, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electrons (resolved) (GeV/c^2)"),

                Plot.make1D("DL_resolved_2b_mumu_nJets", op.rng_len(ak4Jets), DL_resolved_2b_mumu, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_2b_mumu_InvM_leadingJet_pt", ak4Jets[0].pt, DL_resolved_2b_mumu, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_mumu_InvM_subleadingJet_pt", ak4Jets[1].pt, DL_resolved_2b_mumu, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_mumu_InvM_leadingJet_eta", ak4Jets[0].eta, DL_resolved_2b_mumu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_2b_mumu_InvM_subleadingJet_eta", ak4Jets[1].eta, DL_resolved_2b_mumu, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_2b_mumu_DR_jets", op.deltaR(ak4Jets[0].p4, ak4Jets[1].p4), DL_resolved_2b_mumu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_2b_mumu_InvM_mumu", op.invariant_mass(MuMuDileptonPreSel[0][0].p4, MuMuDileptonPreSel[0][1].p4), DL_resolved_2b_mumu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of muons (resolved) (GeV/c^2)"),

                Plot.make1D("DL_resolved_2b_emu_nJets", op.rng_len(ak4Jets), DL_resolved_2b_emu, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("DL_resolved_2b_emu_InvM_leadingJet_pt", ak4Jets[0].pt, DL_resolved_2b_emu, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_emu_InvM_subleadingJet_pt", ak4Jets[1].pt, DL_resolved_2b_emu, EqBin(
                    500, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)"),
                Plot.make1D("DL_resolved_2b_emu_InvM_leadingJet_eta", ak4Jets[0].eta, DL_resolved_2b_emu, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("DL_resolved_2b_emu_InvM_subleadingJet_eta", ak4Jets[1].eta, DL_resolved_2b_emu, EqBin(
                    80, -3, 3), title="eta(j2)", xTitle="eta(j2)"),
                Plot.make1D("DL_resolved_2b_emu_DR_jets", op.deltaR(ak4Jets[0].p4, ak4Jets[1].p4), DL_resolved_2b_emu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)"),
                Plot.make1D("DL_resolved_2b_emu_InvM_emu", op.invariant_mass(ElMuDileptonPreSel[0][0].p4, ElMuDileptonPreSel[0][1].p4), DL_resolved_2b_emu, EqBin(
                    160, 40., 200.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pairs (resolved) (GeV/c^2)")

            ])
        if self.channel == "SL":
            plots.extend([
                # SL boosted plots
                Plot.make1D("SL_boosted_fatJet_pt", ak8BJets[0].pt, SL_boosted, EqBin(
                    400, 200, 1000), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet1_pt", ak8BJets[0].subJet1.pt, SL_boosted, EqBin(
                    250, 0, 500), title=" pT(j1 subjet1)", xTitle="pT(j1 subjet1) (GeV/c)"),
                Plot.make1D("SL_boosted_subjet2_pt", ak8BJets[0].subJet2.pt, SL_boosted, EqBin(
                    250, 0, 500), title=" pT(j1 subjet2)", xTitle="pT(j1 subjet2) (GeV/c)"),
                Plot.make1D("SL_boosted_fatJet_eta", ak8BJets[0].eta, SL_boosted, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)"),
                Plot.make1D("SL_boosted_subjet1_eta", ak8BJets[0].subJet1.eta, SL_boosted, EqBin(
                    80, -3, 3), title="eta(j1 subjet1)", xTitle="eta(j1 subjet1)"),
                Plot.make1D("SL_boosted_subjet2_eta", ak8BJets[0].subJet2.eta, SL_boosted, EqBin(
                    80, -3, 3), title="eta(j1 subjet2)", xTitle="eta(j1 subjet2)"),
                Plot.make1D("SL_boosted_InvM_jj", op.invariant_mass(ak8BJets[0].subJet1.p4, ak8BJets[0].subJet2.p4), SL_boosted, EqBin(
                    160, 40., 200.), title="InvM(jj)", xTitle="Invariant Mass of jets (GeV/c^2)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet1_eta", [op.invariant_mass(ak8BJets[0].subJet1.p4, ak8BJets[0].subJet2.p4), ak8BJets[0].eta], SL_boosted, [
                    EqBin(160, 40., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet1 eta", xTitle="Invariant Mass of jets (GeV/c^2)", yTitle="eta(j1)"),
                Plot.make2D("SL_boosted_InvM_jj_vs_jet2_eta", [op.invariant_mass(ak8BJets[0].subJet1.p4, ak8BJets[0].subJet2.p4), ak8BJets[0].subJet2.eta], SL_boosted, [
                    EqBin(160, 40., 200.), EqBin(-8, -3, 3)], title="InvM(jj) vs jet2 eta", xTitle="Invariant Mass of jets (GeV/c^2)", yTitle="eta(j2)"),

                # SL resolved plots
                Plot.make1D("SL_resolved_nJets", op.rng_len(ak4BJets), SL_resolved, EqBin(
                    15, 0., 15.), xTitle="Number of jets"),
                Plot.make1D("SL_resolved_InvM_leadingJet_pt", ak4BJets[0].pt, SL_resolved, EqBin(
                    500, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)"),
                Plot.make1D("SL_resolved_InvM_leadingJet_eta", ak4BJets[0].eta, SL_resolved, EqBin(
                    80, -3, 3), title="eta(j1)", xTitle="eta(j1)")
            ])

        return plots
