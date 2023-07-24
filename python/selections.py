from bamboo import treefunctions as op

# common variables for DL and SL channels
Zmass = 91.1876

def lowMllCut(dileptons): return op.NOT(op.rng_any(
    dileptons, lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4) < 12.))

def outZ(dileptons): return op.NOT(op.rng_any(
    dileptons, lambda dilep: op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4) - Zmass) < 10.))

# end common variables

def makeDLSelection(self, noSel):
    """Selections for the DL channel
    return the following selections:
    - DL_boosted_ee: boosted selection for ee channel
    - DL_boosted_mumu: boosted selection for mumu channel
    - DL_boosted_emu: boosted selection for emu channel
    - DL_resolved_1b_ee: resolved selection for ee channel with at least one b-tagged ak4 jet
    - DL_resolved_1b_mumu: resolved selection for mumu channel with at least one b-tagged ak4 jet
    - DL_resolved_1b_emu: resolved selection for emu channel with at least one b-tagged ak4 jet
    - DL_resolved_2b_ee: resolved selection for ee channel with at least two b-tagged ak4 jets
    - DL_resolved_2b_mumu: resolved selection for mumu channel with at least two b-tagged ak4 jets
    - DL_resolved_2b_emu: resolved selection for emu channel with at least two b-tagged ak4 jets
    """

    def OSDilepton(dilep): return dilep[0].charge != dilep[1].charge

    # Pt cuts #
    def lowPtCutElEl(dilep): return op.AND(
        self.electron_conept[dilep[0].idx] > 15, self.electron_conept[dilep[1].idx] > 15)  # subleading above 15 GeV

    def lowPtCutMuMu(dilep): return op.AND(
        self.muon_conept[dilep[0].idx] > 15, self.muon_conept[dilep[1].idx] > 15)  # subleading above 15 GeV

    def lowPtCutElMu(dilep): return op.AND(
        self.electron_conept[dilep[0].idx] > 15, self.muon_conept[dilep[1].idx] > 15)  # subleading above 15 GeV

    def leadingPtCutElEl(dilep): return op.OR(
        self.electron_conept[dilep[0].idx] > 25, self.electron_conept[dilep[1].idx] > 25)  # leading above 25 GeV

    def leadingPtCutMuMu(dilep): return op.OR(
        self.muon_conept[dilep[0].idx] > 25, self.muon_conept[dilep[1].idx] > 25)  # leading above 25 GeV

    def leadingPtCutElMu(dilep): return op.OR(
        self.electron_conept[dilep[0].idx] > 25, self.muon_conept[dilep[1].idx] > 25)  # leading above 25 GeV

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

    def leptonOS(l1, l2): return l1.charge != l2.charge

    self.ElElDileptonPreSel = op.combine(self.electrons, N=2)
    self.MuMuDileptonPreSel = op.combine(self.muons, N=2)
    self.ElMuDileptonPreSel = op.combine((self.electrons, self.muons))

    OSElElDileptonPreSel = op.combine(self.electrons, N=2, pred=leptonOS)
    OSMuMuDileptonPreSel = op.combine(self.muons, N=2, pred=leptonOS)
    OSElMuDileptonPreSel = op.combine((self.electrons, self.muons), pred=leptonOS)

    # Dilepton for selection #
    ElElFakeSel = op.combine(self.fakeElectrons, N=2)
    MuMuFakeSel = op.combine(self.fakeMuons, N=2)
    ElMuFakeSel = op.combine((self.fakeElectrons, self.fakeMuons))

    ElElTightSel = op.combine(self.tightElectrons, N=2)
    MuMuTightSel = op.combine(self.tightMuons, N=2)
    ElMuTightSel = op.combine((self.tightElectrons, self.tightMuons))

    # Fake lepton selection #
    # elelSel = noSel.refine('elelSel', cut=[op.rng_len(self.ElElFakeSel) >= 1,
    #                                        op.OR(op.rng_len(self.fakeMuons) == 0,
    #                                              op.AND(op.rng_len(self.fakeMuons) == 1,
    #                                                     self.electron_conept[self.ElElFakeSel[0]
    #                                                                     [0].idx] > self.muon_conept[self.fakeMuons[0].idx],
    #                                                     self.electron_conept[self.ElElFakeSel[0][1].idx] > self.muon_conept[self.fakeMuons[0].idx]),
    #                                              op.AND(op.rng_len(self.fakeMuons) >= 2,
    #                                                     self.electron_conept[self.ElElFakeSel[0]
    #                                                                     [0].idx] > self.muon_conept[self.fakeMuons[0].idx],
    #                                                     self.electron_conept[self.ElElFakeSel[0]
    #                                                                     [1].idx] > self.muon_conept[self.fakeMuons[0].idx],
    #                                                     self.electron_conept[self.ElElFakeSel[0]
    #                                                                     [0].idx] > self.muon_conept[self.fakeMuons[1].idx],
    #                                                     self.electron_conept[self.ElElFakeSel[0][1].idx] > self.muon_conept[self.fakeMuons[1].idx]))])

    # mumuSel = noSel.refine('mumuSel', cut=[op.rng_len(self.MuMuFakeSel) >= 1,
    #                                        op.OR(op.rng_len(self.fakeElectrons) == 0,
    #                                              op.AND(op.rng_len(self.fakeElectrons) == 1,
    #                                                     self.muon_conept[self.MuMuFakeSel[0]
    #                                                                 [0].idx] > self.electron_conept[self.fakeElectrons[0].idx],
    #                                                     self.muon_conept[self.MuMuFakeSel[0][1].idx] > self.electron_conept[self.fakeElectrons[0].idx]),
    #                                              op.AND(op.rng_len(self.fakeElectrons) >= 2,
    #                                                     self.muon_conept[self.MuMuFakeSel[0]
    #                                                                 [0].idx] > self.electron_conept[self.fakeElectrons[0].idx],
    #                                                     self.muon_conept[self.MuMuFakeSel[0]
    #                                                                 [1].idx] > self.electron_conept[self.fakeElectrons[0].idx],
    #                                                     self.muon_conept[self.MuMuFakeSel[0]
    #                                                                 [0].idx] > self.electron_conept[self.fakeElectrons[1].idx],
    #                                                     self.muon_conept[self.MuMuFakeSel[0][1].idx] > self.electron_conept[self.fakeElectrons[1].idx]))])
    # elmuSel = noSel.refine('elmuSel', cut=[op.rng_len(self.ElMuFakeSel) >= 1,
    #                                        op.OR(op.AND(op.rng_len(self.fakeElectrons) == 1,
    #                                                     op.rng_len(self.fakeMuons) == 1),
    #                                              op.AND(op.rng_len(self.fakeElectrons) >= 2,
    #                                                     op.rng_len(
    #                                                  self.fakeMuons) == 1,
    #                                            self.muon_conept[self.ElMuFakeSel[0][1].idx] > self.electron_conept[self.fakeElectrons[1].idx]),
    #     op.AND(op.rng_len(self.fakeMuons) >= 2,
    #                                            op.rng_len(
    #         self.fakeElectrons) == 1,
    #                                            self.electron_conept[self.ElMuFakeSel[0][0].idx] > self.muon_conept[self.fakeMuons[1].idx]),
    #     op.AND(op.rng_len(self.fakeElectrons) >= 2,
    #                                            op.rng_len(
    #         self.fakeMuons) >= 2,
    #                                            self.muon_conept[self.ElMuFakeSel[0]
    #                                                        [1].idx] > self.electron_conept[self.fakeElectrons[1].idx],
    #                                            self.electron_conept[self.ElMuFakeSel[0][0].idx] > self.muon_conept[self.fakeMuons[1].idx]))])

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
    mllCut = [lowMllCut(self.ElElDileptonPreSel), lowMllCut(
        self.MuMuDileptonPreSel), lowMllCut(self.ElMuDileptonPreSel)]
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
        op.rng_len(self.tightElectrons) == 2,
        op.rng_len(self.tightMuons) == 0,
        ElElTightSel[0][0].idx == ElElFakeSel[0][0].idx,
        ElElTightSel[0][1].idx == ElElFakeSel[0][1].idx]
    )
    mumuSel.refine('dileptonCut_mumu', cut=[
        tightpair_MuMu(MuMuFakeSel[0]),
        op.rng_len(self.tightMuons) == 2,
        op.rng_len(self.tightElectrons) == 0,
        MuMuTightSel[0][0].idx == MuMuFakeSel[0][0].idx,
        MuMuTightSel[0][1].idx == MuMuFakeSel[0][1].idx]
    )
    elmuSel.refine('dileptonCut_emu', cut=[
        tightpair_ElMu(ElMuFakeSel[0]),
        op.rng_len(self.tightElectrons) == 1,
        op.rng_len(self.tightMuons) == 1,
        ElMuTightSel[0][0].idx == ElMuFakeSel[0][0].idx,
        ElMuTightSel[0][1].idx == ElMuFakeSel[0][1].idx]
    )

    # boosted -> and at least one b-tagged ak8 jet
    DL_boosted_ee = elelSel.refine(
        'DL_boosted_ee', cut=(op.rng_len(self.ak8BJets) >= 1))
    DL_boosted_mumu = mumuSel.refine(
        'DL_boosted_mumu', cut=(op.rng_len(self.ak8BJets) >= 1))
    DL_boosted_emu = elmuSel.refine(
        'DL_boosted_emu', cut=(op.rng_len(self.ak8BJets) >= 1))

    # resolved -> and at least two ak4 jets with at least one b-tagged and no ak8 jets
    DL_resolved_1b_ee = elelSel.refine('DL_resolved_1b_ee', cut=(op.AND(op.rng_len(self.ak4BJets) >= 1, op.rng_len(self.ak8Jets) == 0)))
    DL_resolved_1b_mumu = mumuSel.refine('DL_resolved_1b_mumu', cut=(op.AND(op.rng_len(self.ak4BJets) >= 1, op.rng_len(self.ak8Jets) == 0)))
    DL_resolved_1b_emu = elmuSel.refine('DL_resolved_1b_emu', cut=(op.AND(op.rng_len(self.ak4BJets) >= 1, op.rng_len(self.ak8Jets) == 0)))

    DL_resolved_2b_ee = DL_resolved_1b_ee.refine('DL_resolved_ee', cut=(op.rng_len(self.ak4BJets) >= 2))
    DL_resolved_2b_mumu = DL_resolved_1b_mumu.refine('DL_resolved_mumu', cut=(op.rng_len(self.ak4BJets) >= 2))
    DL_resolved_2b_emu = DL_resolved_1b_emu.refine('DL_resolved_emu', cut=(op.rng_len(self.ak4BJets) >= 2))

    DL_selections = [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu, DL_resolved_1b_ee, DL_resolved_1b_mumu, DL_resolved_1b_emu, DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu]
        
    return DL_selections

def makeSLSelection(self, noSel):
    """ Selections for the SL channel
    return the following selections:
    - SL_resolved: resolved selection
    - SL_resolved_e: resolved selection for e channel
    - SL_resolved_mu: resolved selection for mu channel
    - SL_boosted: boosted selection"""
    def elPtCut(lep): return self.electron_conept[lep[0].idx] > 32.0

    def muPtCut(lep): return self.muon_conept[lep[0].idx] > 25.0

    def tau_h_veto(taus): return op.rng_len(taus) == 0
    
    def leptonOS(l1, l2): return l1.charge != l2.charge

    ElElDileptonPreSel = op.combine(self.clElectrons, N=2)
    MuMuDileptonPreSel = op.combine(self.muons, N=2)
    ElMuDileptonPreSel = op.combine((self.clElectrons, self.muons))

    OSElElDileptonPreSel = op.combine(self.clElectrons, N=2, pred=leptonOS)
    OSMuMuDileptonPreSel = op.combine(self.muons, N=2, pred=leptonOS)
    
    outZcut = [outZ(OSElElDileptonPreSel), outZ(OSMuMuDileptonPreSel)]

    mllCut = [lowMllCut(ElElDileptonPreSel), lowMllCut(
        MuMuDileptonPreSel), lowMllCut(ElMuDileptonPreSel)]

    SL_resolved = noSel.refine('SL_resolved', cut=[
        mllCut, outZcut, tau_h_veto(self.cleanedTaus),
        op.rng_len(self.ak4Jets) >= 3,
        op.rng_len(self.ak4BJets) >= 1,
        op.rng_len(self.ak8BJets) == 0])

    SL_resolved_e = SL_resolved.refine('SL_resolved_e', cut=[
        elPtCut(self.tightElectrons),
        op.rng_len(self.tightElectrons) == 1,
        op.rng_len(self.tightMuons) == 0])

    SL_resolved_mu = SL_resolved.refine('SL_resolved_mu', cut=[
        muPtCut(self.tightMuons),
        op.rng_len(self.tightElectrons) == 0,
        op.rng_len(self.tightMuons) == 1])

    SL_boosted = noSel.refine('SL_boosted', cut=[
        op.OR(elPtCut, muPtCut), lowMllCut, outZ, tau_h_veto,
        op.rng_len(self.ak8BJets) >= 1,
        op.rng_len(self.ak4JetsCleanedFromAk8b) >= 1])
    
    SL_boosted_e = SL_boosted.refine('SL_boosted_e', cut=[
        op.rng_len(self.tightElectrons) == 1,
        op.rng_len(self.tightMuons) == 0])
    
    SL_boosted_mu = SL_boosted.refine('SL_boosted_mu', cut=[
        op.rng_len(self.tightMuons) == 1,
        op.rng_len(self.tightElectrons) == 0])

    SL_selections = [SL_resolved, SL_resolved_e, SL_resolved_mu, SL_boosted, SL_boosted_e, SL_boosted_mu]

    return SL_selections