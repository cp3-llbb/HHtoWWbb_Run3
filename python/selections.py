from bamboo import treefunctions as op

# common variables for DL and SL channels

Zmass = 91.1876


def lowMllCut(dileptons): return op.NOT(op.rng_any(
    dileptons, lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4) < 12.))


def outZ(dileptons): return op.NOT(op.rng_any(
    dileptons, lambda dilep: op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4) - Zmass) <= 10.))

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
    """

    # Pt cuts : subleading above 15 GeV and leading above 25 GeV
    def ptCutElEl(dilep): return op.AND(
        self.electron_conept[dilep[0].idx] > 15,
        self.electron_conept[dilep[1].idx] > 15,
        op.OR(self.electron_conept[dilep[0].idx] >
              25, self.electron_conept[dilep[1].idx] > 25)
    )

    def ptCutMuMu(dilep): return op.AND(
        self.muon_conept[dilep[0].idx] > 15,
        self.muon_conept[dilep[1].idx] > 15,
        op.OR(self.muon_conept[dilep[0].idx] > 25,
              self.muon_conept[dilep[1].idx] > 25)
    )

    def ptCutElMu(dilep): return op.AND(
        self.electron_conept[dilep[0].idx] > 15,
        self.muon_conept[dilep[1].idx] > 15,
        op.OR(self.electron_conept[dilep[0].idx] >
              25, self.muon_conept[dilep[1].idx] > 25)
    )

    # OS loose lepton pairs of same type to be vetoed around Z peak
    ElElLooseSel = op.combine(
        self.clElectrons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    MuMuLooseSel = op.combine(
        self.muons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    ElMuLooseSel = op.combine(
        (self.clElectrons, self.muons), N=2, pred=lambda el, mu: el.charge != mu.charge)

    # OS tight dilepton collections
    ElElTightSel = op.combine(self.tightElectrons, N=2,
                              pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    MuMuTightSel = op.combine(
        self.tightMuons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    ElMuTightSel = op.combine((self.tightElectrons, self.tightMuons),
                              N=2, pred=lambda el, mu: el.charge != mu.charge)

    self.firstOSElEl = ElElTightSel[0]
    self.firstOSMuMu = MuMuTightSel[0]
    self.firstOSElMu = ElMuTightSel[0]

    # minimum pT cut : at least one lepton pair with leading lepton above 25 GeV
    elelPtSel = noSel.refine('elelptSel', cut=[ptCutElEl(self.firstOSElEl)])
    mumuPtSel = noSel.refine('mumuptSel', cut=[ptCutMuMu(self.firstOSMuMu)])
    elmuPtSel = noSel.refine('elmuptSel', cut=[ptCutElMu(self.firstOSElMu)])

    # low Mll cut : reject events with dilepton mass below 12 GeV
    mllCut = op.AND(lowMllCut(ElElLooseSel), lowMllCut(
        MuMuLooseSel), lowMllCut(ElMuLooseSel))

    # Z-veto : reject events with dileptons of same type with mass around Z peak
    outZCut = op.AND(outZ(ElElLooseSel), outZ(MuMuLooseSel))

    OSoutZelelSel = elelPtSel.refine(
        'OSoutZelelSel', cut=op.AND(mllCut, outZCut))
    OSoutZmumuSel = mumuPtSel.refine(
        'OSoutZmumuSel', cut=op.AND(mllCut, outZCut))
    OSoutZelmuSel = elmuPtSel.refine(
        'OSoutZelmuSel', cut=op.AND(mllCut, outZCut))

    # di-lepton multiplicity cut
    leptonMultiplicityCut_ee = OSoutZelelSel.refine('dileptonCut_ee', cut=[op.AND(
        op.rng_len(ElElTightSel) == 1,
        op.rng_len(MuMuTightSel) == 0,
        op.rng_len(ElMuTightSel) == 0
    )])
    leptonMultiplicityCut_mumu = OSoutZmumuSel.refine('dileptonCut_mumu', cut=[op.AND(
        op.rng_len(ElElTightSel) == 0,
        op.rng_len(MuMuTightSel) == 1,
        op.rng_len(ElMuTightSel) == 0
    )])
    leptonMultiplicityCut_emu = OSoutZelmuSel.refine('dileptonCut_emu', cut=[op.AND(
        op.rng_len(ElElTightSel) == 0,
        op.rng_len(MuMuTightSel) == 0,
        op.rng_len(ElMuTightSel) == 1,
    )])

    # boosted -> at least one b-tagged ak8 jet
    DL_boosted_ee = leptonMultiplicityCut_ee.refine(
        'DL_boosted_ee', cut=(op.rng_len(self.ak8BJets) >= 1))
    DL_boosted_mumu = leptonMultiplicityCut_mumu.refine(
        'DL_boosted_mumu', cut=(op.rng_len(self.ak8BJets) >= 1))
    DL_boosted_emu = leptonMultiplicityCut_emu.refine(
        'DL_boosted_emu', cut=(op.rng_len(self.ak8BJets) >= 1))

    # resolved -> and at least two ak4 jets with at least one b-tagged and no ak8 jets
    DL_resolved_ee = leptonMultiplicityCut_ee.refine(
        'DL_resolved_ee',
        cut=(op.AND(
            op.rng_len(self.ak4Jets) >= 2,
            op.rng_len(self.ak4BJets) >= 1,
            op.rng_len(self.ak8Jets) == 0))
    )

    DL_resolved_mumu = leptonMultiplicityCut_mumu.refine(
        'DL_resolved_mumu',
        cut=(op.AND(
            op.rng_len(self.ak4Jets) >= 2,
            op.rng_len(self.ak4BJets) >= 1,
            op.rng_len(self.ak8Jets) == 0))
    )

    DL_resolved_emu = leptonMultiplicityCut_emu.refine(
        'DL_resolved_emu',
        cut=(op.AND(
            op.rng_len(self.ak4Jets) >= 2,
            op.rng_len(self.ak4BJets) >= 1,
            op.rng_len(self.ak8Jets) == 0))
    )

    DL_selections = [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu,
                     DL_resolved_ee, DL_resolved_mumu, DL_resolved_emu]

    return DL_selections


def makeSLSelection(self, noSel):
    """ Selections for the SL channel
    return the following selections:
    - SL_resolved: resolved selection
    - SL_resolved_e: resolved selection for e channel
    - SL_resolved_mu: resolved selection for mu channel
    - SL_boosted: boosted selection
    - SL_boosted_e: boosted selection for e channel
    - SL_boosted_mu: boosted selection for mu channel"""

    def elPtCut(lep): return self.electron_conept[lep[0].idx] > 32.0

    def muPtCut(lep): return self.muon_conept[lep[0].idx] > 25.0

    def tau_h_veto(taus): return op.rng_len(taus) == 0

    # OS loose lepton pairs of same type to be vetoed around Z peak
    ElElLooseSel = op.combine(
        self.clElectrons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    MuMuLooseSel = op.combine(
        self.muons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    ElMuLooseSel = op.combine(
        (self.clElectrons, self.muons), N=2, pred=lambda el, mu: el.charge != mu.charge)

    # Z-veto : reject events with dileptons of same type with mass around Z peak
    outZCut = op.AND(outZ(ElElLooseSel), outZ(MuMuLooseSel))

    # low Mll cut : reject events with dilepton mass below 12 GeV
    mllCut = op.AND(lowMllCut(ElElLooseSel), lowMllCut(
        MuMuLooseSel), lowMllCut(ElMuLooseSel))

    OSoutZelelSel = noSel.refine('OSoutZsel', cut=op.AND(
        mllCut, outZCut, tau_h_veto(self.cleanedTaus)))

    SL_resolved = OSoutZelelSel.refine('SL_resolved', cut=[
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

    SL_boosted = OSoutZelelSel.refine('SL_boosted', cut=[
        op.rng_len(self.ak8BJets) >= 1,
        op.rng_len(self.ak4JetsCleanedFromAk8b) >= 1])

    SL_boosted_e = SL_boosted.refine('SL_boosted_e', cut=[
        elPtCut(self.tightElectrons),
        op.rng_len(self.tightElectrons) == 1,
        op.rng_len(self.tightMuons) == 0])

    SL_boosted_mu = SL_boosted.refine('SL_boosted_mu', cut=[
        muPtCut(self.tightMuons),
        op.rng_len(self.tightMuons) == 1,
        op.rng_len(self.tightElectrons) == 0])

    SL_selections = [SL_resolved, SL_resolved_e,
                     SL_resolved_mu, SL_boosted, SL_boosted_e, SL_boosted_mu]

    return SL_selections
