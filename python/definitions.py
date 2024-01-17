from bamboo import treefunctions as op

# helper functions

def labeler(label):
    return {'labels': [{'text': label, 'position': [0.235, 0.9], 'size': 24}]}

# Lepton functions

def hasAssociatedJet(lep): return lep.jet.idx != -1


def muon_x(mu): return op.min(
    op.max(0., (0.9*mu.pt*(1+mu.jetRelIso))-20.)/(45.-20.), 1.)


def muon_btagInterpolation(mu): return muon_x(
    mu)*0.0494 + (1-muon_x(mu))*0.2770  # 2018 values


def muon_deepJetInterpIfMvaFailed(mu): return op.OR(op.NOT(
    hasAssociatedJet(mu)), mu.jet.btagDeepFlavB < muon_btagInterpolation(mu))


def lepton_associatedJetLessThanMediumBtag(lep): return op.OR(
    op.NOT(hasAssociatedJet(lep)), lep.jet.btagDeepFlavB <= 0.2770)  # 2018 value


def lepton_associatedJetLessThanTightBtag(lep): return op.OR(
    op.NOT(hasAssociatedJet(lep)), lep.jet.btagDeepFlavB <= 0.7264)  # 2018 value


def electronTightSel(el): return el.mvaTTH >= 0.30

# Object definitions


def muonDef(muons):
    return op.select(muons, lambda mu :op.AND(
        mu.pt >= 5.,
        op.abs(mu.eta) <= 2.4,
        op.abs(mu.dxy) <= 0.05,
        op.abs(mu.dz) <= 0.1,
        mu.miniPFRelIso_all <= 0.4,
        mu.sip3d <= 8,
        mu.looseId
    ))


def muonConePt(muons):
    return op.map(muons, lambda lep: op.multiSwitch(
        (op.AND(op.abs(lep.pdgId) != 11, op.abs(lep.pdgId) != 13), lep.pt),
        (op.AND(op.abs(lep.pdgId) == 13, lep.mediumId, lep.mvaTTH > 0.50), lep.pt),
        0.9*lep.pt*(1.+lep.jetRelIso)
    ))


def muonFakeSel(muons):
    return op.select(muons, lambda mu: op.AND(
        muonConePt(muons)[mu.idx] >= 10.,
        op.OR(lepton_associatedJetLessThanMediumBtag(mu), op.AND(mu.jetRelIso < 0.8, muon_deepJetInterpIfMvaFailed(mu))))
    )


def muonTightSel(muons): return op.select(muons, lambda mu: op.AND(
            muonConePt(muons)[mu.idx] >= 10.,
            lepton_associatedJetLessThanMediumBtag(mu),
            mu.mvaTTH >= 0.50,
            mu.mediumId
            ))


def elDef(electrons):
    return op.select(electrons, lambda el: op.AND(
        el.pt >= 7.,
        op.abs(el.eta) <= 2.5,
        op.abs(el.dxy) <= 0.05,
        op.abs(el.dz) <= 0.1,
        el.sip3d <= 8,
        el.miniPFRelIso_all <= 0.4,
        # el.mvaNoIso > < VALUE TO BE DETERMINED >,
        el.lostHits <= 1
    ))


def elConePt(electrons):
    return op.map(electrons, lambda lep: op.multiSwitch(
        (op.AND(op.abs(lep.pdgId) != 11, op.abs(lep.pdgId) != 13), lep.pt),
        (op.AND(op.abs(lep.pdgId) == 11, lep.mvaTTH > 0.30), lep.pt),
        0.9*lep.pt*(1.+lep.jetRelIso)
    ))


def cleanElectrons(electrons, muons):
    cleanedElectrons = op.select(electrons, lambda el: op.NOT(
        op.rng_any(
            muons, lambda mu: op.deltaR(el.p4, mu.p4) <= 0.3))
    )
    return cleanedElectrons


def elFakeSel(electrons):
    return op.select(electrons, lambda el: op.AND(
        elConePt(electrons)[el.idx] >= 10,
        op.OR(
            op.AND(op.abs(el.eta+el.deltaEtaSC) <= 1.479, el.sieie <= 0.011),
            op.AND(op.abs(el.eta+el.deltaEtaSC) > 1.479, el.sieie <= 0.030)
        ),
        el.hoe <= 0.10,
        el.eInvMinusPInv >= -0.04,
        op.OR(el.mvaTTH >= 0.30, op.AND(el.jetRelIso < 0.7, el.mvaNoIso_WP90)),
        op.switch(
            el.mvaTTH < 0.30,
            lepton_associatedJetLessThanTightBtag(el),
            lepton_associatedJetLessThanMediumBtag(el)),
        el.lostHits == 0,
        el.convVeto
    ))


def elTightSel(electrons): return op.select(electrons, lambda el: op.AND(
            elConePt(electrons)[el.idx] >= 10.,
            op.OR(
                op.AND(op.abs(el.eta+el.deltaEtaSC) <= 1.479, el.sieie <= 0.011),
                op.AND(op.abs(el.eta+el.deltaEtaSC) > 1.479, el.sieie <= 0.030)
            ),
            el.hoe <= 0.10,
            el.eInvMinusPInv >= -0.04,
            el.convVeto,
            el.mvaTTH >= 0.30,
            el.lostHits == 0,
            lepton_associatedJetLessThanMediumBtag(el),
            ))


def ak4jetDef(jets):
    return op.select(jets, lambda jet: op.AND(
        jet.jetId & 2,  # tight
        jet.pt >= 25.,
        op.abs(jet.eta) <= 2.4,
        # op.OR(((jet.puId >> 2) & 1), jet.pt > 50.) # Jet PU ID bit1 is loose # no puId in Run3 so far
    ))


def ak8jetDef(jets):
    return op.select(jets, lambda jet: op.AND(
        jet.pt >= 200.,
        op.abs(jet.eta) <= 2.4,
        jet.jetId & 2,  # tight
        jet.subJet1.isValid,
        jet.subJet2.isValid,
        jet.subJet1.pt >= 20.,
        jet.subJet2.pt >= 20.,
        op.abs(jet.subJet1.eta) <= 2.4,
        op.abs(jet.subJet2.eta) <= 2.4,
        op.AND(jet.msoftdrop >= 30., jet.msoftdrop <= 210.),
        jet.tau2 / jet.tau1 <= 0.75
    ))
    
# bTagging for ak4 jets
def ak4BtagSel(jet): return jet.btagDeepFlavB > 0.2770


def tauDef(taus):
    return op.select(taus, lambda tau: op.AND(
        tau.pt > 20.,
        op.abs(tau.eta) < 2.3,
        op.abs(tau.dxy) <= 1000.0,
        op.abs(tau.dz) <= 0.2,
        tau.idDecayModeOldDMs,
        op.OR(tau.decayMode == 0,
              tau.decayMode == 1,
              tau.decayMode == 2,
              tau.decayMode == 10,
              tau.decayMode == 11),
        (tau.idDeepTau2017v2p1VSjet >> 4 & 0x1) == 1,
        (tau.idDeepTau2017v2p1VSe >> 0 & 0x1) == 1,
        (tau.idDeepTau2017v2p1VSmu >> 0 & 0x1) == 1
    ))

def cleanTaus(taus, electrons, muons):
    return op.select(taus, lambda tau: op.AND(
            op.NOT(op.rng_any(
                electrons, lambda el: op.deltaR(tau.p4, el.p4) <= 0.3)),
            op.NOT(op.rng_any(
                muons, lambda mu: op.deltaR(tau.p4, mu.p4) <= 0.3))
        ))

# remove jets within cone of DR<0.4 of leading leptons at each channel
def cleaningWithRespectToLeadingLepton(electrons, muons, DR):
    return lambda jet: op.multiSwitch(
        (op.AND(op.rng_len(electrons) >= 1, op.rng_len(
            muons) == 0), op.deltaR(jet.p4, electrons[0].p4) >= DR),
        (op.AND(op.rng_len(electrons) == 0, op.rng_len(
            muons) >= 1), op.deltaR(jet.p4, muons[0].p4) >= DR),
        (op.AND(op.rng_len(muons) >= 1, op.rng_len(electrons) >= 1), op.switch(
            elConePt(electrons)[0] >= muonConePt(muons)[0],
            op.deltaR(jet.p4, electrons[0].p4) >= DR,
            op.deltaR(jet.p4, muons[0].p4) >= DR)),
        op.c_bool(True)
    )

def cleaningWithRespectToLeadingLeptons(electrons, muons, DR):
    return lambda j: op.multiSwitch(
        # Only electrons
        (op.AND(op.rng_len(electrons) >= 2, op.rng_len(muons) == 0),
            op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, electrons[1].p4) >= DR)),
        # Only muons
        (op.AND(op.rng_len(electrons) == 0, op.rng_len(muons) >= 2),
            op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(j.p4, muons[1].p4) >= DR)),
        # One electron + one muon
        (op.AND(op.rng_len(electrons) == 1, op.rng_len(muons) == 1),
            op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, muons[0].p4) >= DR)),
        # At least one electron + at least one muon
        (op.AND(op.rng_len(electrons) >= 1, op.rng_len(muons) >= 1),
            op.switch(
            # Electron is the leading lepton
            elConePt(electrons)[0] > muonConePt(muons)[0],
            op.switch(op.rng_len(electrons) == 1,
                        op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(
                            j.p4, muons[0].p4) >= DR),
                        op.switch(elConePt(electrons)[1] > muonConePt(muons)[0],
                                op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(
                                    j.p4, electrons[1].p4) >= DR),
                                op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, muons[0].p4) >= DR))),
            # Muon is the leading lepton
            op.switch(op.rng_len(muons) == 1,
                        op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(
                            j.p4, electrons[0].p4) >= DR),
                        op.switch(muonConePt(muons)[1] > elConePt(electrons)[0],
                                op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(
                                    j.p4, muons[1].p4) >= DR),
                                op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(j.p4, electrons[0].p4) >= DR))))),
        op.c_bool(True)
    )

def defineObjects(self, tree):
    # lepton cone-pt definitions
    self.muon_conept = muonConePt(tree.Muon)
    self.electron_conept = elConePt(tree.Electron)
    
    # lepton definitions sorted by their cone-pt
    self.muons = op.sort(muonDef(tree.Muon), lambda mu: -self.muon_conept[mu.idx])
    self.electrons = op.sort(elDef(tree.Electron), lambda el: -self.electron_conept[el.idx]) # can be liberated from 'self' ?

    # cleaning electrons wrt muons
    self.clElectrons = cleanElectrons(self.electrons, self.muons)

    # Fakeable leptons
    self.fakeMuons = muonFakeSel(self.muons)
    self.fakeElectrons = elFakeSel(self.clElectrons)

    # tight leptons
    self.tightMuons = muonTightSel(self.muons)
    self.tightElectrons = elTightSel(self.clElectrons)

    # Taus
    taus = tauDef(tree.Tau)
    self.cleanedTaus = cleanTaus(taus, self.fakeElectrons, self.fakeMuons)

    # AK4 Jets sorted by their pt
    ak4JetsPreSel = op.sort(ak4jetDef(tree.Jet), lambda jet: -jet.pt)
    
    # clean jets wrt leptons
    if self.channel == 'DL':
        self.cleanAk4Jets = cleaningWithRespectToLeadingLeptons(self.fakeElectrons, self.fakeMuons, 0.4)
        self.cleanAk8Jets = cleaningWithRespectToLeadingLeptons(self.fakeElectrons, self.fakeMuons, 0.8)
        
    if self.channel == 'SL':
        self.cleanAk4Jets = cleaningWithRespectToLeadingLepton(self.fakeElectrons, self.fakeMuons, 0.4)
        self.cleanAk8Jets = cleaningWithRespectToLeadingLepton(self.fakeElectrons, self.fakeMuons, 0.8)
    
    self.ak4Jets = op.select(ak4JetsPreSel, self.cleanAk4Jets)
    self.ak4JetsByBtagScore = op.sort(self.ak4Jets, lambda j: -j.btagDeepFlavB)
    
    self.ak4BJets = op.select(self.ak4Jets, ak4BtagSel)
    
    # AK8 Jets
    self.ak8JetsDef = ak8jetDef(tree.FatJet)

    if self.channel == 'SL': # sorted by btag score
        ak8JetsPreSel = op.sort(self.ak8JetsDef, lambda j: -j.btagDeepB)
    if self.channel == 'DL': # sorted by pt
        ak8JetsPreSel = op.sort(self.ak8JetsDef, lambda j: -j.pt)
        
    self.ak8Jets = op.select(ak8JetsPreSel, self.cleanAk8Jets)
    
    # 2018 DeepJet WP
    def subjetBtag(subjet): return subjet.btagDeepB > 0.4184

    def ak8Btag(fatjet): return op.OR(op.AND(fatjet.subJet1.pt >= 30, subjetBtag(fatjet.subJet1)),
                                        op.AND(fatjet.subJet2.pt >= 30, subjetBtag(fatjet.subJet2)))

    self.ak8BJets = op.select(self.ak8Jets, ak8Btag)

    # Ak4 Jet Collection cleaned from Ak8b #
    def cleanAk4FromAk8b(ak4j): return op.AND(op.rng_len(
        self.ak8BJets) > 0, op.deltaR(ak4j.p4, self.ak8BJets[0].p4) > 1.2)
    
    self.ak4JetsCleanedFromAk8b = op.select(self.ak4Jets, cleanAk4FromAk8b)