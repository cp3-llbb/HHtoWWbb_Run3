from bamboo import treefunctions as op

# Object definitions


def muonDef(mu):
    return op.AND(
        mu.pt > 5.,
        op.abs(mu.eta) < 2.4,
        op.abs(mu.dxy) < 0.05,
        op.abs(mu.dz) < 0.1,
        mu.miniPFRelIso_all < 0.4,
        mu.sip3d < 8,
        mu.looseId
    )


def muonConePt(muons):
    return op.map(muons, lambda lep: op.multiSwitch(
        (op.AND(op.abs(lep.pdgId) != 11, op.abs(lep.pdgId) != 13), lep.pt),
        (op.AND(op.abs(lep.pdgId) == 13, lep.mvaTTH > 0.50), lep.pt),
        0.9*lep.pt*(1.+lep.jetRelIso)
    ))

def muonFakeSel(muons):
    return op.select(muons, lambda mu: op.AND(
        muonConePt(muons)[mu.idx] > 10,
        # self.lambda_lepton_associatedJetLessThanMediumBtag(mu),
        # op.OR(mu.mvaTTH >= 0.50, op.AND(mu.jetRelIso<0.8 , self.lambda_muon_deepJetInterpIfMvaFailed(mu)))) # will implement the second selection in the AND later, instead the following is used
        op.OR(mu.mvaTTH >= 0.50, mu.jetRelIso<0.8)
    ))


def elDef(el):
    return op.AND(
        el.pt >= 7.,
        op.abs(el.eta) < 2.5,
        op.abs(el.dxy) < 0.05,
        op.abs(el.dz) < 0.1,
        el.miniPFRelIso_all <= 0.4,
        el.sip3d < 8,
        # el.mvaNoIso_WPL, # Run3 MC doesn't have mvaNoIso_WPL for electrons
        el.lostHits <= 1
    )


def elConePt(electrons):
    return op.map(electrons, lambda lep: op.multiSwitch(
        (op.AND(op.abs(lep.pdgId) != 11, op.abs(lep.pdgId) != 13), lep.pt),
        # (op.AND(op.abs(lep.pdgId) == 11, lep.mvaTTH > 0.30), lep.pt), # run3 MC doesn't have mvaTTH for electrons
        (op.AND(op.abs(lep.pdgId) == 11), lep.pt),
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
        elConePt(electrons)[el.idx] > 10,
        op.OR(
            op.AND(op.abs(el.eta+el.deltaEtaSC)
                   <= 1.479, el.sieie < 0.011),
            op.AND(op.abs(el.eta+el.deltaEtaSC)
                   > 1.479, el.sieie < 0.030)
        ),
        el.hoe < 0.10,
        el.eInvMinusPInv > -0.04,
        # op.OR(el.mvaTTH >= 0.30, op.AND(el.jetRelIso < 0.7, el.mvaFall17V2noIso_W90)),
        # op.switch(el.mvaTTV, el.mvaTTV >= 0.30, el.mvaTTH < 0.30, self.lambda_lepton_associatedJetLessThanTightBtag(el), self.lambda_lepton_associatedJetLessThanMediumBtag(el)),
        el.lostHits == 0,
        el.convVeto,
        el.jetRelIso < 0.7
    ))


def ak4jetDef(jet):
    return op.AND(
        jet.jetId & 2,  # tight
        jet.pt > 25.,
        op.abs(jet.eta) <= 2.4
    )


def ak8jetDef(jet):
    return op.AND(
        jet.jetId & 2,  # tight
        jet.subJet1.isValid,
        jet.subJet2.isValid,
        jet.subJet1.pt > 20.,
        jet.subJet2.pt > 20.,
        op.abs(jet.subJet1.eta) <= 2.4,
        op.abs(jet.subJet2.eta) <= 2.4,
        jet.msoftdrop >= 30.,
        jet.msoftdrop <= 210.,
        jet.pt > 200.,
        op.abs(jet.eta) <= 2.4,
        jet.tau2 / jet.tau1 <= 0.75
    )
