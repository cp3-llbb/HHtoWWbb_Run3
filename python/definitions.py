from bamboo import treefunctions as op

# Lepton Lambda funtions


def lambda_hasAssociatedJet(lep): return lep.jet.idx != -1


def lambda_lepton_associatedJetLessThanMediumBtag(lep): return op.OR(op.NOT(lambda_hasAssociatedJet(lep)),
                                                                     lep.jet.btagDeepFlavB <= 0.2770)  # 2018 value


def lambda_muon_x(mu): return op.min(
    op.max(0., (0.9*mu.pt*(1+mu.jetRelIso))-20.)/(45.-20.), 1.)


def lambda_muon_btagInterpolation(mu): return lambda_muon_x(
    mu)*0.0494 + (1-lambda_muon_x(mu))*0.2770  # 2018 values


def lambda_muon_deepJetInterpIfMvaFailed(mu): return op.OR(op.NOT(
    lambda_hasAssociatedJet(mu)), mu.jet.btagDeepFlavB < lambda_muon_btagInterpolation(mu))


def lambda_lepton_associatedJetLessThanMediumBtag(lep): return op.OR(op.NOT(lambda_hasAssociatedJet(lep)),
                                                                     lep.jet.btagDeepFlavB <= 0.2770)  # 2018 value


def lambda_lepton_associatedJetLessThanTightBtag(lep): return op.OR(op.NOT(lambda_hasAssociatedJet(lep)),
                                                                    lep.jet.btagDeepFlavB <= 0.7264)  # 2018 value

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
        lambda_lepton_associatedJetLessThanMediumBtag(mu),
        op.OR(mu.mvaTTH > 0.50, op.AND(mu.jetRelIso < 0.8, lambda_muon_deepJetInterpIfMvaFailed(mu))))
    )


def elDef(el):
    return op.AND(
        el.pt > 7.,
        op.abs(el.eta) < 2.5,
        op.abs(el.dxy) < 0.05,
        op.abs(el.dz) < 0.1,
        el.sip3d < 8,
        el.miniPFRelIso_all < 0.4,
        el.mvaNoIso > 0.5,  # check WPL
        el.lostHits <= 1
    )


def elConePt(electrons):
    return op.map(electrons, lambda lep: op.multiSwitch(
        (op.AND(op.abs(lep.pdgId) != 11, op.abs(lep.pdgId) != 13), lep.pt),
        (op.AND(op.abs(lep.pdgId) == 11, lep.mvaTTH > 0.30), lep.pt),
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
        op.OR(el.mvaTTH >= 0.30, op.AND(
            el.jetRelIso < 0.7, el.mvaNoIso_WP90)),
        op.switch(
            el.mvaTTH < 0.30,
            lambda_lepton_associatedJetLessThanTightBtag(el),
            lambda_lepton_associatedJetLessThanMediumBtag(el)),
        el.lostHits == 0,
        el.convVeto,
        el.jetRelIso < 0.7
    ))


def ak4jetDef(jet):
    return op.AND(
        jet.jetId & 2,  # tight
        jet.pt > 25.,
        op.abs(jet.eta) < 2.4
    )


def ak8jetDef(jet):
    return op.AND(
        jet.pt > 200.,
        op.abs(jet.eta) < 2.4,
        jet.jetId & 2,  # tight
        jet.subJet1.isValid,
        jet.subJet2.isValid,
        jet.subJet1.pt > 20.,
        jet.subJet2.pt > 20.,
        op.OR(op.AND(jet.subJet1.pt > 30., jet.subJet1.btagDeepB > 0.4184),
              op.AND(jet.subJet2.pt > 30., jet.subJet2.btagDeepB > 0.4184)),
        op.abs(jet.subJet1.eta) < 2.4,
        op.abs(jet.subJet2.eta) < 2.4,
        op.AND(jet.msoftdrop > 30., jet.msoftdrop < 210.),
        jet.tau2 / jet.tau1 < 0.75
    )


def tauDef(taus):
    return op.select(taus, lambda tau: op.AND(
        tau.pt > 20.,
        op.abs(tau.p4.Eta()) < 2.3,
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
