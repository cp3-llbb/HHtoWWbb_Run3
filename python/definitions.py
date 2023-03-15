from bamboo import treefunctions as op

# Object definitions


def muonDef(mu):
    return op.AND(
        mu.pt >= 5.,
        op.abs(mu.eta) <= 2.4,
        op.abs(mu.dxy) <= 0.05,
        op.abs(mu.dz) <= 0.1,
        mu.miniPFRelIso_all <= 0.4,
        mu.sip3d <= 8,
        mu.tightId
    )


def eleDef(el):
    return op.AND(
        el.pt >= 7.,
        op.abs(el.eta) <= 2.5,
        op.abs(el.dxy) <= 0.05,
        op.abs(el.dz) <= 1.,
        el.miniPFRelIso_all <= 0.4,
        el.sip3d <= 8,
        # el.mvaNoIso_WPL,
        el.lostHits <= 1
    )


def ak8jetDef(jet):
    return op.AND(
        jet.jetId & 2,  # tight
        jet.pt > 200.,
        op.abs(jet.eta) <= 2.4
    )


def ak4jetDef(jet):
    return op.AND(
        jet.jetId & 2,  # tight
        jet.pt > 25.,
        op.abs(jet.eta) <= 2.4
    )
