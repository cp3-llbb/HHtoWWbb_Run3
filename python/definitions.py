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
