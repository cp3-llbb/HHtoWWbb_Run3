
from bamboo.analysismodules import NanoAODHistoModule

from bamboo.treedecorators import NanoAODDescription
description_CMSRun1OpenData_ROOT_H4ltutorial = NanoAODDescription(groups=["PV_"], collections=["nMuon", "nElectron", "nFatJet", "nJet"])

from bamboo.plots import Plot, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

class NanoBaseHHWWbb(NanoAODHistoModule):
    def mergeCounters(self, outF, infileNames, sample=None):
        pass ## disabled because the test file is not a full NanoAOD
    def isMC(self, sampleName):
        return sampleName.split("_")[0] not in ("DoubleMuParked", "DoubleElectron")
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        return super(NanoBaseHHWWbb, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg, description=description_CMSRun1OpenData_ROOT_H4ltutorial)
    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        muons = op.sort(op.select(tree.Muon, lambda mu : op.AND(
            mu.pt >= 5.,
            op.abs(mu.eta) <= 2.4,
            op.abs(mu.dxy) <= 0.05,
            op.abs(mu.dz) <= 0.1,
            mu.miniPFRelIso_all <= 0.4,
            mu.sip3d <= 8,
            mu.tightId
            )), lambda mu : -mu.pt)
        electrons = op.sort(op.select(tree.Electron, lambda el : op.AND(
            el.pt >= 7.,
            op.abs(el.eta) <= 2.5,
            op.abs(el.dxy) <= 0.05,
            op.abs(el.dz) <= 1.,
            el.miniPFRelIso_all <= 0.4,
            el.sip3d <= 8,
            # el.mvaNoIso_WPL,
            el.lostHits <=1
            )), lambda el : -el.pt)
        
        ########################### AK8 Jets ###########################
        ak8Jets = op.sort(op.select(tree.FatJet, lambda j : op.AND(
            j.jetId & 2, # tight
            j.pt > 200.,
            op.abs(j.eta) <= 2.4
        )), lambda jet : -jet.pt)

        ak8BJets = op.select(ak8Jets, lambda fatjet : fatjet.btagDeepB > 0.4184) # 2018 WP

        ########################### AK4 Jets ###########################
        ak4Jets = op.sort(op.select(tree.Jet, lambda j : op.AND(
            j.jetId & 2, # tight
            j.pt >= 25.,
            op.abs(j.eta) < 2.4
        )), lambda jet : -jet.pt)

        ak4BJets = op.select(ak4Jets, lambda jet : jet.btagDeepB > 0.2770) # 2018 WP

        plots += self.controlPlots_2l(noSel, muons, electrons, ak4Jets, ak4BJets)
        return plots

    def controlPlots_2l(self, noSel, muons, electrons, jets, bjets):
        plots = [
            Plot.make1D("nEl", op.rng_len(electrons), noSel, EqBin(10, 0., 10.), xTitle="Number of good electrons"),
            Plot.make1D("nMu", op.rng_len(muons), noSel, EqBin(10, 0., 10.), xTitle="Number of good muons"),
            Plot.make1D("nJet", op.rng_len(jets), noSel, EqBin(10, 0., 10.), xTitle="Number of good jets"),
            ]

        hasOSElEl = noSel.refine("hasOSElEl", cut=[ op.rng_len(electrons) >= 2, # op.rng_len(muons) >= 1, # for MuonEG dataset
            electrons[0].charge != electrons[1].charge, electrons[0].pt > 20., electrons[1].pt > 10. ])

        hasJets = hasOSElEl.refine("hasJets", cut=[ op.rng_len(jets) >= 2, op.rng_len(bjets) >= 2 ])

        plots.append(Plot.make1D("massZto2e", op.invariant_mass(electrons[0].p4, electrons[1].p4),
            hasOSElEl, EqBin(120, 40., 120.), title="mass of Z to 2e",
            xTitle="Invariant Mass of Nelectrons=2 (in GeV/c^2)"))

        plots.append(Plot.make1D("massZto2e_hasJets", op.invariant_mass(electrons[0].p4, electrons[1].p4),
            hasJets, EqBin(120, 40., 120.), title="mass of Z to 2e",
            xTitle="Invariant Mass of Nelectrons=2 (in GeV/c^2)"))
        hasOSMuMu = noSel.refine("hasOSMuMu", cut=[ op.rng_len(muons) >= 2, op.rng_len(jets) >= 2, op.rng_len(bjets) >= 2,
            muons[0].charge != muons[1].charge, muons[0].pt > 20., muons[1].pt > 10. ])
        plots.append(Plot.make1D("massZto2mu", op.invariant_mass(muons[0].p4, muons[1].p4),
            hasOSMuMu, EqBin(120, 40., 120.), title="mass of Z to 2mu",
            xTitle="Invariant Mass of Nmuons=2 (in GeV/c^2)"))
        return plots