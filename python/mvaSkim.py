
from bamboo.plots import Plot, Skim
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

import definitions as defs

from basePlotter import NanoBaseHHWWbb
from selections import makeDLSelection, makeSLSelection


class mvaSkimmer(NanoBaseHHWWbb):
    """HH->WWbb control plots"""

    def __init__(self, args):
        super(mvaSkimmer, self).__init__(args)
        self.channel = self.args.channel
        self.mvaSkim = self.args.mvaSkim

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        # define objects
        defs.defineObjects(self, tree)

        def labeler(label):
            return {'labels': [{'text': label, 'position': [0.23, 0.87], 'size': 25}]}

        if self.channel == 'DL':
            # get DL selections
            DL_boosted_ee, DL_boosted_mumu,\
            DL_boosted_emu, DL_resolved_ee,\
            DL_resolved_mumu, DL_resolved_emu = makeDLSelection(self, noSel)
            
            DLresolvedEE_label = labeler('DL resolved EE')
            DLresolvedMuMu_label = labeler('DL resolved MuMu')
            DLresolvedEMu_label = labeler('DL resolved EMu')

        if self.channel == 'SL':
            # get SL selections
            SL_resolved, SL_resolved_e,\
            SL_resolved_mu, SL_boosted,\
            SL_boosted_e, SL_boosted_mu = makeSLSelection(self, noSel)
            

        #############################################################################
        #                                 Skim                                      #
        #############################################################################
        if self.args.mvaSkim and self.channel == 'DL':
            mvaVars_DL_resolved_ee = {
                "weight": noSel.weight,
                "ak4bjet1_pt": self.ak4BJets[0].pt,
                "ak4bjet1_eta": self.ak4BJets[0].eta,
                "ak4bjet1_phi": self.ak4BJets[0].phi,
                "leadingLepton_pt": self.tightElectrons[0].pt,
                "leadingLepton_eta": self.tightElectrons[0].eta,
                "leadingLepton_phi": self.tightElectrons[0].phi,
                "subleadingLepton_pt": self.tightElectrons[1].pt,
                "subleadingLepton_eta": self.tightElectrons[1].eta,
                "subleadingLepton_phi": self.tightElectrons[1].phi
            }
            
            plots.extend([
                Skim("DL_resolved", mvaVars_DL_resolved_ee, DL_resolved_ee),
                # Invariant mass of leptons
                Plot.make1D("DL_resolved_InvM_ee", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_resolved_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_InvM_mumu", op.invariant_mass(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_resolved_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_InvM_emu", op.invariant_mass(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_resolved_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (GeV/c^{2})", plotopts=DLresolvedEMu_label),
            ])

        #############################################################################
        #                            DNN evaluation                                 #
        #############################################################################
        
        if self.args.mvaEval and self.channel == 'DL':
            
            DNN_model_even = self.resultsdir + '/model_even.onnx'
            DNN_model_odd = self.resultsdir + '/model_odd.onnx'
            
            
        return plots
