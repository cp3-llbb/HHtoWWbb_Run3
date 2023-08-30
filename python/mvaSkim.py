
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

        if self.channel == 'DL':
            # get DL selections
            DL_boosted_ee, DL_boosted_mumu,\
            DL_boosted_emu, DL_resolved_ee,\
            DL_resolved_mumu, DL_resolved_emu = makeDLSelection(self, noSel)

        if self.channel == 'SL':
            # get SL selections
            SL_resolved, SL_resolved_e,\
            SL_resolved_mu, SL_boosted,\
            SL_boosted_e, SL_boosted_mu = makeSLSelection(self, noSel)
            

        #############################################################################
        #                                 Plots                                     #
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
                Plot.make1D("DL_boosted_nfatJet_ee", op.rng_len(self.ak8Jets), DL_boosted_ee, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet"),
            ])

        return plots
