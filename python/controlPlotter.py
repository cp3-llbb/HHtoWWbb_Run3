
from bamboo.plots import Plot, CutFlowReport, Skim
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection

from basePlotter import NanoBaseHHWWbb
from selections import makeDLSelection, makeSLSelection
import definitions as defs

from itertools import chain

class controlPlotter(NanoBaseHHWWbb):
    """ Class to create control plots and skims"""

    def __init__(self, args):
        super(controlPlotter, self).__init__(args)
        self.channel = self.args.channel
        self.mvaModels = self.args.mvaModels

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        # call object definitions
        defs.defineObjects(self, tree)
        
        # cutflow report
        yields = CutFlowReport("yields", printInLog=True, recursive=True)
        plots.append(yields)
        
        yields.add(noSel, 'no selection')

        # Gen Weight and Triggers
        if self.is_MC:
            noSel = noSel.refine('genWeight', weight=tree.genWeight, cut=(
                op.OR(*chain.from_iterable(self.triggersPerPrimaryDataset.values()))))
        else:
            noSel = noSel.refine('trigger', cut=[makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggersPerPrimaryDataset)])
        
        yields.add(noSel, 'trigger sel.')

        if self.channel == 'DL':
            # get DL selections
            DL_boosted_ee, DL_boosted_mumu,\
            DL_boosted_emu, DL_resolved_ee,\
            DL_resolved_mumu, DL_resolved_emu = makeDLSelection(self, noSel)

            # cutflow report for DL channel
            yields.add(DL_boosted_ee, 'DL boosted ee')
            yields.add(DL_boosted_mumu, 'DL boosted mumu')
            yields.add(DL_boosted_emu, 'DL boosted emu')
            yields.add(DL_resolved_ee, 'DL resolved ee')
            yields.add(DL_resolved_mumu, 'DL resolved mumu')
            yields.add(DL_resolved_emu, 'DL resolved emu')
            
            # labels on plots
            DLboostedEE_label = defs.labeler('DL boosted EE')
            DLboostedMuMu_label = defs.labeler('DL boosted MuMu')
            DLboostedEMU_label = defs.labeler('DL boosted EMu')
            DLresolvedEE_label = defs.labeler('DL resolved EE')
            DLresolvedMuMu_label = defs.labeler('DL resolved MuMu')
            DLresolvedEMu_label = defs.labeler('DL resolved EMu')
            
            DLresolvedEEdnnCat1_label = defs.labeler('DL resolved EE DNN cat. 1')
            DLresolvedEEdnnCat2_label = defs.labeler('DL resolved EE DNN cat. 2')
            DLresolvedEEdnnCat3_label = defs.labeler('DL resolved EE DNN cat. 3')
            DLresolvedEEdnnCat4_label = defs.labeler('DL resolved EE DNN cat. 4')
            
        if self.channel == 'SL':
            # get SL selections
            SL_resolved, SL_resolved_e,\
            SL_resolved_mu, SL_boosted,\
            SL_boosted_e, SL_boosted_mu = makeSLSelection(self, noSel)

            # cutflow report for SL channel
            yields.add(SL_boosted_e, 'SL boosted e')
            yields.add(SL_boosted_mu, 'SL boosted mu')
            yields.add(SL_resolved_e, 'SL resolved e')
            yields.add(SL_resolved_mu, 'SL resolved mu')

            # labels on plots
            SLboostedE_label = defs.labeler('SL boosted E')
            SLboostedMu_label = defs.labeler('SL boosted Mu')
            SLresolvedE_label = defs.labeler('SL resolved E')
            SLresolvedMu_label = defs.labeler('SL resolved Mu')
        
        # mva variables
        mvaVars_DL_resolved = {
            "weight": noSel.weight,
            "ak4bjet1_pt": self.ak4BJets[0].pt,
            "ak4bjet1_eta": self.ak4BJets[0].eta,
            "ak4bjet1_phi": self.ak4BJets[0].phi,
            'ak4jet1_pt': self.ak4Jets[0].pt,
            'ak4jet1_eta': self.ak4Jets[0].eta,
            'ak4jet1_phi': self.ak4Jets[0].phi,
            'ak4jet2_pt': self.ak4Jets[1].pt,
            'ak4jet2_eta': self.ak4Jets[1].eta,
            'ak4jet2_phi': self.ak4Jets[1].phi,
            "leadingLepton_pt": self.tightElectrons[0].pt,
            "leadingLepton_eta": self.tightElectrons[0].eta,
            "leadingLepton_phi": self.tightElectrons[0].phi,
            "subleadingLepton_pt": self.tightElectrons[1].pt,
            "subleadingLepton_eta": self.tightElectrons[1].eta,
            "subleadingLepton_phi": self.tightElectrons[1].phi
        }
        
        mvaVars_SL_resolved = {
            "weight": noSel.weight,
            "ak4bjet1_pt": self.ak4BJets[0].pt,
            "ak4bjet1_eta": self.ak4BJets[0].eta,
            "ak4bjet1_phi": self.ak4BJets[0].phi,
            'ak4jet1_pt': self.ak4Jets[0].pt,
            'ak4jet1_eta': self.ak4Jets[0].eta,
            'ak4jet1_phi': self.ak4Jets[0].phi,
            'ak4jet2_pt': self.ak4Jets[1].pt,
            'ak4jet2_eta': self.ak4Jets[1].eta,
            'ak4jet2_phi': self.ak4Jets[1].phi,
            "leadingLepton_pt": self.tightElectrons[0].pt,
            "leadingLepton_eta": self.tightElectrons[0].eta,
            "leadingLepton_phi": self.tightElectrons[0].phi,
        }

        #############################################################################
        #                            MVA evaluation                                 #
        #############################################################################
        if self.args.mvaModels and self.channel == 'DL':
            mvaVars_DL_resolved.pop("weight", None)
            
            # import random
            # split_var = random.randint(0,1)
            split_var = 1

            if split_var == 0:
                model = self.args.mvaModels + "/model_test1_even/model.onnx"
            elif split_var == 1:
                model = self.args.mvaModels + "/model_test1_odd/model.onnx"
            else:
                print("ERROR: split_var is not 0 or 1")
                exit(1)
            
            dnn = op.mvaEvaluator(model, otherArgs = ("predictions"))
            inputs = op.array('float', *[op.c_float(val) for val in mvaVars_DL_resolved.values()])
            output = dnn(inputs)
            
            # DNN cuts
            DNNcat1 = DL_resolved_ee.refine("DNNcat1", cut = op.in_range(0.1, output[0], 0.6))
            DNNcat2 = DL_resolved_ee.refine("DNNcat2", cut = op.in_range(0.6, output[0], 0.8))
            DNNcat3 = DL_resolved_ee.refine("DNNcat3", cut = op.in_range(0.8, output[0], 0.92))
            DNNcat4 = DL_resolved_ee.refine("DNNcat4", cut = op.in_range(0.92, output[0], 1.0))
            
            yields.add(DNNcat1, 'DNNcat1')
            yields.add(DNNcat2, 'DNNcat2')
            yields.add(DNNcat3, 'DNNcat3')
            yields.add(DNNcat4, 'DNNcat4')
            
            plots.extend([
                Plot.make1D("dnn_score", output[0], DL_resolved_ee, EqBin(40, 0, 1.), xTitle="DNN Score", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_InvM_ee_DNNcat1", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DNNcat1, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLresolvedEEdnnCat1_label),
                Plot.make1D("DL_resolved_InvM_ee_DNNcat2", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DNNcat2, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLresolvedEEdnnCat2_label),
                Plot.make1D("DL_resolved_InvM_ee_DNNcat3", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DNNcat3, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLresolvedEEdnnCat3_label),
                Plot.make1D("DL_resolved_InvM_ee_DNNcat4", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DNNcat4, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLresolvedEEdnnCat4_label),
            ])

        #############################################################################
        #                                 Plots                                     #
        #############################################################################
        
        if self.channel == 'DL':
            plots.extend([
                #########################################
                #                 Skims                 #
                #########################################
                
                Skim("DL_resolved_ee", mvaVars_DL_resolved, DL_resolved_ee),
                
                #########################################
                ######                             ######
                ######       DL boosted plots      ######
                ######                             ######
                #########################################

                # number of ak8 b-jets
                Plot.make1D("DL_boosted_nfatJet_ee", op.rng_len(self.ak8Jets), DL_boosted_ee, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_nfatJet_mumu", op.rng_len(self.ak8Jets), DL_boosted_mumu, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_nfatJet_emu", op.rng_len(self.ak8Jets), DL_boosted_emu, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet", plotopts=DLboostedEMU_label),

                # fatjet pt
                Plot.make1D("DL_boosted_fatJet_pt_ee", self.ak8Jets[0].pt, DL_boosted_ee, EqBin(
                    100, 200, 800), title="pT(ak8jet)", xTitle="Fatjet p_{T} (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_fatJet_pt_mumu", self.ak8Jets[0].pt, DL_boosted_mumu, EqBin(
                    100, 200, 800), title="pT(ak8jet)", xTitle="Fatjet p_{T} (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_fatJet_pt_emu", self.ak8Jets[0].pt, DL_boosted_emu, EqBin(
                    100, 200, 800), title="pT(ak8jet)", xTitle="Fatjet p_{T} (GeV/c)", plotopts=DLboostedEMU_label),

                # subjet1 pt
                Plot.make1D("DL_boosted_subjet1_pt_ee", self.ak8Jets[0].subJet1.pt, DL_boosted_ee, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="First sub-jet p_{T} (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subjet1_pt_mumu", self.ak8Jets[0].subJet1.pt, DL_boosted_mumu, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="First sub-jet p_{T} (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subjet1_pt_emu", self.ak8Jets[0].subJet1.pt, DL_boosted_emu, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="First sub-jet p_{T} (GeV/c)", plotopts=DLboostedEMU_label),

                # subjet2 pt
                Plot.make1D("DL_boosted_subjet2_pt_ee", self.ak8Jets[0].subJet2.pt, DL_boosted_ee, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Second sub-jet p_{T} (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subjet2_pt_mumu", self.ak8Jets[0].subJet2.pt, DL_boosted_mumu, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Second sub-jet p_{T} (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subjet2_pt_emu", self.ak8Jets[0].subJet2.pt, DL_boosted_emu, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Second sub-jet p_{T} (GeV/c)", plotopts=DLboostedEMU_label),

                # fatjet eta
                Plot.make1D("DL_boosted_fatJet_eta_ee", self.ak8Jets[0].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_fatJet_eta_mumu", self.ak8Jets[0].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_fatJet_eta_emu", self.ak8Jets[0].eta, DL_boosted_emu, EqBin(
                    30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedEMU_label),

                # Invariant mass of leptons
                Plot.make1D("DL_boosted_InvM_ee", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_boosted_ee, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_InvM_mumu", op.invariant_mass(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_boosted_mumu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons  (GeV/c^{2})", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_InvM_emu", op.invariant_mass(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_boosted_emu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair  (GeV/c^{2})", plotopts=DLboostedEMU_label),

                # total charge of leptons
                Plot.make1D("DL_boosted_totalCharge_ee", op.sum(self.firstOSElEl[0].charge, self.firstOSElEl[1].charge), DL_boosted_ee, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_totalCharge_mumu", op.sum(self.firstOSMuMu[0].charge, self.firstOSMuMu[1].charge), DL_boosted_mumu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons ", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_totalCharge_emu", op.sum(self.firstOSElMu[0].charge, self.firstOSElMu[1].charge), DL_boosted_emu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair ", plotopts=DLboostedEMU_label),

                # invariant mass of subjets
                Plot.make1D("DL_boosted_InvM_jj_ee", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_ee, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_InvM_jj_mumu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_mumu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_InvM_jj_emu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_emu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedEMU_label),
                
                # leading lepton pt
                Plot.make1D("DL_boosted_leadingLepton_pt_ee", self.firstOSElEl[0].pt, DL_boosted_ee, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_leadingLepton_pt_mumu", self.firstOSMuMu[0].pt, DL_boosted_mumu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_leadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].pt, self.firstOSElMu[1].pt), DL_boosted_emu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedEMU_label),
                
                # sub-leading lepton pt
                Plot.make1D("DL_boosted_subleadingLepton_pt_ee", self.firstOSElEl[1].pt, DL_boosted_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subleadingLepton_pt_mumu", self.firstOSMuMu[1].pt, DL_boosted_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subleadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].pt, self.firstOSElMu[0].pt), DL_boosted_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedEMU_label),
                
                # leading lepton eta
                Plot.make1D("DL_boosted_leadingLepton_eta_ee", self.firstOSElEl[0].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_leadingLepton_eta_mumu", self.firstOSMuMu[0].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_leadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].eta, self.firstOSElMu[1].eta), DL_boosted_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedEMU_label),
                
                # sub-leading lepton eta
                Plot.make1D("DL_boosted_subleadingLepton_eta_ee", self.firstOSElEl[1].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subleadingLepton_eta_mumu", self.firstOSMuMu[1].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subleadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].eta, self.firstOSElMu[0].eta), DL_boosted_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedEMU_label),
                
                # DR between leading and sub-leading lepton
                Plot.make1D("DL_boosted_DR_leptons_ee", op.deltaR(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_boosted_ee, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_leptons_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_boosted_mumu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_leptons_emu", op.deltaR(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_boosted_emu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedEMU_label),

                # DR between leading lepton and ak8 jet
                Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_ee", op.deltaR(self.firstOSElEl[0].p4, self.ak8Jets[0].p4), DL_boosted_ee, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.ak8Jets[0].p4), DL_boosted_mumu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_emu", op.deltaR(op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), self.ak8Jets[0].p4), DL_boosted_emu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedEMU_label),

                # DR between subleading lepton and ak8 jet
                Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_ee", op.deltaR(self.firstOSElEl[1].p4, self.ak8Jets[0].p4), DL_boosted_ee, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_mumu", op.deltaR(self.firstOSMuMu[1].p4, self.ak8Jets[0].p4), DL_boosted_mumu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_emu", op.deltaR(op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].p4, self.firstOSElMu[0].p4), self.ak8Jets[0].p4), DL_boosted_emu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLboostedEMU_label),
                
                # number of electrons
                Plot.make1D("DL_boosted_nElectrons_ee", op.rng_len(self.tightElectrons), DL_boosted_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_nElectrons_mumu", op.rng_len(self.tightElectrons), DL_boosted_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_nElectrons_emu", op.rng_len(self.tightElectrons), DL_boosted_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEMU_label),
                
                # number of muons
                Plot.make1D("DL_boosted_nMuons_ee", op.rng_len(self.tightMuons), DL_boosted_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_nMuons_mumu", op.rng_len(self.tightMuons), DL_boosted_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_nMuons_emu", op.rng_len(self.tightMuons), DL_boosted_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEMU_label),
                
                #########################################
                ######                             ######
                ######      DL resolved plots      ######
                ######                             ######
                #########################################
                
                # number of ak4 bjets
                Plot.make1D("DL_resolved_nAK4bJets_ee", op.rng_len(self.ak4BJets), DL_resolved_ee, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nAK4bJets_mumu", op.rng_len(self.ak4BJets), DL_resolved_mumu, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nAK4bJets_emu", op.rng_len(self.ak4BJets), DL_resolved_emu, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolvedEMu_label),
                
                # ak4 bjet pt
                Plot.make1D("DL_resolved_ak4BJet_pt_ee", self.ak4BJets[0].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_ak4BJet_pt_mumu", self.ak4BJets[0].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_ak4BJet_pt_emu", self.ak4BJets[0].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolvedEMu_label),

                # ak4 bjet eta
                Plot.make1D("DL_resolved_ak4BJet_eta_ee", self.ak4BJets[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_ak4BJet_eta_mumu", self.ak4BJets[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_ak4BJet_eta_emu", self.ak4BJets[0].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolvedEMu_label),

                # number of ak4 jets
                Plot.make1D("DL_resolved_nak4Jets_ee", op.rng_len(self.ak4Jets), DL_resolved_ee, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nak4Jets_mumu", op.rng_len(self.ak4Jets), DL_resolved_mumu, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nak4Jets_emu", op.rng_len(self.ak4Jets), DL_resolved_emu, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolvedEMu_label),
                
                # leading jet pt
                Plot.make1D("DL_resolved_leadingJet_pt_ee", self.ak4Jets[0].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingJet_pt_mumu", self.ak4Jets[0].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingJet_pt_emu", self.ak4Jets[0].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolvedEMu_label),
                
                # leading jet eta
                Plot.make1D("DL_resolved_leadingJet_eta_ee", self.ak4Jets[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingJet_eta_mumu", self.ak4Jets[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingJet_eta_emu", self.ak4Jets[0].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolvedEMu_label),
                
                # sub-leading jet pt
                Plot.make1D("DL_resolved_subleadingJet_pt_ee", self.ak4Jets[1].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingJet_pt_mumu", self.ak4Jets[1].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingJet_pt_emu", self.ak4Jets[1].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolvedEMu_label),

                # sub-leading jet eta
                Plot.make1D("DL_resolved_subleadingJet_eta_ee", self.ak4Jets[1].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingJet_eta_mumu", self.ak4Jets[1].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingJet_eta_emu", self.ak4Jets[1].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolvedEMu_label),
                
                # DR between leading and sub-leading jet
                Plot.make1D("DL_resolved_DR_jets_ee", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_jets_mumu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_jets_emu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolvedEMu_label),
                
                # Invariant mass of leptons
                Plot.make1D("DL_resolved_InvM_ee", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_resolved_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_InvM_mumu", op.invariant_mass(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_resolved_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_InvM_emu", op.invariant_mass(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_resolved_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (GeV/c^{2})", plotopts=DLresolvedEMu_label),

                # total charge of leptons
                Plot.make1D("DL_resolved_totalCharge_ee", op.sum(self.firstOSElEl[0].charge, self.firstOSElEl[1].charge), DL_resolved_ee, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_totalCharge_mumu", op.sum(self.firstOSMuMu[0].charge, self.firstOSMuMu[1].charge), DL_resolved_mumu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_totalCharge_emu", op.sum(self.firstOSElMu[0].charge, self.firstOSElMu[1].charge), DL_resolved_emu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair", plotopts=DLresolvedEMu_label),

                # leading lepton pt
                Plot.make1D("DL_resolved_leadingLepton_pt_ee", self.firstOSElEl[0].pt, DL_resolved_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingLepton_pt_mumu", self.firstOSMuMu[0].pt, DL_resolved_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].pt, self.firstOSElMu[1].pt), DL_resolved_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolvedEMu_label),
                
                # sub-leading lepton pt
                Plot.make1D("DL_resolved_subleadingLepton_pt_ee", self.firstOSElEl[1].pt, DL_resolved_ee, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingLepton_pt_mumu", self.firstOSMuMu[1].pt, DL_resolved_mumu, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].pt, self.firstOSElMu[0].pt), DL_resolved_emu, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolvedEMu_label),
                
                # leading lepton eta
                Plot.make1D("DL_resolved_leadingLepton_eta_ee", self.firstOSElEl[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingLepton_eta_mumu", self.firstOSMuMu[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].eta, self.firstOSElMu[1].eta), DL_resolved_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolvedEMu_label),
                
                # sub-leading lepton eta
                Plot.make1D("DL_resolved_subleadingLepton_eta_ee", self.firstOSElEl[1].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingLepton_eta_mumu", self.firstOSMuMu[1].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].eta, self.firstOSElMu[0].eta), DL_resolved_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolvedEMu_label),
                
                # DR between leading and sub-leading lepton
                Plot.make1D("DL_resolved_DR_leptons_ee", op.deltaR(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_leptons_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_leptons_emu", op.deltaR(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolvedEMu_label),

                # DR between leading lepton and ak4 b jet
                Plot.make1D("DL_resolved_DR_leadingleptonANDak4bjet_ee", op.deltaR(self.firstOSElEl[0].p4, self.ak4BJets[0].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_leadingleptonANDak4bjet_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.ak4BJets[0].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_leadingleptonANDak4bjet_emu", op.deltaR(op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), self.ak4BJets[0].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolvedEMu_label),

                # DR between sub-leading lepton and ak4 b jet
                Plot.make1D("DL_resolved_DR_subleadingleptonANDak4bjet_ee", op.deltaR(self.firstOSElEl[1].p4, self.ak4BJets[0].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_subleadingleptonANDak4bjet_mumu", op.deltaR(self.firstOSMuMu[1].p4, self.ak4BJets[0].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_subleadingleptonANDak4bjet_emu", op.deltaR(op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].p4, self.firstOSElMu[0].p4), self.ak4BJets[0].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLresolvedEMu_label),
                
                # number of electrons
                Plot.make1D("DL_resolved_nElectrons_ee", op.rng_len(self.tightElectrons), DL_resolved_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nElectrons_mumu", op.rng_len(self.tightElectrons), DL_resolved_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nElectrons_emu", op.rng_len(self.tightElectrons), DL_resolved_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEMu_label),
                
                # number of muons
                Plot.make1D("DL_resolved_nMuons_ee", op.rng_len(self.tightMuons), DL_resolved_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nMuons_mumu", op.rng_len(self.tightMuons), DL_resolved_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nMuons_emu", op.rng_len(self.tightMuons), DL_resolved_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEMu_label),
            ])
        if self.channel == "SL":
            plots.extend([
                
                #########################################
                #                 Skims                 #
                #########################################
                
                Skim("SL_resolved_e", mvaVars_SL_resolved, SL_resolved_e),
                
                #########################################
                ######                             ######
                ######       SL boosted plots      ######
                ######                             ######
                #########################################
                
                # number of fat b-jets
                Plot.make1D("SL_boosted_nfatJet_e", op.rng_len(self.ak8BJets), SL_boosted_e, EqBin(
                    10, 0, 10), title="N(ak8bjet)", xTitle="Number of fat b-jet", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_nfatJet_mu", op.rng_len(self.ak8BJets), SL_boosted_mu, EqBin(
                    10, 0, 10), title="N(ak8bjet)", xTitle="Number of fat b-jet", plotopts=SLboostedMu_label),

                # fatjet pt
                Plot.make1D("SL_boosted_fatJet_pt_e", self.ak8BJets[0].pt, SL_boosted_e, EqBin(
                    400, 200, 1000), title="pT(j)", xTitle="Fat b-jet p_{T} (GeV/c)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_fatJet_pt_mu", self.ak8BJets[0].pt, SL_boosted_mu, EqBin(
                    400, 200, 1000), title="pT(j)", xTitle="Fat b-jet p_{T} (GeV/c)", plotopts=SLboostedMu_label),
                
                # subjet1 pt
                Plot.make1D("SL_boosted_subjet1_pt_e", self.ak8BJets[0].subJet1.pt, SL_boosted_e, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="Subjet 1 p_{T} (GeV/c)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_subjet1_pt_mu", self.ak8BJets[0].subJet1.pt, SL_boosted_mu, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="Subjet 1 p_{T} (GeV/c)", plotopts=SLboostedMu_label),
                
                # subjet2 pt
                Plot.make1D("SL_boosted_subjet2_pt_e", self.ak8BJets[0].subJet2.pt, SL_boosted_e, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Subjet 2 p_{T} (GeV/c)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_subjet2_pt_mu", self.ak8BJets[0].subJet2.pt, SL_boosted_mu, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Subjet 2 p_{T} (GeV/c)", plotopts=SLboostedMu_label),
                
                # ak8jet eta
                Plot.make1D("SL_boosted_fatJet_eta_e", self.ak8BJets[0].eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(j)", xTitle="eta(j)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_fatJet_eta_mu", self.ak8BJets[0].eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(j)", xTitle="eta(j)", plotopts=SLboostedMu_label),
                
                # subjet1 eta
                Plot.make1D("SL_boosted_subjet1_eta_e", self.ak8BJets[0].subJet1.eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(subjet1)", xTitle="Subjet 1 \eta", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_subjet1_eta_mu", self.ak8BJets[0].subJet1.eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(subjet1)", xTitle="Subjet 1 \eta", plotopts=SLboostedMu_label),
                
                # subjet2 eta
                Plot.make1D("SL_boosted_subjet2_eta_e", self.ak8BJets[0].subJet2.eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(subjet2)", xTitle="Subjet 2 \eta", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_subjet2_eta_mu", self.ak8BJets[0].subJet2.eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(subjet2)", xTitle="Subjet 2 \eta", plotopts=SLboostedMu_label),
                
                # Invariant mass of subjets
                Plot.make1D("SL_boosted_InvM_jj_e", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), SL_boosted_e, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_InvM_jj_mu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), SL_boosted_mu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=SLboostedMu_label),
                
                # lepton pt
                Plot.make1D("SL_boosted_lepton_pt_e", self.tightElectrons[0].pt, SL_boosted_e, EqBin(
                    100, 0., 300.), title="lepton pT", xTitle="p_{T} of the lepton (GeV/c)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_lepton_pt_mu", self.tightMuons[0].pt, SL_boosted_mu, EqBin(
                    100, 0., 300.), title="lepton pT", xTitle="p_{T} of the lepton (GeV/c)", plotopts=SLboostedMu_label),
                
                # lepton eta
                Plot.make1D("SL_boosted_lepton_eta_e", self.tightElectrons[0].eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="lepton eta", xTitle="lepton \eta", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_lepton_eta_mu", self.tightMuons[0].eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="lepton eta", xTitle="lepton \eta", plotopts=SLboostedMu_label),
                
                # DR between lepton and ak8 jet
                Plot.make1D("SL_boosted_DR_leptonANDak8bjet_e", op.deltaR(self.tightElectrons[0].p4, self.ak8Jets[0].p4), SL_boosted_e, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_DR_leptonANDak8bjet_mu", op.deltaR(self.tightMuons[0].p4, self.ak8Jets[0].p4), SL_boosted_mu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=SLboostedMu_label),


                #########################################
                ######                             ######
                ######      SL resolved plots      ######
                ######                             ######
                #########################################

                # number of ak4 jets
                Plot.make1D("SL_resolved_nJets_e", op.rng_len(self.ak4BJets), SL_resolved_e, EqBin(
                    15, 0., 15.), xTitle="Number of b-jets", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_nJets_mu", op.rng_len(self.ak4BJets), SL_resolved_mu, EqBin(
                    15, 0., 15.), xTitle="Number of b-jets", plotopts=SLresolvedMu_label),
                
                # leading jet pt
                Plot.make1D("SL_resolved_leadingJet_pt_e", self.ak4BJets[0].pt, SL_resolved_e, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_leadingJet_pt_mu", self.ak4BJets[0].pt, SL_resolved_mu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)", plotopts=SLresolvedMu_label),
                
                # leading jet eta
                Plot.make1D("SL_resolved_leadingJet_eta_e", self.ak4BJets[0].eta, SL_resolved_e, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="B-jet \eta", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_leadingJet_eta_mu", self.ak4BJets[0].eta, SL_resolved_mu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="eta(j1)", plotopts=SLresolvedMu_label),
                
                # sub-leading jet pt
                Plot.make1D("SL_resolved_subleadingJet_pt_e", self.ak4Jets[1].pt, SL_resolved_e, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_subleadingJet_pt_mu", self.ak4Jets[1].pt, SL_resolved_mu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)", plotopts=SLresolvedMu_label),
                
                # sub-leading jet eta
                Plot.make1D("SL_resolved_subleadingJet_eta_e", self.ak4Jets[1].eta, SL_resolved_e, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="eta(j2)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_subleadingJet_eta_mu", self.ak4Jets[1].eta, SL_resolved_mu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="eta(j2)", plotopts=SLresolvedMu_label),
                
                # DR between leading and sub-leading jet
                Plot.make1D("SL_resolved_DR_jets_e", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), SL_resolved_e, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_DR_jets_mu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), SL_resolved_mu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)", plotopts=SLresolvedMu_label),

                # DR between  lepton and ak4 b jet
                Plot.make1D("SL_resolved_DR_leptonANDak4bjet_e", op.deltaR(self.tightElectrons[0].p4, self.ak4BJets[0].p4), SL_resolved_e, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(lepton, ak8bjet)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_DR_leptonANDak4bjet_mu", op.deltaR(self.tightMuons[0].p4, self.ak4BJets[0].p4), SL_resolved_mu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(lepton, ak8bjet)", plotopts=SLresolvedMu_label),
                
                # lepton pt
                Plot.make1D("SL_resolved_lepton_pt_e", self.tightElectrons[0].pt, SL_resolved_e, EqBin(
                    100, 0., 300.), title="lepton pT", xTitle="p_{T} of the lepton (GeV/c)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_lepton_pt_mu", self.tightMuons[0].pt, SL_resolved_mu, EqBin(
                    100, 0., 300.), title="lepton pT", xTitle="p_{T} of the lepton (GeV/c)", plotopts=SLresolvedMu_label),
                
                # lepton eta
                Plot.make1D("SL_resolved_lepton_eta_e", self.tightElectrons[0].eta, SL_resolved_e, EqBin(
                    30, -3, 3), title="lepton pT", xTitle="p_{T} of the lepton (GeV/c)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_lepton_eta_mu", self.tightMuons[0].eta, SL_resolved_mu, EqBin(
                    30, -3, 3), title="lepton pT", xTitle="p_{T} of the lepton (GeV/c)", plotopts=SLresolvedMu_label),
            ])

        return plots
