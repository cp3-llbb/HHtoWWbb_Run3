import os
import re
import yaml
import logging
from itertools import chain

import utils
import definitions as defs

from bamboo import treefunctions as op
from bamboo import treedecorators as td
from bamboo.analysismodules import NanoAODModule, HistogramsModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection

logger = logging.getLogger(__name__)

JECTagDatabase = {
    "2022": {
        "MC": "Summer22_22Sep2023_V2_MC",
        "C": "Summer22_22Sep2023_RunCD_V2_DATA",
        "D": "Summer22_22Sep2023_RunCD_V2_DATA"},
    "2022EE": {
        "MC": "Summer22EE_22Sep2023_V2_MC",
        "F": "Summer22EE_22Sep2023_RunF_V2_DATA",
        "G": "Summer22EE_22Sep2023_RunG_V2_DATA"},
}

JERTagDatabase = {
    "2022": "Summer22EEPrompt22_JRV1_MC",
    "2022EE": "Summer22EEPrompt22_JRV1_MC",
}

jsonPathBase = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/"

JSONFiles = {
    "2022": {
        "AK4": jsonPathBase + "JME/2022_Summer22/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2022_Summer22/fatJet_jerc.json.gz"},
    "2022EE": {
        "AK4": jsonPathBase + "JME/2022_Summer22EE/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2022_Summer22EE/fatJet_jerc.json.gz"},
}

BTV_SF_JSONFiles = {
    "2022": jsonPathBase + "BTV/2022_Summer22/btagging.json.gz",
    "2022EE": jsonPathBase + "BTV/2022_Summer22EE/btagging.json.gz",
}


def getRunEra(sample):
    """Return run era (A/B/...) for data sample"""
    result = re.search(r'Run20..([A-Z]?)', sample)
    if result is None:
        return "MC"
    else:
        return result.group(1)


class NanoBaseHHWWbb(NanoAODModule, HistogramsModule):
    """ Base module for HH->WWbb analysis """

    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument("-c", "--channel",
                            dest="channel",
                            type=str,
                            required=True,
                            help='Channel to be selected between SL and DL')
        parser.add_argument("--mvaModels",
                            dest="mvaModels",
                            type=str,
                            help="Path to MVA models and Evaluate DNN")
        parser.add_argument("--samples", nargs='*',
                            required=True, help="Sample template YML file")
        parser.add_argument("--backend", type=str, default="dataframe",
                            help="Backend to use, 'dataframe' (default), 'lazy', or 'compiled'")
        parser.add_argument("--postprocessed", action="store_true",
                            help="Run on postprocessed NanoAOD")

    def customizeAnalysisCfg(self, analysisCfg):
        # fill sample template using JSON files
        if self.args.samples:
            eras = self.args.eras[1]
            samples = {}
            # make sure we use absolute paths as this argument will be used by the worker jobs
            self.args.samples = [os.path.abspath(p) for p in self.args.samples]
            for tmpPath in self.args.samples:
                with open(tmpPath) as f_:
                    template = yaml.load(f_, Loader=yaml.SafeLoader)
                    samples.update(utils.fillSampleTemplate(template, eras))
            analysisCfg["samples"] = samples

    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        if self.args.postprocessed:
            return self.prepare_postprocessed(tree, sample=sample, sampleCfg=sampleCfg, backend=backend)
        else:
            return self.prepare_ondemand(tree, sample=sample, sampleCfg=sampleCfg, backend=backend)

    def prepare_ondemand(self, tree, sample=None, sampleCfg=None, backend=None):
        era = sampleCfg["era"] if sampleCfg else None
        isMC = self.isMC(sample)

        # Decorate the tree
        from bamboo.treedecorators import NanoAODDescription, nanoFatJetCalc, CalcCollectionsGroups
        metName = "PuppiMET"
        nanoJetMETCalc_both = CalcCollectionsGroups(
            Jet=("pt", "mass"), changes={metName: (f"{metName}T1", f"{metName}T1Smear")},
            **{metName: ("pt", "phi")})
        nanoJetMETCalc_data = CalcCollectionsGroups(
            Jet=("pt", "mass"), changes={metName: (f"{metName}T1",)},
            **{metName: ("pt", "phi")})
        systVars = (([nanoFatJetCalc])
                    + [nanoJetMETCalc_both if isMC else nanoJetMETCalc_data])
        tree, noSel, be, lumiArgs = super().prepareTree(
            tree, sample=sample, sampleCfg=sampleCfg,
            description=NanoAODDescription.get(
                "v12", year=era[:4], isMC=isMC, systVariations=systVars),
            backend=self.args.backend or backend)

        # MC weight
        if isMC:
            noSel = noSel.refine('genWeight', weight=tree.genWeight)

        # Triggers
        self.triggersPerPrimaryDataset = {}

        def addHLTPath(PD, HLT):
            if PD not in self.triggersPerPrimaryDataset.keys():
                self.triggersPerPrimaryDataset[PD] = []
            try:
                self.triggersPerPrimaryDataset[PD].append(
                    getattr(tree.HLT, HLT))
            except AttributeError:
                print("Couldn't find branch tree.HLT.%s, cross check!" % HLT)

        # Muon
        addHLTPath('Muon', 'IsoMu24')
        addHLTPath('Muon', 'IsoMu27')
        addHLTPath('Muon', 'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8')
        # EGamma
        addHLTPath('EGamma', 'Ele32_WPTight_Gsf')
        addHLTPath('EGamma', 'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL')
        # MuonEG
        addHLTPath('MuonEG', 'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ')
        # SingleMuon
        addHLTPath('SingleMuon', 'IsoMu24')
        addHLTPath('SingleMuon', 'IsoMu27')
        # DoubleMuon
        addHLTPath('DoubleMuon', 'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8')

        if isMC:
            noSel = noSel.refine('trigger',  cut=(
                op.OR(*chain.from_iterable(self.triggersPerPrimaryDataset.values()))))
        else:
            noSel = noSel.refine('trigger', cut=makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggersPerPrimaryDataset))

        # JEC/JER
        runEra = getRunEra(sample)
        jecTag = JECTagDatabase[era]["MC" if isMC else runEra]
        smearTag = JERTagDatabase[era] if isMC else None

        cmJMEArgs = {
            "jsonFile": JSONFiles[era]["AK4"],
            "jec": jecTag,
            # "smear": smearTag,
            # "splitJER": True,
            "jesUncertaintySources": (["Total"] if isMC else None),
            "isMC": isMC,
            "backend": be,
        }
        from bamboo.analysisutils import configureJets, configureType1MET
        configureJets(tree._Jet, jetType="AK4PFPuppi", **cmJMEArgs)
        configureType1MET(
            getattr(tree, f"_{metName}T1"),
            enableSystematics=(
                (lambda v: not v.startswith("jer")) if isMC else None),
            **cmJMEArgs)
        cmJMEArgs.update({"jsonFile": JSONFiles[era]["AK8"], })
        cmJMEArgs.update({"jetAlgoSubjet": "AK4PFPuppi", })
        cmJMEArgs.update({"jecSubjet": jecTag, })
        cmJMEArgs.update({"jsonFileSubjet": JSONFiles[era]["AK4"], })
        configureJets(tree._FatJet, jetType="AK8PFPuppi", **cmJMEArgs)

        # define objects
        defs.defineObjects(self, tree)

        # btagging SF
        if isMC:
            from bamboo.scalefactors import get_bTagSF_itFit, makeBtagWeightItFit

            def btvSF(flav): return get_bTagSF_itFit(
                BTV_SF_JSONFiles[era], "particleNet", "btagDeepFlavB", flav, noSel)
            btvWeight = makeBtagWeightItFit(self.ak4Jets, btvSF)
            noSel = noSel.refine("btag", weight=btvWeight)

        # top pt reweighting
        if isMC and sample.startswith("TT"):
            def top_pt_weight(pt):
                return 0.103 * op.exp(-0.0118 * pt) - 0.000134 * pt + 0.973

            def getTopPtWeight(tree):
                lastCopy = op.select(
                    tree.GenPart, lambda p: (op.static_cast("int", p.statusFlags) >> 13) & 1)
                tops = op.select(lastCopy, lambda p: p.pdgId == 6)
                antitops = op.select(lastCopy, lambda p: p.pdgId == -6)
                weight = op.switch(op.AND(op.rng_len(tops) >= 1, op.rng_len(antitops) >= 1),
                                   op.sqrt(top_pt_weight(
                                       tops[0].pt) * top_pt_weight(antitops[0].pt)),
                                   1.)
                return weight

            logger.info(
                "Applying Top Pt reweighting everywhere, with a 'noTopPt' variation without it")

            noSel = noSel.refine("topPt", weight=op.systematic(
                getTopPtWeight(tree), noTopPt=op.c_float(1.)))

        return tree, noSel, be, lumiArgs

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Postprocess: run plotIt

        The list of plots is created if needed (from a representative file,
        this enables rerunning the postprocessing step on the results files),
        and then plotIt is executed
        """
        if not self.plotList:
            self.plotList = self.getPlotList(
                resultsdir=resultsdir, config=config)
        from bamboo.plots import Plot, DerivedPlot, CutFlowReport
        plotList_cutflowreport = [
            ap for ap in self.plotList if isinstance(ap, CutFlowReport)]
        plotList_plotIt = [ap for ap in self.plotList
                           if (isinstance(ap, Plot) or isinstance(ap, DerivedPlot))
                           and len(ap.binnings) == 1]
        eraMode, eras = self.args.eras
        if eras is None:
            eras = list(config["eras"].keys())
        if plotList_cutflowreport:
            from bamboo.analysisutils import printCutFlowReports
            printCutFlowReports(
                config, plotList_cutflowreport, workdir=workdir, resultsdir=resultsdir,
                readCounters=self.readCounters, eras=(eraMode, eras), verbose=self.args.verbose)
        if plotList_plotIt:
            from bamboo.analysisutils import writePlotIt, runPlotIt
            import os
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(
                config, plotList_plotIt, cfgName, eras=eras, workdir=workdir, resultsdir=resultsdir,
                readCounters=self.readCounters, plotDefaults=self.plotDefaults,
                vetoFileAttributes=self.__class__.CustomSampleAttributes)
            runPlotIt(
                cfgName, workdir=workdir, plotIt=self.args.plotIt, eras=(
                    eraMode, eras),
                verbose=self.args.verbose)

        from bamboo.plots import Skim
        skims = [ap for ap in self.plotList if isinstance(ap, Skim)]

        from bamboo.analysisutils import loadPlotIt
        p_config, samples, _, systematics, legend = loadPlotIt(
            config, [], eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes)

        if skims and not self.args.mvaModels:
            for skim in skims:
                pqoutname = os.path.join(resultsdir, f"{skim.name}.parquet")
                if os.path.isfile(pqoutname):
                    print(
                        f"WARNING: dataframe for skim {skim.name} already exists in {resultsdir}")
                    print("         skipping...")
                    return
                else:
                    from bamboo.root import gbl
                    import pandas as pd
                    frames = []
                    for smp in samples:
                        for cb in (smp.files if hasattr(smp, "files") else [smp]):
                            tree = cb.tFile.Get(skim.treeName)
                            if not tree:
                                print("WARNING: skim tree %s not found in file %s" % (
                                    skim.treeName, cb.tFile.GetName()))
                                print("         skipping...")
                            else:
                                N = tree.GetEntries()
                                cols = gbl.ROOT.RDataFrame(tree).AsNumpy()
                                cols["weight"] *= cb.scale
                                cols["process"] = [smp.name] * \
                                    len(cols["weight"])
                                frames.append(pd.DataFrame(cols))
                    df = pd.concat(frames)
                    df["process"] = pd.Categorical(
                        df["process"], categories=pd.unique(df["process"]))
                    df.to_parquet(pqoutname)
                    print(
                        f"Saved dataframe for skim {skim.name} to {pqoutname}")
