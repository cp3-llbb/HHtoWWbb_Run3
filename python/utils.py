from copy import deepcopy


def fillSampleTemplate(template, selEras=None):
    outTemplate = {}

    # Expand eras
    for name, sample in template.items():
        if "dbs" in sample:
            for era, das in sample["dbs"].items():
                era = str(era)
                if selEras is not None and era not in selEras:
                    continue
                thisSample = deepcopy(sample)
                if "syst" in thisSample:
                    syst, nom = thisSample["syst"]
                    newName = f"{nom}__{era}__{syst}"
                    thisSample["syst"][1] = f"{name}__{era}"
                else:
                    newName = f"{name}__{era}"
                thisSample["db"] = das
                thisSample["era"] = era
                thisSample.pop("dbs")
                outTemplate[newName] = thisSample
        else:
            outTemplate[name] = sample

    return outTemplate


def labeler(label):
    return {'labels': [{'text': label, 'position': [0.235, 0.9], 'size': 24}]}
