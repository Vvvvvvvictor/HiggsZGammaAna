CENTRAL_WEIGHT = "weight_central" # name of the central weight branch
NOMINAL_TAG = "nominal" # name of the nominal events (relevant when we have multiple sets of events corresponding to systematics with independent collections)

#https://indico.cern.ch/event/560226/contributions/2277448/attachments/1324704/1988050/wgm_vfp_change_ebutz.pdf
#pre-VFP  runs: 273150-278800 lumi: 19480.4566773 /pb
#post-VFP runs: 278801-284044 lumi: 16755.0362868 /pb

LUMI = {
    "2016" : 36.31,
    "2016preVFP" : 19.5, 
    "2016postVFP" : 16.8, 
    "2017" : 41.48,
    "2018" : 59.83,
    "2022" : 35.18,
    "2022preEE": 7.98,
    "2022postEE": 26.67,
    "2023" : 27.14,
    "2023preBPix": 17.79,
    "2023postBPix": 9.45
}

GOLDEN_JSON = {
    "2016" : "metadata/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
    "2016preVFP" : "metadata/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
    "2016postVFP" : "metadata/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
    "2017" : "metadata/golden_json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
    "2018" : "metadata/golden_json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
    "2022" : "metadata/golden_json/Cert_Collisions2022_355100_362760_Golden.json",
    "2022preEE" : "metadata/golden_json/Cert_Collisions2022_355100_362760_Golden.json",
    "2022postEE" : "metadata/golden_json/Cert_Collisions2022_355100_362760_Golden.json",
    "2023" : "metadata/golden_json/Cert_Collisions2023_366442_370790_Golden.json",
    "2023preBPix" : "metadata/golden_json/Cert_Collisions2023_366442_370790_Golden.json",
    "2023postBPix" : "metadata/golden_json/Cert_Collisions2023_366442_370790_Golden.json"
}

# nanoAOD branches to always include
BRANCHES = {
    "data" : {
        "2016" : [
            "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"
        ],
        "2016postVFP" : [
            "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"
        ],
        "2016preVFP" : [
            "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"
        ],
        "2017" : [
            "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"
        ],
        "2018" : [
            "HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90"
        ],
        "2022" : [
            "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"
        ],
        "2022preEE" : [],
        "2022postEE" : [],
        "2023" : [
            "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"
        ],
        "2023preBPix" : [],
        "2023postBPix" : [],
        "any" : ["event", "run", "luminosityBlock"]
    },
    "mc" : {
        "2016APV" : [],
        "2016" : [],
        "2016preVFP" : [],
        "2016postVFP" : [],
        "2017" : [],
        "2018" : [],
        "2022" : [],
        "2022preEE" : [],
        "2022postEE" : [],
        "2023" : [],
        "2023preBPix" : [],
        "2023postBPix" : [],
        "any" : ["event", "run", "luminosityBlock"]
    }
}
