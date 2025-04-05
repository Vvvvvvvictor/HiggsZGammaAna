import uproot

file_path = "root://cms-xrd-global.cern.ch//store/data/Run2023D/Muon0/NANOAOD/22Sep2023_v1-v1/50000/6ee1a50c-efb1-4b1c-8b3c-9a5b01376064.root"

try:
    with uproot.open(file_path) as f:
        print(f.keys())  # 列出 TTree
except Exception as e:
    print("Failed to open file with uproot:", e)
