import json
import math
from pathlib import Path
from datetime import datetime

# 在這裡填入四個 fold 的 significance；請將 None 改成數字 (float)
FOLD_SIGNIFICANCE = {
    "ggF": [0.49, 0.61, 0.79, 0.45],  # 例: [2.31, 2.45, 2.28, 2.52]
    "VBF": [0.85, 0.77, 0.61, 0.30],  # 例: [0.85, 0.90, 0.78, 0.88]
}

def check_filled(values):
    return all(v is not None for v in values)

def compute_stats(values):
    # sample mean
    n = len(values)
    mean = sum(values) / n
    # sample standard deviation (n-1)
    if n > 1:
        var = sum((v - mean) ** 2 for v in values) / (n - 1)
        std = math.sqrt(var)
    else:
        std = 0.0
    return mean, std

def main():
    # 檢查是否所有 None 都已替換成數字
    incomplete = {k: v for k, v in FOLD_SIGNIFICANCE.items() if not check_filled(v)}
    if incomplete:
        print("尚有未填入的 fold significance：")
        for k, arr in incomplete.items():
            print(f"  {k}: {arr}")
        print("請編輯此檔案，將 None 改成對應的數值後再執行。")
        return

    output = {}
    for proc, arr in FOLD_SIGNIFICANCE.items():
        mean, std = compute_stats(arr)
        output[proc] = {
            "folds": arr,
            "mean": mean,
            "std": std,
        }

    output["timestamp"] = datetime.utcnow().isoformat() + "Z"

    # 將輸出目錄改為 Path 物件，並以腳本所在資料夾建立 plots 子資料夾
    out_dir = "/afs/cern.ch/work/p/pelai/HZgamma/HiggsZGammaAna/hzgml/plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "fold_significance_uncertainty.json"

    with out_path.open("w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False, indent=2)

    print(f"已寫入: {out_path}")
    for k in ("ggF", "VBF"):
        print(f"{k}: mean={output[k]['mean']:.4f}, std={output[k]['std']:.4f}")

if __name__ == "__main__":
    main()
