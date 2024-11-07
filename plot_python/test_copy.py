import matplotlib.pyplot as plt
import uproot
import os
import numpy as np
import mplhep as hep
plt.style.use(hep.style.CMS)
from pdb import set_trace

# Settings
file_path = "/eos/home-j/jiehan/root"
type_dirs = ["", "_run2_gdr_4", "_run2_jet25"]
signal_samples = ["ggH", "VBF", "WminusH", "WplusH", "ZH", "ttH"]
background_samples = ["ZGToLLG", "DYJetsToLL", "EWKZ2J", "ZG2JToG2L2J", "WGToLNuG", "TT", "TTGJets", "TGJets", "WW", "WZ", "ZZ", "ttZJets", "ttWJets"]
trees = ["zero_to_one_jet", "two_jet", "VH", "ZH", "ttH_lep", "ttH_had"]


# Load the root file
sig_yields, bkg_yields = [], []

for type in type_dirs:
    for sig in signal_samples:
        sig_yield = np.zeros(len(trees))
        for file in os.listdir(f"{file_path}/skimmed_ntuples{type}/{sig}_M125"):
            if file.endswith(".root"):
                temp_yields = np.array([])
                print(f"Opening {file_path}/skimmed_ntuples{type}/{sig}_M125/{file}")
                f = uproot.open(f"{file_path}/skimmed_ntuples{type}/{sig}_M125/{file}")
                for chan in trees:
                    data = f[chan].arrays(["H_mass", "weight"], library="pd")
                    temp_yields = np.append(temp_yields, np.sum(data.query("H_mass > 100")["weight"]))
                sig_yield += temp_yields
        sig_yields.append(sig_yield)

# 假设数据集是6个channel，每个channel有6个样本的yields
channels = np.arange(len(trees))*4  # 6个channel作为x轴
n_samples = len(signal_samples)  # 每个channel有6个样本
n_datasets = len(type_dirs)  # 三类数据集
width = 0.4  # 每个柱状图的宽度
gap_between_sets = 0.1  # 每个数据集之间的间距

# 假设我们有三个数据集，每个有6个channel和6个samples的数据
# 这里随机生成数据，实际使用时替换为你的数据
# np.random.seed(0)
# yields = np.random.rand(n_datasets, n_samples, len(channels))
temp_yields = np.array(sig_yields).reshape(n_datasets, n_samples, len(channels))
# 归一化最后一维
yields = temp_yields / temp_yields[1].sum(axis=0)

# 定义不同samples的颜色
colors = ['blue', 'green', 'red', 'orange', 'purple', 'brown']

# 创建图形和子图
fig, ax = plt.subplots()

# 创建堆叠柱状图
for i, dataset in enumerate(yields):
    # 计算每个数据集的x坐标，留出一定的间距
    x = channels + i * (width + gap_between_sets)
    
    # 堆叠样本的柱状图
    bottom = np.zeros(len(channels))  # 初始化底部为0
    for j in range(n_samples):
        ax.bar(x, dataset[j], width, bottom=bottom, label=signal_samples[j] if i == 0 else "", color=colors[j])
        bottom += dataset[j]  # 更新底部

# 标注数据集类型
types = ["dR3", "dR4", "jet veto"]
for i in range(n_datasets):
    mid_x = channels[0] + i * (0.5*width + gap_between_sets)
    ax.text(mid_x, 1.05 * yields[i][0], types[i], ha='center', va='bottom', fontsize=12)
    # 画一个箭头指向文字
    ax.annotate('', xy=(mid_x, 1.02 * yields[i][0]), xytext=(mid_x, 1.05 * yields[i][0]), arrowprops=dict(facecolor='black', shrink=0.05))

# 设置x轴刻度
ax.set_xticks(channels + (n_datasets - 1) * (width + gap_between_sets) / 2)
ax.set_xticklabels(trees, rotation=45, ha='right')

# 添加图例
ax.legend()

# 添加标签
ax.set_xlabel('Channels')
ax.set_ylabel('Yields')

# 应用 hep 样式
hep.cms.label("Preliminary", ax=ax)

# 显示图形
plt.tight_layout()
plt.savefig("yields_change_with_cond.png")

