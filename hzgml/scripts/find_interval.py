import os
import uproot
import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace

# 指定文件夹路径
folder_path = "/eos/user/j/jiehan/root/skimmed_ntuples/"

# 定义一个函数来获取文件夹中所有的.root文件路径
def get_root_files(folder_path, folders):
    root_files = []
    for folder in folders:
        folder_full_path = os.path.join(folder_path, folder)
        if os.path.exists(folder_full_path):
            for file_name in os.listdir(folder_full_path):
                if file_name.endswith(".root"):
                    root_files.append(os.path.join(folder_full_path, file_name))
    return root_files

# 定义一个函数来找到包含指定百分比的区间
def find_interval(hist, target_percent):
    # total_entries = np.sum(hist)
    # # 找到大于等于目标事例数的最小区间
    # cum_entries = np.cumsum(hist)
    # target_entries = (1-target_percent)/2 * total_entries
    # left_index = np.argmax(cum_entries >= target_entries)
    # target_entries = (1+target_percent)/2 * total_entries
    # right_index = np.argmax(cum_entries >= target_entries)

    peak = max(hist)
    peak_index = np.argmax(hist)
    left_index = np.argmax(hist[:peak_index] >= peak*0.5)
    right_index = peak_index + np.argmax(hist[peak_index:] <= peak*0.5)
    return left_index, right_index

# 指定需要读取的文件夹列表
folders_to_read = ["ggH", "VBF", "WplusH", "WminusH", "ZH", "ttH"]

# 获取文件夹中所有的.root文件路径
root_files = get_root_files(folder_path, folders_to_read)

# 创建一个空的直方图
bins = np.linspace(100, 150, 200)  # 设置直方图的区间和bin数目
total_hist = np.zeros(len(bins) - 1)

# 读取所有的.root文件并累积数据到总的直方图中
for file_path in root_files:
    print(f"Reading file: {file_path}")
    file = uproot.open(file_path)
    tree = file["inclusive"]
    h_mass_values = tree["H_mass"].array()
    weights = tree["weight"].array()
    
    # 更新总的直方图
    hist, _ = np.histogram(h_mass_values, bins=bins, weights=weights)
    total_hist += hist

# 找到总的直方图的峰值
peak_index = np.argmax(total_hist)
peak_value = bins[peak_index]

# 找到包含总事例左右共target_percent的区间
target_percent = 0.68269  # target_percent对应的小数值
# left_index, _ = find_interval(total_hist[:peak_index], target_percent)
# _, right_index = find_interval(total_hist[peak_index+1:], target_percent)
# right_index += peak_index
left_index, right_index = find_interval(total_hist, target_percent)
left_point = bins[left_index]
right_point = bins[right_index]

# 绘制总的直方图并标记峰值两侧的点
plt.hist(bins[:-1], bins=bins, weights=total_hist, alpha=0.7, label="H_mass")
plt.axvline(x=peak_value, color='r', linestyle='--', label="Peak")
plt.axvline(x=left_point, color='g', linestyle='--', label="Left Point")
plt.axvline(x=right_point, color='g', linestyle='--', label="Right Point")
plt.xlabel(r"H_mass(GeV/$c^2$)")
plt.ylabel(r"Events(0.25GeV/$c^2$)")
plt.title('90% Interval of Signal Samples')
plt.legend()

height = max(total_hist) / 2
plt.text(peak_value, height*1.5, f'Peak Value: {peak_value:.2f}', ha='center', va='bottom')
plt.text(left_point, height, f'Left Point: {left_point:.2f}', ha='right', va='bottom')
plt.text(right_point, height, f'Right Point: {right_point:.2f}', ha='left', va='bottom')

plt.savefig("figures/three_sigma_interval.png")

# 输出峰值两侧的点的位置
print(f"Peak Value: {peak_value}")
print(f"Left Point: {left_point}, ratio to peak: {total_hist[left_index]/height/2}")
print(f"Right Point: {right_point}, ratio to peak: {total_hist[right_index]/height/2}")
