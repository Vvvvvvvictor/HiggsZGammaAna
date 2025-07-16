import numpy as np
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
from pdb import set_trace

def plot_weight_distribution(file_name, weight_variable="weight_pu_reweight_sf_central", output_prefix="PU_weight"):
    """
    绘制指定权重变量的分布
    
    参数:
    file_name: ROOT文件路径
    weight_variable: 要绘制的权重变量名
    output_prefix: 输出图片的前缀名
    """
    input_file = uproot.open(file_name)
    
    for tree in ["zero_jet", "one_jet", "two_jet", "ZH", "ttH_had", "ttH_lep"]:
        try:
            # 读取权重数据
            weight_data = input_file[tree].arrays([weight_variable], library="pd")[weight_variable]
            
            # 获取唯一值和对应的事件数
            weight_point = np.unique(weight_data)
            weight_dis = np.array([np.sum(weight_data[weight_data==point]) for point in weight_point])
            
            # 绘制散点图
            plt.scatter(weight_point, weight_dis, label=tree, s=3)
            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel(f"{weight_variable}")
            plt.ylabel("Number of events")
            plt.legend()
            
            # 保存图片
            output_path = f"/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/{output_prefix}_{file_name.split('/')[-1].replace('.root', f'_{tree}.png')}"
            plt.savefig(output_path)
            plt.clf()
            
            print(f"Saved plot for {tree}: {output_path}")
            
        except Exception as e:
            print(f"Error processing tree {tree}: {e}")
            continue

# 使用示例
if __name__ == "__main__":
    # 绘制 weight_pu_reweight_sf_central 分布
    file_path = "/eos/home-j/jiehan/root/skimmed_ntuples_run3/WminusH_M125/2022postEE.root"
    plot_weight_distribution(file_path, "weight_pu_reweight_sf_central", "PU_weight")
    
    # 也可以绘制其他权重变量，例如：
    # plot_weight_distribution(file_path, "weight", "total_weight")
