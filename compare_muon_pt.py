import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 文件路径和树名
input1_file = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/11a9091d-617b-4ac9-bffd-10b526d12068.root"
tree1 = "Events"

input2_file = "/eos/project-h/htozg-dy-privatemc/HiggsDNA_skimmed/ggH_M125_2023postBPix/11a9091d-617b-4ac9-bffd-10b526d12068skimmed.root"
tree2 = "Events"

# 分析参数
DIFF_THRESHOLD = 1.0  # |diff| 阈值
DIFF_RANGE = (-1.0, 1.0)  # diff 可视化范围
BS_PT_MAX = 26.0  # 仅研究 Muon_bsConstrainedPt < 26 GeV 的μ子

def compare_muon_pt_variables(file1, tree1_name, file2, tree2_name):
    """
    比较两个ROOT文件中的Muon pt相关变量，并制作直方图
    """
    print(f"开始比较文件:")
    print(f"文件1: {file1}")
    print(f"树1: {tree1_name}")
    print(f"文件2: {file2}")
    print(f"树2: {tree2_name}")
    print("-" * 80)
    
    try:
        # 打开ROOT文件
        print("正在读取ROOT文件...")
        with uproot.open(file1) as f1, uproot.open(file2) as f2:
            # 获取树
            tree1_obj = f1[tree1_name]
            tree2_obj = f2[tree2_name]
            
            print(f"文件1中的条目数: {tree1_obj.num_entries}")
            print(f"文件2中的条目数: {tree2_obj.num_entries}")
            
            # 确保两个文件有相同的条目数
            min_entries = min(tree1_obj.num_entries, tree2_obj.num_entries)
            print(f"将比较前 {min_entries} 个条目")
            
            print("-" * 80)
            print("正在读取Muon pt相关分支...")
            
            # 读取Muon分支数据（只读取需要的条目数）
            muon_pt_file1 = tree1_obj['Muon_pt'].array(library="ak", entry_stop=min_entries)
            muon_bs_pt_file1 = tree1_obj['Muon_bsConstrainedPt'].array(library="ak", entry_stop=min_entries)
            muon_corrected_pt_file2 = tree2_obj['Muon_corrected_pt'].array(library="ak", entry_stop=min_entries)
            
            print("正在收集Muon数据...")
            
            # 收集数据
            diff1_data = []  # Muon_corrected_pt (file2) - Muon_pt (file1)
            diff2_data = []  # Muon_corrected_pt (file2) - Muon_bsConstrainedPt (file1)
            
            muon_pt_values = []
            muon_bs_pt_values = []
            muon_corrected_pt_values = []
            
            for i in range(min_entries):
                # 获取该事件的muon数据
                pt_event = muon_pt_file1[i]
                bs_pt_event = muon_bs_pt_file1[i]
                corrected_pt_event = muon_corrected_pt_file2[i]
                
                # 确保两个文件中的muon数量相同
                if len(pt_event) != len(corrected_pt_event) or len(bs_pt_event) != len(corrected_pt_event):
                    continue
                
                # 逐个muon比较
                for j in range(len(corrected_pt_event)):
                    pt = pt_event[j]
                    bs_pt = bs_pt_event[j]
                    corrected_pt = corrected_pt_event[j]
                    
                    # 仅保留 Muon_bsConstrainedPt < 26 GeV 的μ子
                    if bs_pt >= BS_PT_MAX:
                        continue
                    
                    # 计算差值
                    diff1 = corrected_pt - pt
                    diff2 = corrected_pt - bs_pt
                    
                    diff1_data.append(diff1)
                    diff2_data.append(diff2)
                    
                    muon_pt_values.append(pt)
                    muon_bs_pt_values.append(bs_pt)
                    muon_corrected_pt_values.append(corrected_pt)
            
            print(f"总共收集了 {len(diff1_data)} 个muon数据点")
            
            if len(diff1_data) == 0:
                print("错误: 没有收集到任何muon数据!")
                return False
            
            # 转换为numpy数组
            diff1_data = np.array(diff1_data)
            diff2_data = np.array(diff2_data)
            muon_pt_values = np.array(muon_pt_values)
            muon_bs_pt_values = np.array(muon_bs_pt_values)
            muon_corrected_pt_values = np.array(muon_corrected_pt_values)
            
            # 分析 ±1 附近的事例
            outliers1_mask = np.abs(diff1_data) >= DIFF_THRESHOLD
            outliers2_mask = np.abs(diff2_data) >= DIFF_THRESHOLD
            
            outliers1_count = np.sum(outliers1_mask)
            outliers2_count = np.sum(outliers2_mask)
            
            # 打印统计信息
            print("-" * 80)
            print("统计信息:")
            print(f"Muon_corrected_pt (file2) - Muon_pt (file1):")
            print(f"  总事例数: {len(diff1_data)}")
            print(f"  平均差值: {np.mean(diff1_data):.6f} GeV")
            print(f"  标准差: {np.std(diff1_data):.6f} GeV")
            print(f"  最小差值: {np.min(diff1_data):.6f} GeV")
            print(f"  最大差值: {np.max(diff1_data):.6f} GeV")
            print(f"  |差值| >= {DIFF_THRESHOLD:.0f} GeV的事例数: {outliers1_count} ({outliers1_count/len(diff1_data)*100:.2f}%)")
            if outliers1_count > 0:
                print(f"  异常值平均: {np.mean(diff1_data[outliers1_mask]):.6f} GeV")
            
            print(f"Muon_corrected_pt (file2) - Muon_bsConstrainedPt (file1):")
            print(f"  总事例数: {len(diff2_data)}")
            print(f"  平均差值: {np.mean(diff2_data):.6f} GeV")
            print(f"  标准差: {np.std(diff2_data):.6f} GeV")
            print(f"  最小差值: {np.min(diff2_data):.6f} GeV")
            print(f"  最大差值: {np.max(diff2_data):.6f} GeV")
            print(f"  |差值| >= {DIFF_THRESHOLD:.0f} GeV的事例数: {outliers2_count} ({outliers2_count/len(diff2_data)*100:.2f}%)")
            if outliers2_count > 0:
                print(f"  异常值平均: {np.mean(diff2_data[outliers2_mask]):.6f} GeV")
            
            # 制作直方图和二维分析
            print("-" * 80)
            print("正在制作直方图和二维分析...")
            
            # 设置图像参数
            try:
                plt.style.use('seaborn-v0_8')
            except:
                try:
                    plt.style.use('seaborn')
                except:
                    pass  # 使用默认样式
            
            # 创建更大的图像布局：2x3
            fig = plt.figure(figsize=(20, 12))
            
            # 创建子图
            ax1 = plt.subplot(2, 3, 1)
            ax2 = plt.subplot(2, 3, 2) 
            ax3 = plt.subplot(2, 3, 3)
            ax4 = plt.subplot(2, 3, 4)
            ax5 = plt.subplot(2, 3, 5)
            ax6 = plt.subplot(2, 3, 6)
            
            # 直方图1: Muon_corrected_pt - Muon_pt 差值分布（限制X轴范围为±1）
            n1, bins1, patches1 = ax1.hist(diff1_data, bins=100, alpha=0.7, color='blue', edgecolor='black', range=DIFF_RANGE)
            ax1.set_xlabel('Muon_corrected_pt - Muon_pt [GeV]')
            ax1.set_ylabel('Events')
            ax1.set_title(f'Difference: Muon_corrected_pt (file2) - Muon_pt (file1)\nTotal events: {len(diff1_data)}, |diff| >= {DIFF_THRESHOLD:.0f} GeV: {outliers1_count}')
            ax1.grid(True, alpha=0.3)
            ax1.axvline(np.mean(diff1_data), color='red', linestyle='--', 
                       label=f'Mean: {np.mean(diff1_data):.3f} GeV')
            ax1.axvline(DIFF_THRESHOLD, color='orange', linestyle=':', alpha=0.7, label=f'±{DIFF_THRESHOLD:.0f} GeV')
            ax1.axvline(-DIFF_THRESHOLD, color='orange', linestyle=':', alpha=0.7)
            ax1.legend()
            ax1.set_xlim(DIFF_RANGE)
            
            # 在直方图上标注事例数
            visible_mask1 = (diff1_data >= DIFF_RANGE[0]) & (diff1_data <= DIFF_RANGE[1])
            visible_count1 = np.sum(visible_mask1)
            ax1.text(0.02, 0.98, f'Visible: {visible_count1}/{len(diff1_data)}', 
                    transform=ax1.transAxes, verticalalignment='top', 
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # 直方图2: Muon_corrected_pt - Muon_bsConstrainedPt 差值分布（限制X轴范围为±1）
            n2, bins2, patches2 = ax2.hist(diff2_data, bins=100, alpha=0.7, color='green', edgecolor='black', range=DIFF_RANGE)
            ax2.set_xlabel('Muon_corrected_pt - Muon_bsConstrainedPt [GeV]')
            ax2.set_ylabel('Events')
            ax2.set_title(f'Difference: Muon_corrected_pt (file2) - Muon_bsConstrainedPt (file1)\nTotal events: {len(diff2_data)}, |diff| >= {DIFF_THRESHOLD:.0f} GeV: {outliers2_count}')
            ax2.grid(True, alpha=0.3)
            ax2.axvline(np.mean(diff2_data), color='red', linestyle='--', 
                       label=f'Mean: {np.mean(diff2_data):.3f} GeV')
            ax2.axvline(DIFF_THRESHOLD, color='orange', linestyle=':', alpha=0.7, label=f'±{DIFF_THRESHOLD:.0f} GeV')
            ax2.axvline(-DIFF_THRESHOLD, color='orange', linestyle=':', alpha=0.7)
            ax2.legend()
            ax2.set_xlim(DIFF_RANGE)
            
            # 在直方图上标注事例数
            visible_mask2 = (diff2_data >= DIFF_RANGE[0]) & (diff2_data <= DIFF_RANGE[1])
            visible_count2 = np.sum(visible_mask2)
            ax2.text(0.02, 0.98, f'Visible: {visible_count2}/{len(diff2_data)}', 
                    transform=ax2.transAxes, verticalalignment='top', 
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # 散点图3: corrected_pt vs pt
            ax3.scatter(muon_pt_values, muon_corrected_pt_values, alpha=0.3, s=1)
            ax3.plot([np.min(muon_pt_values), np.max(muon_pt_values)], 
                    [np.min(muon_pt_values), np.max(muon_pt_values)], 
                    'r--', label='y=x')
            ax3.set_xlabel('Muon_pt (file1) [GeV]')
            ax3.set_ylabel('Muon_corrected_pt (file2) [GeV]')
            ax3.set_title(f'Muon_corrected_pt vs Muon_pt\n{len(muon_pt_values)} points')
            ax3.grid(True, alpha=0.3)
            ax3.legend()
            
            # 散点图4: corrected_pt vs bsConstrainedPt
            ax4.scatter(muon_bs_pt_values, muon_corrected_pt_values, alpha=0.3, s=1, color='green')
            ax4.plot([np.min(muon_bs_pt_values), np.max(muon_bs_pt_values)], 
                    [np.min(muon_bs_pt_values), np.max(muon_bs_pt_values)], 
                    'r--', label='y=x')
            ax4.set_xlabel('Muon_bsConstrainedPt (file1) [GeV]')
            ax4.set_ylabel('Muon_corrected_pt (file2) [GeV]')
            ax4.set_title(f'Muon_corrected_pt vs Muon_bsConstrainedPt\n{len(muon_bs_pt_values)} points')
            ax4.grid(True, alpha=0.3)
            ax4.legend()
            
            # 二维分布图5: diff1 vs muon_bsConstrainedPt 的热图
            print("制作二维分布分析...")
            # 创建二维直方图
            hist2d, xedges, yedges = np.histogram2d(muon_bs_pt_values, diff1_data, bins=[50, 100], 
                                                   range=[[0, BS_PT_MAX], [DIFF_RANGE[0], DIFF_RANGE[1]]])
            X, Y = np.meshgrid(xedges, yedges)
            im1 = ax5.pcolormesh(X, Y, hist2d.T, cmap='viridis', shading='auto')
            ax5.set_xlabel('Muon_bsConstrainedPt [GeV]')
            ax5.set_ylabel('corrected_pt - pt [GeV]')
            ax5.set_title('2D Distribution: (corrected_pt - pt) vs bsConstrainedPt')
            ax5.axhline(0, color='red', linestyle='-', alpha=0.7, linewidth=1)
            ax5.axhline(DIFF_THRESHOLD, color='orange', linestyle='--', alpha=0.7)
            ax5.axhline(-DIFF_THRESHOLD, color='orange', linestyle='--', alpha=0.7)
            plt.colorbar(im1, ax=ax5, label='Count')
            ax5.grid(True, alpha=0.3)
            
            # 二维分布图6: diff2 vs muon_bsConstrainedPt 的热图
            hist2d2, xedges2, yedges2 = np.histogram2d(muon_bs_pt_values, diff2_data, bins=[50, 100], 
                                                      range=[[0, BS_PT_MAX], [DIFF_RANGE[0], DIFF_RANGE[1]]])
            X2, Y2 = np.meshgrid(xedges2, yedges2)
            im2 = ax6.pcolormesh(X2, Y2, hist2d2.T, cmap='plasma', shading='auto')
            ax6.set_xlabel('Muon_bsConstrainedPt [GeV]')
            ax6.set_ylabel('corrected_pt - bsConstrainedPt [GeV]')
            ax6.set_title('2D Distribution: (corrected_pt - bsConstrainedPt) vs bsConstrainedPt')
            ax6.axhline(0, color='red', linestyle='-', alpha=0.7, linewidth=1)
            ax6.axhline(DIFF_THRESHOLD, color='orange', linestyle='--', alpha=0.7)
            ax6.axhline(-DIFF_THRESHOLD, color='orange', linestyle='--', alpha=0.7)
            plt.colorbar(im2, ax=ax6, label='Count')
            ax6.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            # 保存图像
            output_filename = 'muon_pt_comparison.png'
            plt.savefig(output_filename, dpi=300, bbox_inches='tight')
            print(f"图像已保存为: {output_filename}")
            
            # 额外分析：详细查看±100附近的事例
            print("-" * 80)
            print("±1 GeV附近异常值详细分析:")
            
            # 分析不同bsConstrainedPt区间的差值分布
            pt_bins = [0, 10, 15, 20, 26]  # 仅在 <26 GeV 区间内做分段
            print(f"\n不同bsConstrainedPt区间的差值统计:")
            print("bsConstrainedPt区间     事例数    diff1平均    diff1标准差    |diff1|>=1数      diff2平均    diff2标准差    |diff2|>=1数")
            print("-" * 120)
            
            for i in range(len(pt_bins)-1):
                pt_min, pt_max = pt_bins[i], pt_bins[i+1]
                mask = (muon_bs_pt_values >= pt_min) & (muon_bs_pt_values < pt_max)
                
                if np.sum(mask) > 0:
                    diff1_bin = diff1_data[mask]
                    diff2_bin = diff2_data[mask]
                    
                    outliers1_bin = np.sum(np.abs(diff1_bin) >= DIFF_THRESHOLD)
                    outliers2_bin = np.sum(np.abs(diff2_bin) >= DIFF_THRESHOLD)
                    
                    print(f"[{pt_min:3.0f}, {pt_max:3.0f}) GeV    {np.sum(mask):8d}  {np.mean(diff1_bin):10.3f}  {np.std(diff1_bin):11.3f}  {outliers1_bin:12d}  {np.mean(diff2_bin):10.3f}  {np.std(diff2_bin):11.3f}  {outliers2_bin:12d}")
            
            if outliers1_count > 0:
                outlier_values1 = diff1_data[outliers1_mask]
                outlier_bs_pt1 = muon_bs_pt_values[outliers1_mask]
                print(f"\ncorrected_pt - pt 异常值 (|diff| >= {DIFF_THRESHOLD:.0f}):")
                print(f"  数量: {len(outlier_values1)}")
                print(f"  差值范围: [{np.min(outlier_values1):.2f}, {np.max(outlier_values1):.2f}] GeV")
                print(f"  差值平均值: {np.mean(outlier_values1):.3f} GeV")
                print(f"  差值中位数: {np.median(outlier_values1):.3f} GeV")
                print(f"  对应bsConstrainedPt范围: [{np.min(outlier_bs_pt1):.2f}, {np.max(outlier_bs_pt1):.2f}] GeV")
                print(f"  对应bsConstrainedPt平均值: {np.mean(outlier_bs_pt1):.3f} GeV")
                
            if outliers2_count > 0:
                outlier_values2 = diff2_data[outliers2_mask]
                outlier_bs_pt2 = muon_bs_pt_values[outliers2_mask]
                print(f"\ncorrected_pt - bsConstrainedPt 异常值 (|diff| >= {DIFF_THRESHOLD:.0f}):")
                print(f"  数量: {len(outlier_values2)}")
                print(f"  差值范围: [{np.min(outlier_values2):.2f}, {np.max(outlier_values2):.2f}] GeV")
                print(f"  差值平均值: {np.mean(outlier_values2):.3f} GeV")
                print(f"  差值中位数: {np.median(outlier_values2):.3f} GeV")
                print(f"  对应bsConstrainedPt范围: [{np.min(outlier_bs_pt2):.2f}, {np.max(outlier_bs_pt2):.2f}] GeV")
                print(f"  对应bsConstrainedPt平均值: {np.mean(outlier_bs_pt2):.3f} GeV")
            
            plt.show()
            
            return True
                
    except Exception as e:
        print(f"比较过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    # 执行比较
    result = compare_muon_pt_variables(input1_file, tree1, input2_file, tree2)
    
    if result:
        print("\n🎉 比较完成!")
    else:
        print("\n⚠️  比较过程中出现问题!")
