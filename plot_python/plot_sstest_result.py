import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
import os
plt.style.use(hep.style.CMS)

# limits = [105, 45, 12, 5]
path = "/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/SSTest/outputs/"
for dir_name in os.listdir(path):
    dir_path = path+dir_name
    if os.path.isdir(dir_path):
        for file_name in os.listdir(dir_path):
            file_path = dir_path + "/" + file_name
            # load data
            categories, measurements, stat_errors, data_errors, total_errors = [], [], [], [], []
            with open(file_path, 'r') as f:
                for line in f:
                    if 'status = Pass' in line:
                        data = line.split('\t')
                        category = data[1]
                        measurement = float(data[3].split(":")[1])
                        stat_error = float(data[4].split("=")[1])
                        data_error = float(data[6].split("=")[1])
                        total_error = float(data[7].split("=")[1])
                        categories.append(category)
                        measurements.append(measurement)
                        stat_errors.append(stat_error)
                        data_errors.append(data_error)
                        total_errors.append(total_error)

            # 定义y位置
            y_pos = np.arange(len(categories))

            # 创建图
            fig, ax = plt.subplots()

            print(measurements, y_pos, stat_errors, data_errors, total_errors)

            for i in range(len(categories)):
                # 绘制统计误差的阴影区域（橙色）
                ax.fill_betweenx([y_pos[i]-0.1, y_pos[i]+0.1], -data_errors[i], 
                                data_errors[i], 
                                color='orange', alpha=0.3, label=(r'$\pm\Delta_{SS}$' if i == 0 else None))
                
                # 绘制系统误差的阴影区域（蓝色）
                ax.fill_betweenx([y_pos[i]-0.1, y_pos[i]+0.1], -0.2*data_errors[i], 
                                0.2*data_errors[i], 
                                color='blue', alpha=0.3, label=(r'$\pm0.2\times\Delta_{SS}$' if i == 0 else None))
            
            
            # 绘制系统误差条（用蓝色线条）
            ax.errorbar(measurements, y_pos, xerr=2*np.array(stat_errors), fmt='o', color='black', 
                        ecolor='green', capsize=5, elinewidth=2, label=r'SS$\pm2\Delta_{MC}$', linestyle='none', marker='s', markersize=7)

            # 绘制统计误差条（用蓝色方块）
            ax.errorbar(measurements, y_pos, xerr=stat_errors, fmt='o', color='black', 
                        ecolor='red', capsize=5, elinewidth=2, label=r'SS$\pm1\Delta_{MC}$', linestyle='none', marker='s', markersize=7)

            ax.plot(total_errors, y_pos, 'v--', color='black', label=r'$\sqrt{SS^2+\Delta_{SS}^2}$', markersize=7)

            # 设置分类标签
            ax.set_yticks(y_pos)
            ax.set_yticklabels(categories)

            # # 添加水平线
            # ax.axvline(x=0, color='red', linestyle='--', label='SM')

            # 添加标题和标签
            ax.set_xlabel(r'a.u.')
            ax.set_title(f'Spurious Signal Test Result(cat{file_name.split("_")[-2]})')

            # 显示图例 and add a shadow
            ax.legend(ncol=3, fontsize=20)
            ax.set_ylim(-0.5, len(categories)/0.875-0.5)
            
            cat = int(file_name.split("_")[-2])
            # ax.set_xlim(-limits[cat], limits[cat])

            # 显示图表
            plt.tight_layout()
            plt.savefig(f'/afs/cern.ch/user/j/jiehan/private/HiggsZGammaAna/plot_python/pic/sstest/{dir_name}_{cat}.png')
