import uproot
import numpy as np
import pandas as pd
from collections import defaultdict
import sys

# 文件路径和树名
input1_file = "/eos/project-h/htozg-dy-privatemc/rzou/bdt/BDT_input_new/SM1_2018_pinnacles_all_fixedjet.root"
tree1 = "TreeB"

input2_file = "/eos/project-h/htozg-dy-privatemc/rzou/bdt/Output_ggF_rui_commonparam/SM1_2018_output.root"
tree2 = "outtree"

def compare_root_trees(file1, tree1_name, file2, tree2_name):
    """
    基于event分支比较两个ROOT文件中的树，检查相同事件ID的所有分支数值是否一致
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
            if tree1_name not in f1:
                print(f"错误: 文件1中找不到树 '{tree1_name}'")
                print(f"可用的键: {list(f1.keys())}")
                return False
                
            if tree2_name not in f2:
                print(f"错误: 文件2中找不到树 '{tree2_name}'")
                print(f"可用的键: {list(f2.keys())}")
                return False
                
            tree1_obj = f1[tree1_name]
            tree2_obj = f2[tree2_name]
            
            print(f"文件1中的条目数: {tree1_obj.num_entries}")
            print(f"文件2中的条目数: {tree2_obj.num_entries}")
            
            # 获取分支名称
            branches1 = set(tree1_obj.keys())
            branches2 = set(tree2_obj.keys())
            
            print(f"文件1中的分支数: {len(branches1)}")
            print(f"文件2中的分支数: {len(branches2)}")
            
            # 检查是否存在event分支
            if 'event' not in branches1:
                print("错误: 文件1中找不到'event'分支!")
                print(f"可用分支: {sorted(branches1)}")
                return False
                
            if 'event' not in branches2:
                print("错误: 文件2中找不到'event'分支!")
                print(f"可用分支: {sorted(branches2)}")
                return False
            
            # 找到共同的分支
            common_branches = branches1.intersection(branches2)
            only_in_1 = branches1 - branches2
            only_in_2 = branches2 - branches1
            
            print(f"共同分支数: {len(common_branches)}")
            
            if only_in_1:
                print(f"只在文件1中存在的分支 ({len(only_in_1)}): {sorted(only_in_1)}")
            if only_in_2:
                print(f"只在文件2中存在的分支 ({len(only_in_2)}): {sorted(only_in_2)}")
            
            # 移除event分支，因为它用作索引
            comparison_branches = common_branches - {'event'}
            
            if not comparison_branches:
                print("错误: 除了event分支外，没有找到其他共同的分支!")
                return False
            
            print(f"将比较的分支数 (除event外): {len(comparison_branches)}")
            print("-" * 80)
            
            # 读取event分支作为事件标识符
            print("正在读取event分支...")
            events1 = tree1_obj['event'].array(library="np")
            events2 = tree2_obj['event'].array(library="np")
            
            print(f"文件1中的唯一事件数: {len(np.unique(events1))}")
            print(f"文件2中的唯一事件数: {len(np.unique(events2))}")
            
            # 创建事件到索引的映射
            event_to_idx1 = {event: idx for idx, event in enumerate(events1)}
            event_to_idx2 = {event: idx for idx, event in enumerate(events2)}
            
            # 找到共同的事件
            common_events = set(events1).intersection(set(events2))
            only_in_file1 = set(events1) - set(events2)
            only_in_file2 = set(events2) - set(events1)
            
            print(f"共同事件数: {len(common_events)}")
            if only_in_file1:
                print(f"只在文件1中的事件数: {len(only_in_file1)}")
                if len(only_in_file1) <= 10:
                    print(f"  事件ID: {sorted(only_in_file1)}")
            if only_in_file2:
                print(f"只在文件2中的事件数: {len(only_in_file2)}")
                if len(only_in_file2) <= 10:
                    print(f"  事件ID: {sorted(only_in_file2)}")
            
            if not common_events:
                print("错误: 没有找到共同的事件!")
                return False
            
            print("-" * 80)
            print("开始逐分支比较共同事件的数据...")
            
            # 读取所有需要比较的分支数据
            branch_data1 = {}
            branch_data2 = {}
            
            print("正在读取分支数据...")
            for branch in comparison_branches:
                try:
                    branch_data1[branch] = tree1_obj[branch].array(library="np")
                    branch_data2[branch] = tree2_obj[branch].array(library="np")
                except Exception as e:
                    print(f"  ⚠️  读取分支 {branch} 时出错: {e}")
                    continue
            
            # 比较每个分支的共同事件
            differences_found = False
            total_mismatches = 0
            
            for branch in sorted(branch_data1.keys()):
                print(f"正在比较分支: {branch}")
                
                data1 = branch_data1[branch]
                data2 = branch_data2[branch]
                
                branch_mismatches = 0
                mismatch_details = []
                
                for event_id in sorted(common_events):
                    idx1 = event_to_idx1[event_id]
                    idx2 = event_to_idx2[event_id]
                    
                    val1 = data1[idx1]
                    val2 = data2[idx2]
                    
                    # 比较数值
                    is_equal = False
                    if isinstance(val1, (int, np.integer)) and isinstance(val2, (int, np.integer)):
                        is_equal = val1 == val2
                    elif isinstance(val1, (float, np.floating)) and isinstance(val2, (float, np.floating)):
                        # 对浮点数使用容差比较
                        is_equal = np.isclose(val1, val2, rtol=1e-9, atol=1e-12, equal_nan=True)
                    else:
                        # 其他类型直接比较
                        is_equal = val1 == val2
                    
                    if not is_equal:
                        branch_mismatches += 1
                        if len(mismatch_details) < 5:  # 只记录前5个不匹配的详情
                            mismatch_details.append((event_id, val1, val2))
                
                if branch_mismatches > 0:
                    print(f"  ❌ 发现 {branch_mismatches} 个事件的数值不匹配")
                    for event_id, val1, val2 in mismatch_details:
                        print(f"    事件 {event_id}: {val1} vs {val2}")
                    if branch_mismatches > 5:
                        print(f"    ... 还有 {branch_mismatches - 5} 个不匹配的事件")
                    differences_found = True
                    total_mismatches += branch_mismatches
                else:
                    print(f"  ✅ 所有共同事件的数值都匹配")
            
            print("-" * 80)
            print(f"比较完成!")
            print(f"共同事件数: {len(common_events)}")
            print(f"比较的分支数: {len(branch_data1)}")
            print(f"总不匹配数: {total_mismatches}")
            
            if differences_found:
                print("❌ 发现差异: 存在相同事件ID但数值不同的情况")
                return False
            else:
                print("✅ 所有共同事件在所有共同分支中的数值都完全匹配!")
                return True
                
    except Exception as e:
        print(f"比较过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    # 执行比较
    result = compare_root_trees(input1_file, tree1, input2_file, tree2)
    
    if result:
        print("\n🎉 结论: 两个ROOT文件中的事件数据完全相同!")
    else:
        print("\n⚠️  结论: 两个ROOT文件中的事件数据存在差异!")
        sys.exit(1)