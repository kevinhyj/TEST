import numpy as np
import os
import torch

# 定义文件路径
input_dir = '/home/zaitpub04/hyj/RhoDesign/RhoDesign/data'
output_dir = '/home/zaitpub04/hyj/RhoDesign/RhoDesign/data/test_ss'

# 确保输出目录存在
os.makedirs(output_dir, exist_ok=True)

# 遍历输入目录中的所有 .txt 文件
for filename in os.listdir(input_dir):
    if filename.endswith('.txt'):
        # 构建完整的文件路径
        txt_file_path = os.path.join(input_dir, filename)
        npy_file_path = os.path.join(output_dir, filename.replace('.txt', '.npy'))
        
        # 读取 .txt 文件内容
        with open(txt_file_path, 'r') as file:
            content = file.read().strip()
        
        # 将内容转换为布尔类型的数组
        data = np.array([char != '.' for char in content], dtype=bool)
        
        # 保存为 .npy 文件
        np.save(npy_file_path, data)

print("转换完成！")