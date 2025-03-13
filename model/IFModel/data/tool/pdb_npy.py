import os
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

def pdb_to_secondary_structure(pdb_file, output_file):
    """
    将 PDB 文件解析为二级结构并保存为 .npy 格式文件。
    
    :param pdb_file: 输入的 PDB 文件路径
    :param output_file: 输出的 .npy 文件路径
    """
    # 初始化 PDB 解析器
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # 获取第一个模型（通常 PDB 文件中只有一个模型）
    model = structure[0]

    # 使用 DSSP 计算二级结构
    dssp = DSSP(model, pdb_file)

    # 提取二级结构信息
    residues = [key[1] for key in dssp.keys()]  # 提取残基编号
    secondary_structure = [dssp[key][2] for key in dssp.keys()]  # 提取二级结构类型

    # 将二级结构类型转换为数字（便于保存为数组）
    ss_mapping = {
        "H": 0,  # α-螺旋
        "B": 1,  # β-桥
        "E": 2,  # β-链
        "G": 3,  # 310 螺旋
        "I": 4,  # π 螺旋
        "T": 5,  # 转角
        "S": 6,  # 弯曲
        "-": 7   # 无定义
    }
    secondary_structure_numeric = [ss_mapping.get(ss, 7) for ss in secondary_structure]

    # 保存为 .npy 格式
    np.save(output_file, np.array(secondary_structure_numeric))
    print(f"二级结构数据已保存到 {output_file}")

# 示例用法
pdb_file = "/home/zaitpub04/hyj/RhoDesign/RhoDesign/data/test/af_bro.pdb"  # 输入 PDB 文件路径
output_file = "/home/zaitpub04/hyj/RhoDesign/RhoDesign/data/test_ss"  # 输出 .npy 文件路径

if os.path.exists(pdb_file):
    pdb_to_secondary_structure(pdb_file, output_file)
else:
    print(f"文件 {pdb_file} 不存在，请检查路径。")
