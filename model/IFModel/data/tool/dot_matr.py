import numpy as np
import os

def dbn_to_adjacency_matrix(dbn_text):
    """
    将dot-bracket notation (dbn)格式的文本转换为相应的邻接矩阵,并存为.npy文件.
    
    参数:
    dbn_text (str): 包含dot-bracket notation格式的RNA二级结构的文本.
    
    返回:
    numpy.ndarray: 表示RNA二级结构的邻接矩阵.
    """
    # 去除文本中的换行符和空格
    dbn_text = dbn_text.replace('\n', '').replace(' ', '')
    
    # 创建一个全零的邻接矩阵
    n = len(dbn_text)
    adjacency_matrix = np.zeros((n, n), dtype=int)
    
    # 填充邻接矩阵
    stack = []
    for i, char in enumerate(dbn_text):
        if char == '(':
            stack.append(i)
        elif char == ')':
            j = stack.pop()
            adjacency_matrix[i, j] = 1
            adjacency_matrix[j, i] = 1
    
    # 保存邻接矩阵为.npy文件
    current_dir = os.path.dirname(os.path.abspath(__file__))
    output_file = os.path.join(current_dir, 'mango3.npy')
    np.save(output_file, adjacency_matrix)
    
    return adjacency_matrix

# 示例用法
dbn_text = '((((((((......((......))......))))))))'
adjacency_matrix = dbn_to_adjacency_matrix(dbn_text)
print(adjacency_matrix)