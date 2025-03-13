import numpy as np

# 定义两个.npy文件的路径
file_path_1 = '/home/zaitpub04/hyj/RhoDesign/RhoDesign/data/npy/pepper_af.npy'
file_path_2 = '/home/zaitpub04/hyj/RhoDesign/RhoDesign/data/npy/rho_pepper.npy'

# 加载.npy文件
array_1 = np.load(file_path_1)
array_2 = np.load(file_path_2)

# 比较两个数组是否相同
are_equal = np.array_equal(array_1, array_2)

# 打印结果
if are_equal:
    print("两个.npy文件完全相同。")
else:
    print("两个.npy文件不相同。")