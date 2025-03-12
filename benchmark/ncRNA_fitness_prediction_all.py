import matplotlib.pyplot as plt
import numpy as np
import json
from matplotlib import rcParams
from matplotlib import font_manager
import matplotlib.colors as mcolors
# custom_font_manager = font_manager.FontManager()

# print(font_manager.findSystemFonts(fontpaths=None, fontext='ttf')[10:20])

font_dirs = ['/home/zaitpub04/hyj/conda/conda/envs/RhoDesign2/lib/python3.9/site-packages/matplotlib/mpl-data/fonts']
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
for font_file in font_files:
    font_manager.fontManager.addfont(font_file)
rcParams['font.family'] = 'Arial'
rcParams['font.sans-serif'] = ['Arial']
# 从文件读取 JSON 数据
file_path = "/home/zaitpub04/hyj/RhoDesign/benchmark/paint/1/data copy.json"  # 确保文件路径正确
with open(file_path, "r") as file:
    data = json.load(file)

# 参数选择：选择 "Spearman" 或 "Pearson"
selection = "Spearman"  # 可以更改为 "Pearson"

# 选定数据
selected_data = data[selection]

# 定义颜色
colors = ["#FF6D15", "#FE9D3E", "#FFAE60", "#FFCA9B", "#FFE0C1", "#FFEDDB", "#77A9D3", "#98BBDE", "#BED6EA", "#D8E7F3"]

# 提取分组和键值


groups = list(selected_data.keys())
keys = list(selected_data[groups[0]].keys())
values_per_group = [[abs(selected_data[group][key]) for key in keys] for group in groups]

# 调整柱状图布局
plt.figure(figsize=(20, 5))  # 图形宽度调整为 12
bar_width = 0.85 / len(keys)  # 调整柱子宽度
x_positions = np.linspace(0, len(groups) - 1, len(groups))  # 缩小左右两端间距

# 绘制柱状图
for i, key in enumerate(keys):
    group_values = [values[i] for values in values_per_group]
    plt.bar(
        x_positions + i * bar_width - bar_width * len(keys) / 2,
        group_values,
        width=bar_width,
        color=colors[i],
        edgecolor='black',  # 添加柱子外边框
        label=key
    )

# 将键名设置为 x 轴标签
key_name=['RhoDesign','Rinalmo','AIDO','Rnafm','Rnabert','Rnamsm',"Grover","GENA", "NT-Transformer","Evo"]
key_labels = [key for _ in groups for key in key_name]
xtick_positions = [x + i * bar_width - bar_width * len(keys) / 2 for x in x_positions for i in range(len(keys))]
plt.xticks(xtick_positions, key_labels, rotation=45, fontsize=12, ha='right',fontweight='bold')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 设置分组名称
print('groups:',groups)
for idx, group in enumerate(groups):
    plt.text(
        x_positions[idx], 0.63, group, fontsize=16, ha='center',fontweight='bold'
    )

# 轴标签和标题
if selection == 'Spearman':
    plt.ylabel("|Spearman r|", fontsize=18,fontweight='semibold')
else:
    plt.ylabel("|Pearson r|", fontsize=18,fontweight='semibold')
plt.title(f"Zero-shot Performance on Different Dataset",fontsize=20,fontweight='semibold')

# 去掉图例
plt.legend().set_visible(False)

# 调整 x 轴显示范围，缩小两端空白
plt.xlim(-0.58, len(groups) - 0.5)  # 左右两端稍微裁剪掉

plt.yticks(fontsize=16,fontweight='bold') 
# 调整布局
plt.subplots_adjust(bottom=0.4)

# 显示图形
plt.show()
plt.savefig("/home/zaitpub04/hyj/RhoDesign/benchmark/paint/figure_2.png",dpi=300)