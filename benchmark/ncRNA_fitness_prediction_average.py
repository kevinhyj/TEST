import json
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import font_manager
import matplotlib.colors as mcolors
from matplotlib import rcParams

font_dirs = ['/home/zaitpub04/hyj/conda/conda/envs/RhoDesign2/lib/python3.9/site-packages/matplotlib/mpl-data/fonts']
font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
for font_file in font_files:
    font_manager.fontManager.addfont(font_file)
rcParams['font.family'] = 'Arial'
rcParams['font.sans-serif'] = ['Arial']
def plot_bar_chart(data_file, focus, bar_colors):
    """
    Plots a bar chart based on the selected focus ("Spearman" or "Pearson") and overlays scatter points in a beeswarm pattern on top of the bars.

    Parameters:
        data_file (str): Path to the JSON file containing the data.
        focus (str): Either "Spearman" or "Pearson".
        bar_colors (list): List of hex color codes for the bars.
    """
    # Set random seed for reproducibility
    np.random.seed(42)

    # Load the data from the JSON file
    with open(data_file, 'r') as file:
        data = json.load(file)

    # Ensure focus is valid
    if focus not in data:
        raise ValueError("Focus must be either 'Spearman' or 'Pearson'")

    # Define the x-axis categories
    x_categories = [
        "rho_score", "AIDO_score", "rinalmo_score", "rnafm_score",
        "rnamsm_score", "rnabert_score", "evo_score", "nt_score", "gena_score", "grover_score"
    ]

    # Ensure bar_colors length matches x_categories
    if len(bar_colors) != len(x_categories):
        raise ValueError("The number of bar colors must match the number of selected x-axis categories.")

    # Calculate the average values for each category using absolute values
    averages = []
    scatter_data = {}
    for category in x_categories:
        values = [
            abs(data[focus][key][category]) for key in data[focus].keys()
        ]
        averages.append(np.mean(values))
        scatter_data[category] = values

    # Plot the bar chart
    fig, ax = plt.subplots(figsize=(6, 6))
    bars = plt.bar(x_categories, averages, color=bar_colors, label="Average", zorder=1)

    # Overlay scatter points in beeswarm pattern on top of the bars
    for i, category in enumerate(x_categories):
        y_values = scatter_data[category]
        x_values = np.random.normal(i, 0.05, len(y_values))  # Add slight jitter around bar centers
        plt.scatter(x_values, y_values, color="black", s=8, alpha=1, zorder=2)

    # Add labels and title
    if focus == 'Spearman':
        plt.ylabel("|Spearman r|", fontsize=12,fontweight='semibold')
    else:
        plt.ylabel("|Pearson r|", fontsize=12,fontweight='semibold')

    plt.title("Zero-shot ncRNA Fitness Prediction", fontsize=14,fontweight='semibold')

    # Customize x-axis labels
    x_label_mapping = {
        "rho_score": "RhoDesign\n2024 | 37.1M",
        "AIDO_score": "AIDO\n2024 | 1.6B",
        "rinalmo_score": "Rinalmo\n2024 | 651M",
        "rnafm_score": "RNAFM\n2022 | 99.5M",
        "rnamsm_score": "RNAMSM\n2023 | 96.5M",
        "rnabert_score": "RNABERT\n2022 | 0.53M",
        "evo_score": "Evo\n2024 | 6.45B",
        "nt_score": "Nucl. Trans.\n2023 | 500M",
        "gena_score": "GENA\n2023 | 110M",
        "grover_score": "Grover\n2023 | 86.5M"
    }
    # Process labels to add formatting
    mapped_labels = []
    for label in x_categories:
        original_label = x_label_mapping.get(label, label)
        parts = original_label.split('\n')  # Split into model name and year/size
        if len(parts) == 2:
            # Add bold formatting to the first part and normal formatting to the second part
            formatted_label = f"$\\bf{{{parts[0]}}}$\n{parts[1]}"
            mapped_labels.append(formatted_label)
        else:
            # If the label doesn't match the expected format, just use it as is
            mapped_labels.append(label)

    # 原始的短刻度线位置
    original_positions = np.arange(len(x_categories))

    # 设置短刻度线位置
    ax.set_xticks(original_positions)  # 保持刻度线的位置
    ax.set_xticklabels([])  # 移除默认的标签

    # 自定义右移的标签位置
    shifted_positions = original_positions + 0.5  # 每个标签右移 0.2 单位

    # 绘制自定义的标签
    for x, label in zip(shifted_positions, mapped_labels):
        ax.text(
            x,  # 标签新的 x 坐标
            -0.22,  # 标签的 y 坐标（略低于 x 轴）
            label,  # 标签内容
            fontsize=9,
            ha='right',  # 标签右对齐
            rotation=45,
            fontweight='normal'  # 正常字体
        )

    # Add "better" label and arrow outside the plot (use Axes coordinates)
    ax.text(1.015, 0.1, "Better", transform=ax.transAxes,
            rotation=90, fontsize=12, color="black", ha='left', va='center')

    ax.annotate('', xy=(1.03, 1.0), xytext=(1.03, 0.2),
                xycoords=ax.transAxes,
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

    # Add gridlines in y direction
    ax.yaxis.grid(True, linestyle='--', alpha=0.7, zorder=0)
    # Add legend below the chart
    legend_colors = ['#FF6D15', '#FE9D3E', '#77A9D3']
    legend_labels = ['Structure model', 'RNA Language model', 'DNA language model']
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in legend_colors]
    ax.legend(legend_handles, legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.35), ncol=3, frameon=False,fontsize=10)
    ax.text(0.35, -0.53, "● Experimental Study", transform=ax.transAxes,
            fontsize=10, color="black", ha='left', va='center')

    # Adjust figure layout to create space for the arrow
    plt.subplots_adjust(bottom=0.45, right=0.9)
    plt.yticks(fontsize=14,fontweight='semibold')
    plt.savefig("/home/zaitpub04/hyj/RhoDesign/benchmark/paint/1/figure_1_0307.png",dpi=300)
    plt.show()

# Example usage
def main():
    plot_bar_chart(
         #"/home/zaitpub04/hyj/RhoDesign/benchmark/paint/1/data.json",
        "/home/zaitpub04/hyj/RhoDesign/benchmark/paint/1/data_6.json",
        "Spearman",
        # "Pearson",
        ["#FF6D15", "#FE9D3E", "#FFAE60", "#FFCA9B", "#FFE0C1", "#FFEDDB", "#77A9D3", "#98BBDE", "#BED6EA", "#D8E7F3"]
    )

if __name__ == "__main__":
    main()
