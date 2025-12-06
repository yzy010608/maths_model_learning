import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']                 # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False                   # 用来正常显示负号


# 设定文件路径
excel_path = "gene_expression_data.xlsx"

# 读取Excel文件
df = pd.read_excel(
    excel_path,
    index_col=0,
    header=0,
    sheet_name=0,              # 读取工作表
    na_values=["", "NA", "缺失"],  # 定义哪些值视为空值

)

# 填充空值
df = df.fillna(0)

# 相关变量设置
samples = df.index
n_genes = 5
genes = ['Gene_1', 'Gene_2', 'Gene_3', 'Gene_4', 'Gene_5']
cell_line_labels = df['Cell_Line']

# 获得中心化后的数据-->ndarray形式
df_mean = df[genes].mean()
new_data = np.zeros((len(samples), n_genes))

for x in range(1 , n_genes + 1):
    data_x = df[f'Gene_{x}']
    for i in range(len(samples)):
        new_data[i,x-1] = data_x[i] - df_mean[f'Gene_{x}']

# 计算协方差矩阵
cov = (1/149) * np.dot(new_data.T, new_data)

# 计算特征值和特征向量
eigenvalue,eigenvector = np.linalg.eigh(cov)

# 排序
sorted_idx = np.argsort(eigenvalue)[::-1]
sorted_eigenvec = eigenvector[:, sorted_idx]
sorted_eigenval = eigenvalue[sorted_idx]

# 投影，完善数据类型
proj_data = np.dot(new_data, sorted_eigenvec[:,:2])
df_proj_data = pd.DataFrame(proj_data, index=samples, columns=['PCA1','PCA2'])
df_proj_data['Cell_Line'] = cell_line_labels

# 可视化
# 标签映射为数值
cell_labels = {
    'Cell_Line_A':0,
    'Cell_Line_B':1,
    'Cell_Line_C':2
}
# 区分不同类型细胞的颜色
colors_map = {
    0: 'r',
    1: 'y',
    2: 'b'
}

plt.figure(figsize=(8,4))
plt.title("my_PCA")
for x in samples:
    plt.scatter(
        df_proj_data.loc[x,'PCA1'],
        df_proj_data.loc[x,'PCA2'],
        color = colors_map[cell_labels[df_proj_data.loc[x,'Cell_Line']]]
    )
plt.show()