import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']                 # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False                   # 用来正常显示负号

# 设定文件路径
excel_path = "data_1.xlsx"

# 读取Excel文件
df = pd.read_excel(
    excel_path,
    sheet_name=0,              # 读取工作表
    na_values=["", "NA", "缺失"],  # 定义哪些值视为空值
    dtype={"X1": float, "X2": float}  # 指定列的数据类型
)

# 填充空值（避免转换后出现nan）
df = df.fillna(0)

# 仅保留数值类型的列（过滤文本列）
new_df = df.select_dtypes(include=[np.number])

# 转换为ndarray
dataset = new_df.to_numpy()

# 数据预处理
# 取到变量1和2对应的数据
X1 = dataset[:,0]
X2 = dataset[:,1]
X1_mean = np.mean(X1)
X2_mean = np.mean(X2)
data_mean = np.array([X1_mean, X2_mean])
X_mean = np.tile(data_mean,(100,1))

# 中心化处理
data = dataset - X_mean

# 计算协方差矩阵-->实对称矩阵
cov = (1/99) * np.dot(data.T, data)

# 计算特征值和特征向量
eigenvalue,eigenvector = np.linalg.eigh(cov)  # eigenvalue 从小到大排列

# 计算变量（中心化）在主成分中的投影
new_data = np.dot(data, eigenvector)

# 可视化
# 第一幅图：原始数据与中心化数据对比
plt.figure(figsize=(8,4))
# 绘制原始数据
plt.scatter(dataset[:, 0], dataset[:, 1], c='blue', alpha=0.6, label='原始数据')
# 绘制中心化数据
plt.scatter(data[:, 0], data[:, 1], c='orange', alpha=0.6, label='中心化数据')
# 标记均值点
plt.scatter(data_mean[0], data_mean[1], c='red', s=100, marker='X', label='均值点')
# 绘制坐标轴参考线
plt.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
plt.xlabel('X1')
plt.ylabel('X2')
plt.title('原始数据与中心化数据')
plt.legend()
plt.grid(alpha=0.3)
plt.show()

# 第二幅图：主成分方向与投影结果
plt.figure(figsize=(8,8))
# 绘制中心化数据
plt.scatter(data[:, 0], data[:, 1], c='orange', alpha=0.6, label='中心化数据')

# 绘制主成分方向（特征向量），按特征值大小排序后绘制
# 注意：eigh返回的特征值是从小到大排列的，这里将其反转
sorted_idx = np.argsort(eigenvalue)[::-1]
sorted_eigenvec = eigenvector[:, sorted_idx]
sorted_eigenval = eigenvalue[sorted_idx]

# 主成分方向线（以均值为起点，按特征值大小缩放显示）
for i in range(2):
    vec = sorted_eigenvec[:, i]
    val = sorted_eigenval[i]
    # 按特征值开方缩放向量（体现方差大小）
    plt.quiver(0, 0, vec[0]*np.sqrt(val), vec[1]*np.sqrt(val),
               color=['red', 'green'][i], scale=1, scale_units='xy',
               label=f'主成分{i+1} (λ={val:.2f})')

# 绘制数据在第一主成分上的投影
proj_on_pc1 = np.outer(new_data[:, sorted_idx[0]], sorted_eigenvec[:, 0])
plt.scatter(proj_on_pc1[:, 0], proj_on_pc1[:, 1], c='purple', alpha=0.4, label='第一主成分投影')

plt.xlabel('X1（中心化）')
plt.ylabel('X2（中心化）')
plt.title('主成分方向与数据投影')
plt.legend()
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()
