# 导入必要的库
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
plt.rcParams['font.sans-serif'] = ['SimHei']                 # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False                   # 用来正常显示负号

def ODE_model(constant,y,t):
    '''
    constant --> 各个反应对应的速率常数和诱导物加入的速率，字典形式
    y --> 因变量（各物质的浓度），列表形式，顺序：DNA , LD , IPTG , IL, ILD, LacI, mRNA, pro
    '''

    # 以下状态变量均为对应符号含义的浓度
    DNA , LD , IPTG , IL, ILD, LacI, mRNA, pro = y

    # 反应速率计算式
    # 1. 反应调控体系
    # (1) 阻遏蛋白与操纵子结合
    v_1_f = constant['k_LacI_on'] * LacI * DNA
    v_1_r = constant['k_LacI_off'] * LD

    # (2) 诱导物与阻遏蛋白结合
    v_2_f = constant['k_IPTG_on_1'] * IPTG * LacI
    v_2_r = constant['k_IPTG_off_1'] * IL

    # (3) 诱导物与LD结合
    v_3_f = constant['k_IPTG_on_2'] * IPTG * LD
    v_3_r = constant['k_IPTG_off_2'] * ILD

    # (4) ILD解离
    v_4 = constant['k_ILD_out'] * ILD

    # 2. 表达反应体系
    # (5) 转录
    v_5 = constant['k_tc'] * DNA

    # (6) 翻译
    v_6 = constant['k_tl'] * mRNA

    # 3. 降解反应体系
    # (7) mRNA降解
    v_7 = constant['k_d_mRNA'] * mRNA

    # (8) β-半乳糖苷酶（pro）降解
    v_8 = constant['k_d_pro'] * pro

    # ODE模型建立
    # 1. DNA
    dDNA_dt = v_4 - (v_1_f - v_1_r)
    dLD_dt = (v_1_f - v_1_r) - (v_3_f - v_3_r)
    dILD_dt = (v_3_f - v_3_r) - v_4

    # 2. 阻遏蛋白
    dLacI_dt = 0 - (v_1_f - v_1_r) - (v_2_f - v_2_r)

    # 3. 诱导物
    dIPTG_dt = 0 - (v_2_f - v_2_r) - (v_3_f - v_3_r)

    dIL_dt = v_4 + (v_2_f - v_2_r)

    # 4. mRNA
    dmRNA_dt = v_5 - v_7

    # 5. β-半乳糖苷酶（pro）表达
    dpro_dt = v_6 - v_8

    # 返回：DNA , LD , IPTG , IL, ILD, LacI, mRNA, pro
    ODE = [
        dDNA_dt,        # DNA
        dLD_dt,         # LacI_DNA
        dIPTG_dt,       # IPTG
        dIL_dt,         # IPTG_LacI
        dILD_dt,        # IPTG_LacI_DNA
        dLacI_dt,       # LacI
        dmRNA_dt,       # mRNA
        dpro_dt,        # β-半乳糖苷酶（pro）
    ]

    return ODE

# 定义常数
constant = {
    # 结合/解离参数
    # 反应1
    'k_LacI_on': 1e7,      # M⁻¹s⁻¹
    'k_LacI_off': 0.01,    # s⁻¹

    # 反应2
    'k_IPTG_on_1': 1e5,      # M⁻¹s⁻¹
    'k_IPTG_off_1': 1.0,     # s⁻¹
    # 反应3
    'k_IPTG_on_2': 1e5,      # M⁻¹s⁻¹
    'k_IPTG_off_2': 1.0,     # s⁻¹

    # 反应4
    'k_ILD_out': 10.0,      # s⁻¹

    # 表达参数
    'k_tc': 0.3,           # s⁻¹
    'k_tl': 2.0,           # s⁻¹

    # 降解参数
    'k_d_mRNA': 0.05,        # s⁻¹ (半衰期~14s)
    'k_d_pro': 1e-5,     # s⁻¹

    # 总浓度
    'DNA_total': 1e-9,     # M (假设1nM)
    'LacI_total': 2e-8,    # M (20nM)
    'RNAP': 3e-7,          # M (300nM)

}

# 初始参数
# DNA , LD , IPTG , IL, ILD, LacI, mRNA, pro
LacI_0 = constant['LacI_total'] - constant['DNA_total']
y0 = [
    0 ,                       # DNA
    constant['DNA_total'],    # LacI_DNA 即DNA总数
    100e-6,                   # IPTG
    0 ,                       # IPTG_LacI，初始时为0
    0 ,                       # IPTG_LacI_DNA = 0
    LacI_0 ,                  # LacI 少量
    0 ,                       # mRNA
    0                         # β-半乳糖苷酶（pro）
]

# 设置时间
t_start = 0
t_end = 1000
t_eval = np.linspace(t_start, t_end, 1000)

# 解ODE
solution = solve_ivp(
    lambda t ,y : ODE_model(constant, y, t),
    t_span=(t_start, t_end),
    y0=y0,
    t_eval=t_eval,
    method='BDF'
)


# 绘制动态响应（单一诱导物浓度）

plt.figure(figsize=(12, 8))

# β-半乳糖苷酶动态响应
plt.subplot(2, 2, 1)
plt.plot(solution.t, solution.y[7], 'b-', linewidth=2)
plt.xlabel('时间 (s)')
plt.ylabel('β-半乳糖苷酶浓度 (M)')
plt.title('酶动态响应')
plt.grid()

# mRNA动态响应
plt.subplot(2, 2, 2)
plt.plot(solution.t, solution.y[6], 'g-', linewidth=2)
plt.xlabel('时间 (s)')
plt.ylabel('mRNA浓度 (M)')
plt.title('mRNA动态响应')
plt.grid()

# DNA结合状态
plt.subplot(2, 2, 3)
plt.plot(solution.t, solution.y[0], label='游离DNA')
plt.plot(solution.t, solution.y[1], label='LD复合物')
plt.plot(solution.t, solution.y[4], label='ILD复合物')
plt.xlabel('时间 (s)')
plt.ylabel('浓度 (M)')
plt.title('DNA结合状态')
plt.legend()
plt.grid()

# 阻遏蛋白状态
plt.subplot(2, 2, 4)
plt.plot(solution.t, solution.y[5], label='游离LacI')
plt.plot(solution.t, solution.y[1], label='LD复合物')
plt.plot(solution.t, solution.y[3], label='IL复合物')
plt.plot(solution.t, solution.y[4], label='ILD复合物')
plt.xlabel('时间 (s)')
plt.ylabel('浓度 (M)')
plt.title('阻遏蛋白状态')
plt.legend()
plt.grid()

plt.tight_layout()


# 绘制不同IPTG浓度的IPTG_pro曲线
# 1. 设置IPTG范围
IPTG_range = np.logspace(-9, -3, 20)

# 2. 定义酶蛋白浓度列表
pro_range = []

# 2. 求解
for x in range (len(IPTG_range)) :
    IPTG = IPTG_range[x]
    y0 =[
        0 ,                       # DNA
        constant['DNA_total'],    # LacI_DNA 即DNA总数
        IPTG,                   # IPTG
        0 ,                       # IPTG_LacI，初始时为0
        0 ,                       # IPTG_LacI_DNA = 0
        LacI_0 ,                  # LacI 少量
        0 ,                       # mRNA
        0                         # β-半乳糖苷酶（pro）
    ]

    solution = solve_ivp(
        lambda t ,y : ODE_model(constant, y, t),
        t_span=(t_start, t_end),
        y0=y0,
        t_eval=t_eval,
        method='BDF'
    )

    pro_range.append(solution.y[7][-1])

plt.figure(figsize=(8, 5))
plt.semilogx(IPTG_range, pro_range, 'ro-', label='β-半乳糖苷酶稳态浓度')

# 负指数的正常显示设置
xticks = [10**i for i in range(-9, -2)]             # 刻度位置：10^-9 到 10^-3
xlabels = [f'10$^{{{i}}}$' for i in range(-9, -2)]  # 用LaTeX语法渲染指数
plt.xticks(xticks, xlabels)

plt.xlabel('IPTG浓度 (M)')
plt.ylabel('酶稳态浓度 (M)')
plt.title('不同诱导物浓度下的酶表达量')
plt.legend()
plt.show()

