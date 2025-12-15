# 文件说明
这个文件是我学习质量作用定律的笔记和学习案例，案例模型尚不完美，欢迎大佬指导！！！

# 
## 质量作用定律（mass_action_law）
1. **反应速率：**
单位体积反应系统中反应进度随时间的变化率，常用单位时间内反应物或生成物的浓度变化来表示。

2. **核心思想：**
化学反应的速率，与各反应物浓度的幂成正比

3. **适用范围：**
基元反应，即一步反应

4. **形式：**

    如基元反应： $aA + bB ⇌ cC + dD $
    
    则**正向反应速率**：$v_f = k_f * [A]^a * [B]^b $
    - $k_f$：正向反应速率常数
    - $[A]，[B]$：物质A，B的瞬时浓度
    
    反向类似，C，D作为反应物，一般用v_r，k_r表示

    其中，**净反应速率**：$v = v_f - v_r $

5. **应用举例：**
    - 连接动力学与热力学：
    $$
    平衡时有：v_f = v_r
    $$
    则
    $$
    k_f * [A]^a * [B]^b  = k_r * [C]^c * [D]^d
    $$
    $$
    \frac{[C]^c * [D]^d}{[A]^a * [B]^b} = \frac{k_f}{k_r} = K(平衡常数)
    $$

    - 反应微分方程建立
    
    生物化学反应中，有不少是单生成物反应，这些反应就可以根据质量作用定律构建浓度与时间的微分方程式

    如：$ A + B → C$

    则
    $$
    v = \frac{d[C]}{dt} = k_f * [A] * [B]
    $$

# 案例说明
## 1. 题目
请基于乳糖操纵子的调控机理，预测在不同诱导物（IPTG）浓度下，pro基因表达产物β-半乳糖苷酶（pro）的稳态浓度和动态响应。

## 2. 参考假设
1. 细胞体积恒定，分子均匀分布
2. lacZYA认为是同时转录翻译
3. 忽略葡萄糖的调控（CAP）
4. 认为调控反应、表达反应和降解反应分为三个步骤
5. 忽略协同反应
6. 把DNA转录看做是一步反应
7. 将蛋白质的翻译看做是一步反应
8. 将DNA转录和蛋白质翻译看做是先后进行
9. 忽略翻译出的阻遏蛋白的稀释
10. 认为IPTG不被分解


## 二、参数说明

| 参数 | 意义 |  | 参数 | 意义 |
|------|---------|-|------|---------|
|LacI|阻遏蛋白| | IPTG|诱导剂|
|LacI_DNA|阻碍蛋白与DNA结合物| |IPTG_LacI|诱导剂与阻遏蛋白结合物|
|IPTG_LacI_DNA|诱导剂、阻遏蛋白、DNA三者结合物||RNAP|RNA聚合酶|
|RNAP_DNA|RNA聚合酶与DNA结合物||pro|β-半乳糖苷酶|
| k_LacI_on | LacI_DNA结合速率常数 || k_LacI_off | LacI_DNA解离速率常数 |
| k_IPTG_on_1 | IPTG_LacI结合速率常数 || k_IPTG_off_1 | IPTG_LacI解离速率常数 |
|k_IPTG_on_2| IPTG_(LacI_DNA)结合速率常数||k_IPTG_off_2|IPTG_(LacI_DNA)解离速率常数|
|k_ILD_out|IPTG_LacI_DNA|| k_tx | 转录速率常数 |
| k_tl | 翻译速率常数 ||k_d_mRNA|mRNA降解速率常数|
|k_d_pro|β-半乳糖苷酶降解速率常数||$v_{imput}$|诱导物输入速率|

## 三、分析系统的反应组成
### 1. 调控反应体系
#### (1) **阻遏蛋白与操纵子结合**
```
LacI + DNA ⇌ LacI_DNA(LD)
```
- 正向速率：`v₁⁺ = k_LacI_on × [LacI] × [DNA]`
- 逆向速率：`v₁⁻ = k_LacI_off × [LD]`
- 净速率：`v₁ = v₁⁺ - v₁⁻`

#### (2) **诱导物与阻遏蛋白结合**
```
IPTG + LacI ⇌ IPTG_LacI(IL)
```
- 净速率：`v₂ = k_IPTG_on_1 × [IPTG] × [LacI] - k_IPTG_off_1 × [IL]`

#### (3) **诱导物与LD结合**
```
IPTG + LD ⇌ IPTG_LacI_DNA(ILD)
```
- 净速率：`v₃ = k_IPTG_on_2 × [IPTG] × [LD] - k_IPTG_off_2 × [ILD]`

#### (4) **ILD解离**
```
ILD → DNA + IL
```
- 净速率：`v₄ = k_ILD_out × [ILD]`

### 2. 表达反应体系
#### (5) **转录**
```
DNA + RNAP → RNAP_DNA(RD) → DNA + RNAP + mRNA
```
- 净速率：`v₅ = k_tx × [DNA]`

#### (6) **翻译**
```
mRNA (一系列反应) → pro + mRNA
```
- `v₆ = k_tl × [mRNA]`

### 3. 降解
#### (7) **mRNA的降解**
```
mRNA → ∅
```
- 速率：`v₇ = k_d_mRNA × [mRNA]`

#### (8) **β-半乳糖苷酶（pro）的降解**
```
pro → ∅
```
- `v₈ = k_d_pro × [pro]`

## 四、建立方程
### 1. DNA
$$
\frac{d[DNA]}{dt} = v₄ - v₁
$$
$$
\frac{d[LD]}{dt} = v₁ - v₃
$$
$$
\frac{d[ILD]}{dt} = v₃ - v₄
$$

### 2. 阻遏蛋白
$$
\frac{d[LacI]}{dt} = - v₁ - v₂
$$

### 3. 诱导物
$$
\frac{d[IPTG]}{dt} =  - v₂ - v₃ 
$$

$$
\frac{d[IL]}{dt} = v₄ + v₂
$$

### 4. 基因表达
$$
\frac{d[mRNA]}{dt} = v₅ - v₇
$$

### 5. β-半乳糖苷酶（pro）表达
$$
\frac{d[pro]}{dt} = v₆ - v₈
$$


### 6. 约束条件
1. `[DNA_total] = [DNA] + [LD] + [ILD]`
2. `[LacI_total] = [LacI] + [LD] + [IL] + [ILD]`




