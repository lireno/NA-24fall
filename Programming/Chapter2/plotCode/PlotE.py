import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 定义指数函数 (exp_func) = a * exp(b * x) + c
def exp_func(x, a, b, c):
    return a * np.exp(b * x) + c

# 输入天数和体重数据
days = [0, 6, 10, 13, 17, 20, 28]
sp1_weights = [6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7]
sp2_weights = [6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89]

# 选取最后5天的数据
days_last5 = np.array([10, 13, 17, 20, 28])
sp1_weights_last5 = np.array([42.7, 37.3, 30.1, 29.3, 28.7])
sp2_weights_last5 = np.array([18.9, 15.0, 10.6, 9.44, 8.89])

# 拟合指数模型到两个物种的体重数据
params_sp1, _ = curve_fit(exp_func, days_last5, sp1_weights_last5)
params_sp2, _ = curve_fit(exp_func, days_last5, sp2_weights_last5, maxfev=2000)

# 生成拟合的值用于绘图
x_range = np.linspace(min(days_last5), 43, 1000)
fitted_curve_sp1 = exp_func(x_range, *params_sp1)
fitted_curve_sp2 = exp_func(x_range, *params_sp2)

# 绘制原始数据和拟合的曲线
plt.figure(figsize=(10, 6))

# 绘制物种1的数据和拟合曲线
plt.scatter(days_last5, sp1_weights_last5, color='blue', label='Species 1 Data (last 5 days)')
plt.plot(x_range, fitted_curve_sp1, label='Sample 1 Fitted Exp. Function', linestyle='--', color='blue')

# 绘制物种2的数据和拟合曲线
plt.scatter(days_last5, sp2_weights_last5, color='green', label='Species 2 Data (last 5 days)')
plt.plot(x_range, fitted_curve_sp2, label='Sample 2 Fitted Exp. Function', linestyle='--', color='green')

# 添加标题和标签
plt.title("Exponential Fit")
plt.xlabel("Days")
plt.ylabel("Weight (g)")
plt.legend()
plt.grid(True)
plt.show()

# 输出拟合的参数
print("Species 1 Fitted Parameters (a, b, c):", params_sp1)
print("Species 2 Fitted Parameters (a, b, c):", params_sp2)
