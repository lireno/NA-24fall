import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

# 输入数据
days = np.array([0, 6, 10, 13, 17, 20, 28])
sp1_weights = np.array([6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7])
sp2_weights = np.array([6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89])

# 生成三次样条插值
cs_sp1 = CubicSpline(days, sp1_weights)
cs_sp2 = CubicSpline(days, sp2_weights)

# 预测未来15天后的值 (x = 43)
x_future = 43
sp1_future_value = cs_sp1(x_future)
sp2_future_value = cs_sp2(x_future)

# 打印结果
print(f"Sample 1 at x=43: {sp1_future_value}")
print(f"Sample 2 at x=43: {sp2_future_value}")

# 绘制插值结果
x_new = np.linspace(0, 43, 200)
y_sp1 = cs_sp1(x_new)
y_sp2 = cs_sp2(x_new)

plt.plot(days, sp1_weights, 'o', label='Sample 1 weight')
plt.plot(x_new, y_sp1, label='Sample 1 Cubic Spline')
plt.plot(days, sp2_weights, 'o', label='Sample 2 weight')
plt.plot(x_new, y_sp2, label='Sample 2 Cubic Spline')
plt.axvline(x=43, color='gray', linestyle='--', label='Prediction x=43')
plt.legend()
plt.xlabel('Days')
plt.ylabel('Weight')
plt.title('Cubic Spline Interpolation')
plt.show()
