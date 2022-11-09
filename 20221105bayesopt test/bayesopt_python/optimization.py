'''
Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
Date: 2022-11-09 11:47:20
LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
LastEditTime: 2022-11-09 20:56:01
FilePath: \bayesopt_python\optimization.py
Description: 测试python高斯过程回归

Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
'''

#%% 1. 导入库
import time
import matlab
import matlab.engine
from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from sklearn.gaussian_process.kernels import Matern
from sklearn.gaussian_process import GaussianProcessRegressor
#%% 2. 定义函数
eng = matlab.engine.start_matlab()
def black_box_function(x, y):
    z = eng.ras(x, y)
    return z

#%% 3. 定义优化器
optimizer = BayesianOptimization(
    f=black_box_function,
    pbounds={"x": (-2, 2), "y": (-2, 2)},
    verbose=2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
    random_state=1,
)

optimizer.maximize(
    init_points=100,
    n_iter=100,
    acq="ucb", 
    kappa=2.576,
    # What follows are GP regressor parameters
    kernel=Matern(nu=2.5),
    alpha=1e-6,
    normalize_y=True,
    n_restarts_optimizer=5,
)

#%% 4. 绘制结果
# x= np.linspace(-2, 2, 1000).reshape(-1, 1)
# y= np.linspace(-2, 2, 1000).reshape(-1, 1)
# xv,yv=np.meshgrid(x,y)
# all=np.hstack((xv.reshape(-1,1),yv.reshape(-1,1)))

# xx=all[:,0].reshape(-1,1)
# yy=all[:,1].reshape(-1,1)

# zz=np.zeros((xx.size,1))
# # ax.scatter(xx, yy)
# # for k1 in range(xx.size):
# #     zz[k1]=black_box_function(xx[k1],yy[k1])
# def target_function(xx, yy):
#     zz=-(20 + np.square(xx) + np.square(yy) - 10*(np.cos(2*np.pi*xx) + np.cos(2*np.pi*yy)))
#     return zz

# zz=target_function(xx,yy)
# # c=-(20 + xx^2 + yy^2 - 10*(np.cos(2*np.pi*xx) + np.cos(2*np.pi*yy)))
# # z=black_box_function(xx,yy)    
# # z=black_box_function(xx,yy)
# # z = black_box_function(x, y)
# fig = plt.figure(figsize=(12, 12))
# ax = fig.add_subplot(projection='3d')
# zz=np.asarray(zz)
# # plt.scatter(x, z)
# ax.scatter(xx, yy, zz,c=zz)
# print("finish")

# %%

# %%


labels_top = ["target"]  + ["masked target"]
labels_bot = ["target estimate"] +  ["acqusition function"]
labels = [labels_top, labels_bot]

def target_function(xx, yy):
    zz=-(20 + np.square(xx) + np.square(yy) - 10*(np.cos(2*np.pi*xx) + np.cos(2*np.pi*yy)))
    return zz

pbounds = {'x': (-2, 2), 'y': (-2, 2)}
x = np.linspace(pbounds['x'][0], pbounds['x'][1], 1000)
y = np.linspace(pbounds['y'][0], pbounds['y'][1], 1000)
xy = np.array([[x_i, y_j] for y_j in y for x_i in x])
X, Y = np.meshgrid(x, y)
# Evaluate the actual functions on the grid
Z = target_function(X, Y)
n_plots_per_row=2
fig, axs = plt.subplots(2, n_plots_per_row, constrained_layout=True, figsize=(12,8))
for i in range(2):
    for j in range(n_plots_per_row):
        axs[i, j].set_aspect("equal")
        axs[i, j].set_title(labels[i][j])
    

# Extract & unpack the optimization results
max_ = optimizer.max
res = optimizer.res
x_ = np.array([r["params"]['x'] for r in res])
y_ = np.array([r["params"]['y'] for r in res])

Z_est = optimizer._gp.predict(xy).reshape(Z.shape)
target_vbounds = np.min([Z, Z_est]), np.max([Z, Z_est])

axs[0,0].contourf(X, Y, Z, cmap=plt.cm.coolwarm, vmin=target_vbounds[0], vmax=target_vbounds[1])
axs[1,0].contourf(X, Y, Z_est, cmap=plt.cm.coolwarm, vmin=target_vbounds[0], vmax=target_vbounds[1])



for i in range(2):
    for j in range(n_plots_per_row):
        axs[i,j].scatter(x_, y_, c='white', s=80, edgecolors='black')
        axs[i,j].scatter(max_["params"]['x'], max_["params"]['y'], s=80, c='green', edgecolors='black')


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf=ax.plot_surface(X, Y, Z, cmap=plt.cm.coolwarm, vmin=target_vbounds[0], vmax=target_vbounds[1])

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf=ax.plot_surface(X, Y, Z_est, cmap=plt.cm.coolwarm, vmin=target_vbounds[0], vmax=target_vbounds[1])

