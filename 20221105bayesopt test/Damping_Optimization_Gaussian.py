'''
Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
Date: 2022-11-09 11:47:20
LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
LastEditTime: 2022-11-11 15:35:01
FilePath: \20221105bayesopt test\Damping_Optimization_Gaussian.py
Description: 测试python高斯过程回归

Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
'''

#%% 1. 导入库
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'qt')
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
from sklearn.model_selection import cross_val_score,KFold,cross_validate
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import r2_score
import scipy.io
from scipy.io import savemat


#%% 2. 定义函数
eng = matlab.engine.start_matlab()
def black_box_function(fTMD1, xTMD1):
    fTMD1=fTMD1*0.5+0.7
    xTMD1=xTMD1*660
    mode_numbers=1.
    numberofTMD=1.
    calmodes_all=1.
    mu=0.02
    min_damping=eng.Optim_Damping_for_n_foces_n_modes_bayesopt2(mode_numbers,numberofTMD,fTMD1,xTMD1,calmodes_all,mu)
    return min_damping

#%% 3. 定义优化器
optimizer = BayesianOptimization(
    f=black_box_function,
    pbounds={"fTMD1": (0, 1), "xTMD1": (0, 1)},
    verbose=2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
    random_state=1,
)

optimizer.maximize(
    init_points=100,
    n_iter=100,
    acq="ucb", 
    # kappa=0.5112018360898254,
    # What follows are GP regressor parameters
    # kernel=Matern(nu=0.5),
    # alpha=3.554704749334242e-05,
    # alpha=1e-5,
    # normalize_y=True,
    # n_restarts_optimizer=5,
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

# %% 提取计算数据
res = optimizer.res
x_ = np.array([r["params"]['fTMD1'] for r in res])
y_ = np.array([r["params"]['xTMD1'] for r in res])
z_ = np.array([r["target"] for r in res])
xy_ = [x_,y_]
xy_=np.transpose(xy_)
xy_=np.asarray(xy_)
X_, Y_ = np.meshgrid(x_, y_)
Z_est = optimizer._gp.predict(xy_)
# cv=KFold(n_splits=5,shuffle=True,random_state=0)
# result=cross_val_score(optimizer._gp,xy_,z_,cv=cv,scoring='neg_mean_squared_error')
# # loss=result["test_rmse"]
# loss=np.mean(result)
#交叉验证

mat = scipy.io.loadmat('results.mat')
cases=mat['cases']
results=mat['results']
x_cases=(cases[:,0]-0.7)/0.5
y_cases=cases[:,1]/660
xy_ = [x_cases,y_cases]
xy_=np.transpose(xy_)
xy_=np.asarray(xy_)
Z_est_randon = optimizer._gp.predict(xy_).reshape(-1)
z_cal=results.reshape(-1)*-1
r2=r2_score(Z_est_randon,z_cal)

# max_loc=z_.argmax()
# x_opt_=x_[max_loc]
# y_opt_=y_[max_loc]
# variable_percent=1/100
# np.random.seed(1)
# x_randon=x_opt_+(np.random.rand(100)-0.5)*variable_percent
# y_randon=y_opt_+(np.random.rand(100)-0.5)*variable_percent
# # y_randon=np.random.rand(100)
# xy_randon = [x_randon,y_randon]
# xy_randon=np.transpose(xy_randon)
# Z_est_randon = optimizer._gp.predict(xy_randon)
# z_randon=np.zeros((100,1))
# for k1 in range(xy_randon.shape[0]):
#     z_randon[k1]=black_box_function(xy_[k1,0],xy_[k1,1])
    
# Z_est_randon=Z_est_randon.reshape(-1,1)
# z_randon=z_randon.reshape(-1,1)   
# r2=r2_score(z_randon,Z_est_randon)


# %%


labels_top = ["target"]  + ["masked target"]
labels_bot = ["target estimate"] +  ["acqusition function"]
labels = [labels_top, labels_bot]




# Evaluate the actual functions on the grid

n_plots_per_row=2
fig, axs = plt.subplots(2, n_plots_per_row, constrained_layout=True, figsize=(12,8))
for i in range(2):
    for j in range(n_plots_per_row):
        axs[i, j].set_aspect("equal")
        axs[i, j].set_title(labels[i][j])
    

# Extract & unpack the optimization results
max_ = optimizer.max
res = optimizer.res
x_ = np.array([r["params"]['fTMD1'] for r in res])
y_ = np.array([r["params"]['xTMD1'] for r in res])
z_ = np.array([r["target"] for r in res])
xy_ = [x_,y_]
xy_=np.transpose(xy_)
xy_=np.asarray(xy_)

pbounds = {'fTMD1': (0, 1), 'xTMD1': (0, 1)}
x = np.linspace(pbounds['fTMD1'][0], pbounds['fTMD1'][1], 100)
y = np.linspace(pbounds['xTMD1'][0], pbounds['xTMD1'][1], 100)
xy = np.array([[x_i, y_j] for y_j in y for x_i in x])
X, Y = np.meshgrid(x, y)
Z_est = optimizer._gp.predict(xy).reshape(X.shape)


# test=np.array([[0.811,383.7]])


for i in range(2):
    for j in range(n_plots_per_row):
        axs[i,j].scatter(x_, y_, c='white', s=80, edgecolors='black')
        axs[i,j].scatter(max_["params"]['fTMD1'], max_["params"]['xTMD1'], s=80, c='green', edgecolors='black')


Z_est_min=np.min(Z_est)
z_min=np.min(z_)
Z_est_max=np.max(Z_est)
z_max=np.max(z_)

target_vbounds = np.min([Z_est_min, z_min]), np.max([Z_est_max, z_max])


# %%
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf=ax.scatter(x_, y_, z_,c='white', s=80, edgecolors='black', vmin=target_vbounds[0], vmax=target_vbounds[1])



fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf2=ax.plot_surface(X, Y, Z_est,cmap=plt.cm.coolwarm, vmin=target_vbounds[0], vmax=target_vbounds[1])

# %%


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.scatter(cases[:,0],cases[:,1],results[:,0]*-1,c=results[:,0],alpha=0.5)

surf2=ax.plot_surface(X*0.5+0.7, Y*660, Z_est,cmap=plt.cm.coolwarm, vmin=target_vbounds[0], vmax=target_vbounds[1],alpha=0.5)

plt.show()


# %%

# z_cal=z_cal.reshape(-1)
# Z_est_randon=Z_est_randon.reshape(-1)
plt.figure()
plt.scatter(z_cal,Z_est_randon,color='b',s= 0.5)
# plt.scatter(z_cal,Z_est_randon,color='b',s= 0.5)
np.savetxt("z_cal.txt",z_cal)
np.savetxt("Z_est_randon.txt",Z_est_randon)


# %%
