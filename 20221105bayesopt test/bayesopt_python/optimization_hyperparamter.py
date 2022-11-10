'''
Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
Date: 2022-11-09 11:47:20
LastEditors: Shengyi xushengyichn@outlook.com
LastEditTime: 2022-11-10 01:07:47
FilePath: \20221105bayesopt test\bayesopt_python\optimization_hyperparamter.py
Description: 测试python高斯过程回归，优化超参数

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
from sklearn.model_selection import cross_val_score,KFold,cross_validate
from sklearn.linear_model import LogisticRegression
import optuna

# You can use Matplotlib instead of Plotly for visualization by simply replacing `optuna.visualization` with
# `optuna.visualization.matplotlib` in the following examples.
from optuna.visualization import plot_contour
from optuna.visualization import plot_edf
from optuna.visualization import plot_intermediate_values
from optuna.visualization import plot_optimization_history
from optuna.visualization import plot_parallel_coordinate
from optuna.visualization import plot_param_importances
from optuna.visualization import plot_slice

#%% 1. 定义函数
def objective(mu,kappa,alpha,init_points,n_iter):
    #评估器 
    eng = matlab.engine.start_matlab()
    def black_box_function(x, y):
        z = eng.ras(x, y)
        return z
    optimizer = BayesianOptimization(
        f=black_box_function,
        pbounds={"x": (-2, 2), "y": (-2, 2)},
        verbose=2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state=1,
    )
    optimizer.maximize(
        init_points=init_points,
        n_iter=n_iter,
        acq="ucb", 
        kappa=kappa,
        # What follows are GP regressor parameters
        kernel=Matern(nu=mu),
        alpha=alpha,
        normalize_y=True,
        n_restarts_optimizer=5,
    )
    res = optimizer.res
    x_ = np.array([r["params"]['x'] for r in res])
    y_ = np.array([r["params"]['y'] for r in res])
    z_ = np.array([r["target"] for r in res])
    xy_ = [x_,y_]
    xy_=np.transpose(xy_)
    xy_=np.asarray(xy_)
    X_, Y_ = np.meshgrid(x_, y_)
    Z_est = optimizer._gp.predict(xy_)
    #交叉验证
    cv=KFold(n_splits=5,shuffle=True,random_state=0)
    result=cross_val_score(optimizer._gp,xy_,z_,cv=cv,scoring='neg_mean_squared_error')
    # loss=result["test_rmse"] 
    #交叉验证的结果
    loss=np.mean(result)    
    return loss
    

#%% 2. 调用函数
mu=2.5
kappa=2.576
alpha=1e-6
init_points=100
n_iter=100

loss=objective(mu,kappa,alpha,init_points,n_iter)
# %% optuna
def obj(trial):
    mu = trial.suggest_uniform('mu', 0.5, 2.5)
    kappa = trial.suggest_uniform('kappa', 0.5, 10)
    alpha = trial.suggest_uniform('alpha', 1e-10, 1e-2)
    init_points = trial.suggest_int('init_points', 10, 1000)
    n_iter = trial.suggest_int('n_iter', 10, 1000)
    loss = objective(mu,kappa,alpha,init_points,n_iter)
    return loss

study = optuna.create_study(direction='maximize')
study.optimize(obj, n_trials=500)

print(study.best_params)
print(study.best_value)
# %%
