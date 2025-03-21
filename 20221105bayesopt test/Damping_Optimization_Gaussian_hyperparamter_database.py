'''
Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
Date: 2022-11-09 11:47:20
LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
LastEditTime: 2022-11-11 15:38:56
FilePath: \20221105bayesopt test\Damping_Optimization_Gaussian_hyperparamter_database.py
Description: 测试python高斯过程回归，优化超参数，数据库

Copyright (c) 2022 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved. 
'''

#%% 1. 导入库
import logging
import sys
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
from sklearn.metrics import r2_score
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
import scipy.io

#%% 1. 定义函数
def objective(mu,kappa,alpha,init_points,n_iter):
    #评估器 
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

    optimizer = BayesianOptimization(
        f=black_box_function,
        pbounds={"fTMD1": (0, 1), "xTMD1": (0, 1)},
        verbose=0,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
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
    x_ = np.array([r["params"]['fTMD1'] for r in res])
    y_ = np.array([r["params"]['xTMD1'] for r in res])
    z_ = np.array([r["target"] for r in res])
    xy_ = [x_,y_]
    xy_=np.transpose(xy_)
    xy_=np.asarray(xy_)
    X_, Y_ = np.meshgrid(x_, y_)
    Z_est = optimizer._gp.predict(xy_)
    #交叉验证
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
        
    # r2=r2_score(z_randon,Z_est_randon)  
    
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

    loss=r2 
    return loss
    

#%% 2. 调用函数
# mu=2.5
# kappa=2.576
# alpha=1e-6
# init_points=100
# n_iter=100

# loss=objective(mu,kappa,alpha,init_points,n_iter)
# %% optuna
def obj(trial):
    mu_ = trial.suggest_float('mu_', 0.0,2.0,step=1.0)
    kappa = trial.suggest_float('kappa', 0.1, 10)
    alpha = trial.suggest_float('alpha', 1e-10, 1e-2,log=True)
    init_points = 100
    # init_points = trial.suggest_int('init_points', 10, 100)
    n_iter = 100
    # n_iter = trial.suggest_int('n_iter', 10, 100)
    if mu_==0:
        mu=0.5
    elif mu_==1:
        mu=1.5
    else:
        mu=2.5
    loss = objective(mu,kappa,alpha,init_points,n_iter)
    return loss

optuna.logging.get_logger("optuna").addHandler(logging.StreamHandler(sys.stdout))
study_name = "example-study"  # Unique identifier of the study.
storage_name = "sqlite:///{}.db".format(study_name)
study = optuna.create_study(direction='maximize',study_name=study_name, storage=storage_name, load_if_exists=True)
study.optimize(obj, n_trials=50)

print(study.best_params)
print(study.best_value)
# %%
df = study.trials_dataframe(attrs=("number", "value", "params", "state"))
print(df)
print("Best params: ", study.best_params)
print("Best value: ", study.best_value)
print("Best Trial: ", study.best_trial)
print("Trials: ", study.trials)


plot_optimization_history(study)
plot_parallel_coordinate(study)

# plot_parallel_coordinate(study, params=["bagging_freq", "bagging_fraction"])
plot_contour(study)
plot_param_importances(study)
# %%
