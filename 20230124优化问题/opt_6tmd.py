# %%
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

# %%
def objective(mu,xi,alpha,init_points,n_iter):
    #%% 2. 定义函数
    def black_box_function(mTMD1,mTMD2,mTMD3,mTMD4,mTMD5,fTMD1,fTMD2,fTMD3,fTMD4,fTMD5,fTMD6,dTMD1,dTMD2,dTMD3,dTMD4,dTMD5,dTMD6,xTMD1,xTMD2,xTMD3,xTMD4,xTMD5,xTMD6):
        eng = matlab.engine.start_matlab()
        total_tmd_mass_ratio = 0.02 # 总质量比 The total mass ratio
        mass_six_span = 10007779.7 # 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
        total_tmd_mass = total_tmd_mass_ratio * mass_six_span # 总质量 The total mass
        single_tmd_mass=total_tmd_mass/6 # 单个TMD质量 The mass of a single TMD
        mTMD1=0.5*single_tmd_mass+mTMD1*0.5*single_tmd_mass # 质量 The mass mTMD1
        mTMD2=0.5*single_tmd_mass+mTMD2*0.5*single_tmd_mass # 质量 The mass mTMD2
        mTMD3=0.5*single_tmd_mass+mTMD3*0.5*single_tmd_mass # 质量 The mass mTMD3
        mTMD4=0.5*single_tmd_mass+mTMD4*0.5*single_tmd_mass # 质量 The mass mTMD4
        mTMD5=0.5*single_tmd_mass+mTMD5*0.5*single_tmd_mass # 质量 The mass mTMD5
        mTMD6= total_tmd_mass - mTMD1 - mTMD2 - mTMD3 - mTMD4 - mTMD5 # 质量 The mass mTMD6
        
        fTMD1=0.7+fTMD1*1.3 # 频率 The frequency fTMD1
        fTMD2=0.7+fTMD2*1.3 # 频率 The frequency fTMD2
        fTMD3=0.7+fTMD3*1.3 # 频率 The frequency fTMD3
        fTMD4=0.7+fTMD4*1.3 # 频率 The frequency fTMD4
        fTMD5=0.7+fTMD5*1.3 # 频率 The frequency fTMD5
        fTMD6=0.7+fTMD6*1.3 # 频率 The frequency fTMD6
        
        dTMD1=0.05+dTMD1*0.15 # 阻尼比 The damping ratio dTMD1
        dTMD2=0.05+dTMD2*0.15 # 阻尼比 The damping ratio dTMD2
        dTMD3=0.05+dTMD3*0.15 # 阻尼比 The damping ratio dTMD3
        dTMD4=0.05+dTMD4*0.15 # 阻尼比 The damping ratio dTMD4
        dTMD5=0.05+dTMD5*0.15 # 阻尼比 The damping ratio dTMD5
        dTMD6=0.05+dTMD6*0.15 # 阻尼比 The damping ratio dTMD6
        
        
        xTMD1=xTMD1*660 # TMD1的x坐标 The x-coordinate of TMD1
        xTMD2=xTMD2*660 # TMD2的x坐标 The x-coordinate of TMD2
        xTMD3=xTMD3*660 # TMD3的x坐标 The x-coordinate of TMD3
        xTMD4=xTMD4*660 # TMD4的x坐标 The x-coordinate of TMD4
        xTMD5=xTMD5*660 # TMD5的x坐标 The x-coordinate of TMD5
        xTMD6=xTMD6*660 # TMD6的x坐标 The x-coordinate of TMD6
        
        
        t_length=matlab.double(100) # 时间长度 The time length
        number_of_modes_to_control=matlab.double([1,2,3,4,5,6]) # 控制模态 The controlled modes
        number_of_modes_to_consider=10 # 考虑模态 The considered modes
        number_of_tmds=6 # TMD数量 The number of TMDs
        modal_damping_ratios=np.ones((1,number_of_modes_to_consider))*0.003 # 模态阻尼比 The modal damping ratios
        
        result = -eng.b_0_6_tmd(number_of_modes_to_control,number_of_modes_to_consider,number_of_tmds,modal_damping_ratios,t_length,mTMD1,mTMD2,mTMD3,mTMD4,mTMD5,fTMD1,fTMD2,fTMD3,fTMD4,fTMD5,fTMD6,dTMD1,dTMD2,dTMD3,dTMD4,dTMD5,dTMD6,xTMD1,xTMD2,xTMD3,xTMD4,xTMD5,xTMD6,total_tmd_mass)
        return result
        #%% 3. 定义优化器
        #%% 3. 定义优化器
        
    optimizer = BayesianOptimization(
    f=black_box_function,
    pbounds={"mTMD1": (0,1), 
            "mTMD2": (0,1),
            "mTMD3": (0,1),
            "mTMD4": (0,1),
            "mTMD5": (0,1),
            "fTMD1": (0,1),
            "fTMD2": (0,1),
            "fTMD3": (0,1),
            "fTMD4": (0,1),
            "fTMD5": (0,1),
            "fTMD6": (0,1),
            "dTMD1": (0,1),
            "dTMD2": (0,1),
            "dTMD3": (0,1),
            "dTMD4": (0,1),
            "dTMD5": (0,1),
            "dTMD6": (0,1),
            "xTMD1": (0,1),
            "xTMD2": (0,1),
            "xTMD3": (0,1),
            "xTMD4": (0,1),
            "xTMD5": (0,1),
            "xTMD6": (0,1),
            },
    verbose=2,  # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
    random_state=1,
    )
    optimizer.probe(
    params={ "mTMD1": 1, 
             "mTMD2": 1,
             "mTMD3": 1,
             "mTMD4": 1,
             "mTMD5": 1,
             "fTMD1": 0.1,
             "fTMD2": 0.153846,
             "fTMD3": 0.276923,
             "fTMD4": 0.446153,
             "fTMD5": 0.623077,
             "fTMD6": 0.769231,
             "dTMD1": 0.2,
             "dTMD2": 0.2,
             "dTMD3": 0.2,
             "dTMD4": 0.2,
             "dTMD5": 0.2,
             "dTMD6": 0.2,
             "xTMD1": 0.4181818,
             "xTMD2": 0.9181818,
             "xTMD3": 0.9212121,
             "xTMD4": 0.5969697,
             "xTMD5": 0.2515152,
             "xTMD6": 0.5848485,},
    lazy=True,
)

    optimizer.maximize(
        init_points=init_points,
        n_iter=n_iter,
        acq="ei",  # Expected Improvement.
        # kappa=7.2088954429132,
        xi=xi,
        # What follows are GP regressor parameters
        kernel=Matern(nu=mu),
        # alpha=1,
        # alpha=2.2492215663634233e-10,
        alpha=alpha,
        normalize_y=True,
        n_restarts_optimizer=5,
    )
    
    loss=optimizer.max['target']
    return loss

        
        
    

# %%
#%% 2. 调用函数

def obj(trial):
    mu_ = trial.suggest_float('mu_', 0.0,2.0,step=1.0)
    kappa = trial.suggest_float('xi', 1e-4, 1e-1,log=True)
    alpha = trial.suggest_float('alpha', 1e-10, 1e-2,log=True)
    init_points = 10
    # init_points = trial.suggest_int('init_points', 10, 100)
    n_iter = 500
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
study_name = "opt_6tmd"  # Unique identifier of the study.
storage_name = "sqlite:///{}.db".format(study_name)
study = optuna.create_study(direction='maximize',study_name=study_name, storage=storage_name, load_if_exists=True)
study.optimize(obj, n_trials=50)

print(study.best_params)
print(study.best_value)

# %%
optuna.logging.get_logger("optuna").addHandler(logging.StreamHandler(sys.stdout))
study_name = "example-study"  # Unique identifier of the study.
storage_name = "sqlite:///{}.db".format(study_name)
study = optuna.create_study(direction='maximize',study_name=study_name, storage=storage_name, load_if_exists=True)
print(study.best_params)
print(study.best_value)
# print(study.trials)

print(study.directions)



