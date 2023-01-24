'''
Author: Shengyi xushengyichn@outlook.com
Date: 2022-11-06 23:29:11
LastEditors: Shengyi xushengyichn@outlook.com
LastEditTime: 2022-11-07 00:08:19
FilePath: \20221105bayesopt test\test.py
# Description: 在python中调用matlab的函数，注意需要将vs code的工作目录设置为当前路径

Copyright (c) 2022 by Shengyi xushengyichn@outlook.com, All Rights Reserved. 
'''
import time
import matlab
import matlab.engine
eng = matlab.engine.start_matlab()


# 启动一个新的MATLAB进程，并返回Python的一个变量，它是一个MatlabEngine对象，用于与MATLAB过程进行通信。
d = eng.multiplication_matlab(3.,2.) # 可以调用matlab写的脚本函数




mode_numbers=1.
numberofTMD=1.
fTMD1=0.83
xTMD1=55.
calmodes_all=1.
mu=0.02


def test(trial):
    min_damping=eng.Optim_Damping_for_n_foces_n_modes_bayesopt2(mode_numbers,numberofTMD,fTMD1,xTMD1,calmodes_all,mu)
    return min_damping

def objective(trial):
    fTMD1 = trial.suggest_float("fTMD1", 0.7, 10)
    xTMD1 = trial.suggest_float("xTMD1", 0, 660)
    min_damping=eng.Optim_Damping_for_n_foces_n_modes_bayesopt2(mode_numbers,numberofTMD,fTMD1,xTMD1,calmodes_all,mu)
    return min_damping
    


# min_damping=eng.Optim_Damping_for_n_foces_n_modes_bayesopt2(mode_numbers,numberofTMD,fTMD1,xTMD1,calmodes_all,mu)
trial=1
start =time.perf_counter()
min_damping=test(trial)
end = time.perf_counter()
print('Running time: %s Seconds'%(end-start))
print(min_damping)
eng.exit()