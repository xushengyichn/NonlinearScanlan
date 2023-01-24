# %%
import optuna
import time
import matlab
import matlab.engine

# You can use Matplotlib instead of Plotly for visualization by simply replacing `optuna.visualization` with
# `optuna.visualization.matplotlib` in the following examples.
from optuna.visualization import plot_contour
from optuna.visualization import plot_edf
from optuna.visualization import plot_intermediate_values
from optuna.visualization import plot_optimization_history
from optuna.visualization import plot_parallel_coordinate
from optuna.visualization import plot_param_importances
from optuna.visualization import plot_slice



# %%

eng = matlab.engine.start_matlab()

 
mode_numbers=1.
numberofTMD=1.
# fTMD1=0.83
# xTMD1=55.
calmodes_all=1.
mu=0.02

def objective(trial):
    fTMD1 = trial.suggest_float("fTMD1", 0.7, 1.2)
    xTMD1 = trial.suggest_float("xTMD1", 0, 660)
    min_damping=eng.Optim_Damping_for_n_foces_n_modes_bayesopt2(mode_numbers,numberofTMD,fTMD1,xTMD1,calmodes_all,mu)
    return min_damping
 
study = optuna.create_study(direction='minimize')
study.optimize(objective, n_trials=500)
 
print(study.best_params)
print(study.best_value)

# %%


# %%
optuna.visualization.matplotlib.plot_contour(study, params=['fTMD1', 'xTMD1'])

# %%
optuna.visualization.matplotlib.plot_param_importances(study)

# %%
optuna.visualization.matplotlib.plot_optimization_history(study)

# %%
optuna.visualization.matplotlib.plot_intermediate_values(study)

# %%
optuna.visualization.matplotlib.plot_parallel_coordinate(study)

# %%
optuna.visualization.matplotlib.plot_contour(study)

# %%
optuna.visualization.matplotlib.plot_slice(study)



# %%
