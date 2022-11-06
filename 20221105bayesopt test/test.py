import time
import matlab
import matlab.engine
eng = matlab.engine.start_matlab()



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
min_damping=test(trial)

print(min_damping)
eng.exit()