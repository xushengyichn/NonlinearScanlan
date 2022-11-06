mode_numbers=1.;
numberofTMD=1.;
fTMD1=0.83;
xTMD1=55.;
calmodes_all=1.;
mu=0.02;
min_damping=Optim_Damping_for_n_foces_n_modes_bayesopt2(mode_numbers,numberofTMD,fTMD1,xTMD1,calmodes_all,mu);
disp(min_damping)