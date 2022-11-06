clc
clear
close all

ras = @(x, y) 20 + x.^2 + y.^2 - 10*(cos(2*pi*x) + cos(2*pi*y));
ras2 = @(x, y) 20 + x^2 + y^2 - 10*(cos(2*pi*x) + cos(2*pi*y));
rf3 = @(x, y) ras(x/10, y/10);
fsurf(rf3,[-30 30],"ShowContours","on")
title("rastriginsfcn([x/10,y/10])")
xlabel("x")
ylabel("y")


fun=@(z)ras3(z.xx,z.yy);

xx = optimizableVariable('xx',[-30 30],'Type','real');
yy = optimizableVariable('yy',[-30 30],'Type','real'); 


% results = bayesopt(fun,[xx yy],'AcquisitionFunctionName','expected-improvement-per-second-plus','ExplorationRatio',0.5,'MaxObjectiveEvaluations',10000,'NumSeedPoints',1000);

results = bayesopt(fun,[xx yy]);

% function fval = myrastrig(in)
% x(1) = in.xx;
% x(2) = in.yy;
% fval = rastriginsfcn(x/10);
% end


function result=ras3(xx,yy)
    result=20 + xx.^2 + yy.^2 - 10*(cos(2*pi*xx) + cos(2*pi*yy));
end

% fun=@(x)Optim_Damping_for_n_foces_n_modes_bayesopt2(mode_numbers,numberofTMD,x.fTMD1,x.xTMD1,calmodes_all,mu);