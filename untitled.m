clc
clear

theta = optimizableVariable('theta',[1,10],'Type','real');
test = optimizableVariable('test',[1,10],'Type','real');
fun = @(x)objfunxx(x.theta,x.test); %Create Function Handles
results = bayesopt(fun,[theta test])

theta1=results.XAtMinObjective.theta;
test=results.XAtMinObjective.test;
yyy=objfunxx(theta1,test)
disp(yyy)

function y = objfunxx(x,test)
y=x+test;
end