function [u udot u2dot] = NewmarkInt(t,M,C,K,P,gam,beta,u0,udot0)

%Newmark's Direct Integration Method
%--------------------------------------------------------------------------
% OUTPUT
% u =       Displacemente Response   [n,2]
% (n = number of time steps)
% (ndof = number degrees of freedom)

% INPUT
% t =       Time vector         [1,n]
% M =       mass matrix         [ndof,ndof]
% C =       damping matrix      [ndof,ndof]
% K =       stiffness matrix    [ndof,ndof]
% P =       load vs. time       [ndof,n]
% gam =     gamma (constant)
% beta =    beta  (constant)
%u0 =       Initial displacements
%udot0 =    Initial velocity

%--------------------------------------------------------------------------
% beta = 0,     gamma = 1/2 -> explicit central difference method
% beta = 1/4,   gamma = 1/2 -> undamped trapezoidal rule (implicit)

%% 1.0 Initial calculations
%1.1
u=zeros(size(M,1),length(t));
u(:,1)=u0;
udot=udot0;
u2dot = M\(P(:,1)-C*udot0-K*u0); 

%1.2
dt = t(2) - t(1);                                   
%1.3
kgor = K + gam/(beta*dt)*C + M*1/(beta*dt^2);      
%1.4
a = M/(beta*dt) + gam/beta*C;                      
b = 0.5*M/beta + dt*(0.5*gam/beta - 1)*C;


%% 2.0 Calculations for each time step i

for i = 1:(length(t)-1)
    %2.1
    deltaP = P(:,i+1)-P(:,i) + a*udot(:,i) + b*u2dot(:,i);          
 
    %2.2
    du_i = kgor\deltaP;
    %2.3
    dudot_i = gam/(beta*dt)*du_i - gam/beta*udot(:,i) + dt*(1-0.5*gam/beta)*u2dot(:,i);
    %2.4
    du2dot_i = 1/(beta*dt^2)*du_i - 1/(beta*dt)*udot(:,i) - 0.5/beta*u2dot(:,i);
    %2.5
    u(:,i+1) = du_i + u(:,i);
    udot(:,i+1) = dudot_i + udot(:,i);
    u2dot(:,i+1) = du2dot_i + u2dot(:,i);
end

% u = u';
% udot=udot';
% u2dot=u2dot';