function [f,G]=fft_function(x,dt);

%Fourier transfom of discrete-time signal
%--------------------------------------------------------------------------
% (n1 = number of signals)
% (n2 = number of time steps)
%
% OUTPUT
% G = FFT [n1,n2]
% f = Frequency vector in Hz [1,n2]
%
%
% INPUT
% x = Time signal of row vectors [n1,n2]
% dt = Discrete time step of x

% check dimension of x
[n1 n2]=size(x);
if n2<n1 % make x row vectors
	x=x.';
	transform=true;
else
	transform=false;
end	

% make length of x odd numbered
% if mod(length(x),2)==0
%     x=x(:,1:end-1);
% end
[n1 n2]=size(x);


% compute fft
G=zeros(size(x));
for k=1:n1
G(k,:)=fftshift(fft(x(k,:),n2))/n2;
end

% create frequency vector
Fs=1/dt;
f = Fs/2*linspace(-1,1,n2);
