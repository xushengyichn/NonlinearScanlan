function [ amp,fre ] = ee( data,dt )
%EE Summary of this function goes here
%   Detailed explanation goes here
[normalizeddata,amp]=splinenormalizeep(data);
frecarrier=gradient(normalizeddata,dt);
[normalizedfrecarrier,ampfrecarrier]=splinenormalizeep(frecarrier);
fre=ampfrecarrier/2/pi;
end