clc
clear
close all
nodeondeck = importdata('nodeondeck.txt');
nodegap = importdata('nodegap.txt');
nodeondecknew=[];
intalledlocation=[];
for k1=1:length(nodeondeck)
    if mod(k1,5)==0
        nodeondecknew=[nodeondecknew nodeondeck(k1)];
        intalledlocation=[intalledlocation nodegap(k1)];
    end
end
