function Master(i)
% % Master is used to gather the data for SD, CV, and PP analysis
% % ParticleID - number denoting which particle we are running
% % potentially can move parfor loop out if not possible to do it in here
%% make directories
if exist('PairPulse')~=7
    mkdir('PairPulse')
end
if exist('SD')~=7
    mkdir('SD')
end
if exist('CV')~=7
    mkdir('CV')
end
if exist('VoltageTrace')~=7
    mkdir('VoltageTrace')
end
%% Set parameters
% list of pulse widths in ms
PW = [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10];
% list of diameters in um
D = [0.5 1.5];
ParticleID=1;
% list of IPI in ms
IPI = [1:2:20 22:10:102 120:80:520 1000];
% modelType 1:Sundt 2:Tigerholm 3:Rattay 4:Schild97 5:Schild94
modelType=i;
BatchCV(modelType,ParticleID,D);
BatchSD(PW,modelType,ParticleID,D);
BatchPP(IPI,modelType,ParticleID,D);
end
