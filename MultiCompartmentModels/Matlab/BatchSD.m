function BatchSD(PW,modelType,ParticleID,Diams)
% Parameters
tstop = 70;%simulation duration
len=5000;%Stable[um]
segdensity=50/6;%segment length [um]
initialdel=0; % pulse delay
dt=0.005;
type =modelType; % fiber type 1:Sundt 2:Tigerholm 3:Rattay 4:Schild
for i=1:length(Diams)
    D=Diams(i);
    for j=1:length(PW)
        dura=PW(j);% duration of pulse [ms]
        call_neuron_MATLAB('SD.hoc',ParticleID,segdensity,dt,tstop,len,D,type,dura,initialdel);
    end
end