function BatchCV(type,ParticleID,Diams)
% Parameters
tstop = 70;%simulation duration
len=5000;%Stable[um]
segdensity=50/6;%segment length [um]
initialdel=0; % pulse delay
dura=0.1; % duration of pulse [ms]
for D=Diams
    call_neuron_MATLAB('CV.hoc',ParticleID,segdensity,tstop,len,D,type,dura,initialdel);
end
