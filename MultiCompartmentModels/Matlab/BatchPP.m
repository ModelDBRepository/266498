function BatchPP(IPIs,modelTypes,ParticleID,Diams)
% Parameters
len=5000;%Stable[um]
segdensity=50/6;%segment length [um]
initialdel=0; % pulse delay
dura=0.1;
tstop=50;
dt=0.005;
for i=1:length(IPIs)
    for j=1:length(modelTypes)
        for k=1:length(Diams)
            IPI = IPIs(i);
            D = Diams(k);
            modelType = modelTypes(j);
            tstop=IPI+tstop;
            delay=IPI; % duration of pulse [ms]
            type =modelType; % fiber type 1:Sundt 2:Tigerholm 3:Rattay
            start=load(['SD/c_fiber_Thresh_' num2str(D*1000) '_D_' num2str(100) '_Duration_' num2str(modelType) '_Type_' num2str(ParticleID) '_Intra.dat']);
            call_neuron_MATLAB('PairPulseTrue.hoc',ParticleID,segdensity,dt,tstop,len,D,type,dura,initialdel,delay,start);
        end
    end
end