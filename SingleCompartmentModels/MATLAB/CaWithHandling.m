% Author: David Catherall; Duke University 
% Created: September 2017
% Description: Matlab Voltage Clamp Implementation for Schild 1994 for
% calcium currents with dynamic calcium concentrations, intracellular
% calcium buffering, and perineural calcium diffucsion into the
% extracellular bath

% This file has been modified to model C-Fiber Voltage Clamps at 37C and 
% with C-fiber conductances and shift values.

clear
close all

%% Constants
T=311;
gCat=.35; % nS
gCan=3.0; % nS
shift=-7;
dt=0.005;
R=8314;
F=96500;
ku=100;
kr=0.238;
Bi=0.001;
Voli=0.0127;
Vols=0.00146;
tca=4511.0;
%% Cat w/ buffering
V=-60:10:-20;
for j=1:length(V)
    Cait=0.000117;
    Caot=2.0;
    Cab=2.0;
    nb=4;
    Ecat(1)=R*T/(2*F)*log(Caot(1)/Cait(1))-78.7;
    dt0=1.0/(1.0+exp((-80+54.00+shift)/-5.75));
    ft0=1.0/(1.0+exp((-80+68.00+shift)/6.0));
    dtinf=1.0/(1.0+exp((V(j)+54.00+shift)/-5.75));
    taudt=22*exp(-(0.052)^2*(V(j)+68.0)^2)+2.5;
    ftinf=1.0/(1.0+exp((V(j)+68.00+shift)/6.0));
    tauft=103*exp(-(0.050)^2*(V(j)+58.0)^2)+12.5;
    Q10dt=1.90;
    Q10ft=2.20;
    taudt=taudt*Q10dt^((22-T+274)/10);
    tauft=tauft*Q10ft^((22-T+274)/10);

    t=0:0.005:200-0.005;

    dtime=t(2);
    dtnum(1)=dt0;
    ftnum(1)=ft0;
    Oc(1)=0.05;
    dOc(1)=0;
    ddtnum(1)=0;
    dftnum(1)=0;
    dCai(1)=0;
    dCa0(1)=0;

    for i=2:length(t)
        ddtnum(i)=(dtinf-dtnum(i-1))/taudt;
        dftnum(i)=(ftinf-ftnum(i-1))/tauft;
        dtnum(i)=dtnum(i-1)+ddtnum(i)*dtime;
        ftnum(i)=ftnum(i-1)+dftnum(i)*dtime;
        dOc(i)=ku*Cait(i-1)*(1-Oc(i-1))-kr*Oc(i-1);
        Oc(i)=Oc(i-1)+dOc(i)*dtime;
        Ecat(i)=R*T/(2*F)*log(Caot(i-1)/Cait(i-1))-78.7;
        Icatnum(j,i)=gCat.*dtnum(i).*ftnum(i).*(V(j)-Ecat(i));
        dCai(i)=-Icatnum(j,i)/10^3/(2*Voli*F)-nb*Bi*dOc(i);
        dCao(i)=(Caot(i-1)-Cab)/tca+(Icatnum(j,i)/10^3/(2*Vols*F));
        Cait(i)=Cait(i-1)+dCai(i)*dtime;
        Caot(i)=Caot(i-1)+dCao(i)*dtime;
    end

    
end

%% Cat no buffering
V=-60:10:-20;
for j=1:length(V)
    Caitnb=0.000117;
    Caotnb=2.0;
    Cab=2.0;
    nb=0;
    Ecat(1)=R*T/(2*F)*log(Caotnb(1)/Caitnb(1))-78.7;
    dt0=1.0/(1.0+exp((-80+54.00+shift)/-5.75));
    ft0=1.0/(1.0+exp((-80+68.00+shift)/6.0));
    dtinf=1.0/(1.0+exp((V(j)+54.00+shift)/-5.75));
    taudt=22*exp(-(0.052)^2*(V(j)+68.0)^2)+2.5;
    ftinf=1.0/(1.0+exp((V(j)+68.00+shift)/6.0));
    tauft=103*exp(-(0.050)^2*(V(j)+58.0)^2)+12.5;
    Q10dt=1.90;
    Q10ft=2.20;
    taudt=taudt*Q10dt^((22-T+274)/10);
    tauft=tauft*Q10ft^((22-T+274)/10);

    t=0:0.005:200-0.005;

    dtime=t(2);
    dtnum(1)=dt0;
    ftnum(1)=ft0;
    Oc(1)=0.05;
    dOc(1)=0;
    ddtnum(1)=0;
    dftnum(1)=0;
    dCai(1)=0;
    dCa0(1)=0;

    for i=2:length(t)
        ddtnum(i)=(dtinf-dtnum(i-1))/taudt;
        dftnum(i)=(ftinf-ftnum(i-1))/tauft;
        dtnum(i)=dtnum(i-1)+ddtnum(i)*dtime;
        ftnum(i)=ftnum(i-1)+dftnum(i)*dtime;
        dOc(i)=ku*Caitnb(i-1)*(1-Oc(i-1))-kr*Oc(i-1);
        Oc(i)=Oc(i-1)+dOc(i)*dtime;
        Ecat(i)=R*T/(2*F)*log(Caotnb(i-1)/Caitnb(i-1))-78.7;
        Icatnumnb(j,i)=gCat.*dtnum(i).*ftnum(i).*(V(j)-Ecat(i));
        dCai(i)=(-Icatnumnb(j,i))/10^3/(2*Voli*F)-nb*Bi*dOc(i);
        dCao(i)=(Caotnb(i-1)-Cab)/tca+(Icatnumnb(j,i)/10^3/(2*Vols*F));
        Caitnb(i)=Caitnb(i-1)+dCai(i)*dtime;
        Caotnb(i)=Caotnb(i-1)+dCao(i)*dtime;
    end 
end

%% Can buffering

V=-20:10:20;

for j=1:length(V)
    Cai=0.000117;
    Cao=2.0;
    nb=4;
    Eca=R*T/(2*F)*log(Cao/Cai)-78.7;
    dn0=1.0/(1.0+exp((-80+20.0+shift)/-4.5));
    fn10=1.0/(1.0+exp((-80+20.0+shift)/25.0));
    rn0=0.2/(1.0+exp((-80+5.0+shift)/-10.0));
    fn20=rn0+1.0/(1.0+exp((-80+40.0+shift)/10.0));
    dninf=1.0/(1.0+exp((V(j)+20.0+shift)/-4.5));
    taudn=3.25*exp(-(0.042)^2*(V(j)+31.0)^2)+0.395;
    fn1inf=1.0/(1.0+exp((V(j)+20.0+shift)/25.0));
    taufn1=33.5*exp(-(0.0395)^2*(V(j)+30.0)^2)+5.0;
    rn=0.2/(1.0+exp((V(j)+5.0+shift)/-10.0));
    fn2inf=rn+1.0/(1.0+exp((V(j)+40.0+shift)/10.0));
    taufn2=225.0*exp(-(0.0275)^2*(V(j)+40.0)^2)+75.0;
    Q10can=4.30;
    taudn=taudn*Q10can^((22-T+274)/10);
    taufn1=taufn1*Q10can^((22-T+274)/10);
    taufn2=taufn2*Q10can^((22-T+274)/10);

    t=0:0.005:200-0.005;



    dtime=t(2);
    ddnnum(1)=dn0;
    fn1num(1)=fn10;
    fn2num(1)=fn20;
    ddnnum(1)=0;
    dfn1num(1)=0;
    dfn2num(1)=0;
    Oc(1)=0.05;
    dOc(1)=0;
    dCai(1)=0;
    dCa0(1)=0;

    for i=2:length(t)
        ddnnum(i)=(dninf-ddnnum(i-1))/taudn;
        dfn1num(i)=(fn1inf-fn1num(i-1))/taufn1;
        dfn2num(i)=(fn2inf-fn2num(i-1))/taufn2;
        ddnnum(i)=ddnnum(i-1)+ddnnum(i)*dtime;
        fn1num(i)=fn1num(i-1)+dfn1num(i)*dtime;
        fn2num(i)=fn2num(i-1)+dfn2num(i)*dtime;
        dOc(i)=ku*Cai(i-1)*(1-Oc(i-1))-kr*Oc(i-1);
        Oc(i)=Oc(i-1)+dOc(i)*dtime;
        Eca(i)=R*T/(2*F)*log(Cao(i-1)/Cai(i-1))-78.7;
        Icannum(j,i)=gCan.*ddnnum(i).*(0.55*fn1num(i)+0.45*fn2num(i)).*(V(j)-Eca(i));
        dCai(i)=-Icannum(j,i)/10^3/(2*Voli*F)-nb*Bi*dOc(i);
        dCao(i)=(Cao(i-1)-Cab)/tca+(Icannum(j,i)/10^3/(2*Vols*F));
        Cai(i)=Cai(i-1)+dCai(i)*dtime;
        Cao(i)=Cao(i-1)+dCao(i)*dtime;
    end


end

%% Can no buffering

V=-20:10:20;

for j=1:length(V)
    Cainb=0.000117;
    Caonb=2.0;
    Eca=R*T/(2*F)*log(Caonb/Cainb)-78.7;
    nb=0;
    dn0=1.0/(1.0+exp((-80+20.0+shift)/-4.5));
    fn10=1.0/(1.0+exp((-80+20.0+shift)/25.0));
    rn0=0.2/(1.0+exp((-80+5.0+shift)/-10.0));
    fn20=rn0+1.0/(1.0+exp((-80+40.0+shift)/10.0));
    dninf=1.0/(1.0+exp((V(j)+20.0+shift)/-4.5));
    taudn=3.25*exp(-(0.042)^2*(V(j)+31.0)^2)+0.395;
    fn1inf=1.0/(1.0+exp((V(j)+20.0+shift)/25.0));
    taufn1=33.5*exp(-(0.0395)^2*(V(j)+30.0)^2)+5.0;
    rn=0.2/(1.0+exp((V(j)+5.0+shift)/-10.0));
    fn2inf=rn+1.0/(1.0+exp((V(j)+40.0+shift)/10.0));
    taufn2=225.0*exp(-(0.0275)^2*(V(j)+40.0)^2)+75.0;
    Q10can=4.30;
    taudn=taudn*Q10can^((22-T+274)/10);
    taufn1=taufn1*Q10can^((22-T+274)/10);
    taufn2=taufn2*Q10can^((22-T+274)/10);

    t=0:0.005:200-0.005;
    
    dtime=t(2);
    ddnnum(1)=dn0;
    fn1num(1)=fn10;
    fn2num(1)=fn20;
    ddnnum(1)=0;
    dfn1num(1)=0;
    dfn2num(1)=0;
    Oc(1)=0.05;
    dOc(1)=0;
    dCai(1)=0;
    dCa0(1)=0;

    for i=2:length(t)
        ddnnum(i)=(dninf-ddnnum(i-1))/taudn;
        dfn1num(i)=(fn1inf-fn1num(i-1))/taufn1;
        dfn2num(i)=(fn2inf-fn2num(i-1))/taufn2;
        ddnnum(i)=ddnnum(i-1)+ddnnum(i)*dtime;
        fn1num(i)=fn1num(i-1)+dfn1num(i)*dtime;
        fn2num(i)=fn2num(i-1)+dfn2num(i)*dtime;
        dOc(i)=ku*Cainb(i-1)*(1-Oc(i-1))-kr*Oc(i-1);
        Oc(i)=Oc(i-1)+dOc(i)*dtime;
        Eca(i)=R*T/(2*F)*log(Caonb(i-1)/Cainb(i-1))-78.7;
        Icannumnb(j,i)=gCan.*ddnnum(i).*(0.55*fn1num(i)+0.45*fn2num(i)).*(V(j)-Eca(i)); %pA
        dCai(i)=-Icannumnb(j,i)/10^3/(2*Voli*F)-nb*Bi*dOc(i);
        dCao(i)=((Caonb(i-1)-Cab)/tca)+(Icannumnb(j,i)/10^3/(2*Vols*F));
        Cainb(i)=Cainb(i-1)+dCai(i)*dtime;
        Caonb(i)=Caonb(i-1)+dCao(i)*dtime;
    end


end

%save('33CMatlabClamp/MatlabCa','Cai','Cainb','Cait','Caitnb','Icannum','Icannumnb','Icatnum','Icatnumnb','t')