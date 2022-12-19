% Author: David Catherall; Duke University
% Created: September 2017
% Description: Matlab Voltage Clamp Implementation for Schild 1994 for
% primary current, using constant calcium ion concentrations

clear
close all
%% Constants

% Schild 1994 Voltage Clamps Figure Values
% T=296;
% shiftca=0;
% shiftnaf=0;
% shiftnas=0;
% shiftk=0;
% gCat=2.75e-9;
% gNaf=29.5e-9;
% gNas=25.5e-9;
% gKd=30.0e-9;
% gKa=28.0e-9;
% gKds=11.5e-9;
% gCan=22.5e-9;

% C-Fiber Values
T=311;
shiftca=-7;
shiftnaf=-17.5;
shiftnas=-20.0;
shiftk=3.0;
gNaf=1.95e-6;
gNas=0.0295e-6;
gCat=0.00035e-6;
gCan=0.003e-6;
gKd=0.0051e-6;
gKa=0.004e-6;
gKds=0.003e-6;

R=8314;
F=96500;
Nai=8.9;
Nao=154;
Ena=R*T/(F)*log(Nao/Nai);
Ko=5.4;
Ki=145.0;
Ek=R*T/(F)*log(Ko/Ki);
Cai=0.000117;
Cao=2.0;
Eca=R*T/(2*F)*log(Cao/Cai)-78.7;

%% Cat
V=-60:10:-20;
%Analytical
for j=1:length(V)
    Q10dt=1.90;
    Q10ft=2.20;
    dt0=1.0/(1.0+exp((-80+54.00+shiftca)/-5.75));
    ft0=1.0/(1.0+exp((-80+68.00+shiftca)/6.0));
    dtinf=1.0/(1.0+exp((V(j)+54.00+shiftca)/-5.75));
    taudt=22*exp(-(0.052)^2*(V(j)+68.0)^2)+2.5;
    taudt=taudt*Q10dt^((22-T+274)/10);
    ftinf=1.0/(1.0+exp((V(j)+68.00+shiftca)/6.0));
    tauft=103*exp(-(0.050)^2*(V(j)+58.0)^2)+12.5;
    tauft=tauft*Q10ft^((22-T+274)/10);

    t=0:0.005:200-0.005;

    dt=(dt0-dtinf)*exp(-t/taudt)+dtinf;
    ft=(ft0-ftinf)*exp(-t/tauft)+ftinf;

    Icat(j,:)=gCat.*dt.*ft.*(V(j)-Eca)*10^9;

    dtime=t(2);
    dtnum(1)=dt0;
    ftnum(1)=ft0;
    ddtnum(1)=0;
    dftnum(1)=0;
% Numerical
    for i=2:length(t)
        ddtnum(i)=(dtinf-dtnum(i-1))/taudt;
        dftnum(i)=(ftinf-ftnum(i-1))/tauft;
        dtnum(i)=dtnum(i-1)+ddtnum(i)*dtime;
        ftnum(i)=ftnum(i-1)+dftnum(i)*dtime;
    end

    Icatnum(j,:)=gCat.*dtnum.*ftnum.*(V(j)-Eca)*10^9;
end

%% Naf

V=[-40:10:0];
% Analytical
for j=1:length(V)
    Q10m=2.30;
    Q10h=1.50;
    mf0=1.0/(1.0+exp((-80+41.35+shiftnaf)/-4.75));
    hf0=1.0/(1.0+exp((-80+62.00+shiftnaf)/4.5));
    j0=1.0/(1.0+exp((-80+40.00)/1.50));
    mfinf=1.0/(1.0+exp((V(j)+41.35+shiftnaf)/-4.75));
    taumf=0.75*exp(-(0.0635)^2*(V(j)+40.35)^2)+0.12;
    taumf=taumf*Q10m^((22-T+274)/10);
    hfinf=1.0/(1.0+exp((V(j)+62.00+shiftnaf)/4.5));
    tauhf=6.5*exp(-(0.0295)^2*(V(j)+75.00)^2)+0.55;
    tauhf=tauhf*Q10h^((22-T+274)/10);
    jinf=1.0/(1.0+exp((V(j)+40.00)/1.50));
    tauj=(25.0/(1.0+exp((V(j)-20.00)/4.50)))+0.01;

    t=0:0.005:10-0.005;

    mf=(mf0-mfinf)*exp(-t/taumf)+mfinf;
    hf=(hf0-hfinf)*exp(-t/tauhf)+hfinf;
    jf=(j0-jinf)*exp(-t/tauj)+jinf;

    Inaf(j,:)=gNaf.*mf.^3.*hf.*jf.*(V(j)-Ena)*10^9;

    dtime=t(2);
    mfnum(1)=mf0;
    hfnum(1)=hf0;
    jnum(1)=j0;
    dmfnum(1)=0;
    dhfnum(1)=0;
    djnum(1)=0;
% Numerical
    for i=2:length(t)
        dmfnum(i)=(mfinf-mfnum(i-1))/taumf;
        dhfnum(i)=(hfinf-hfnum(i-1))/tauhf;
        djnum(i)=(jinf-jnum(i-1))/tauj;
        mfnum(i)=mfnum(i-1)+dmfnum(i)*dtime;
        hfnum(i)=hfnum(i-1)+dhfnum(i)*dtime;
        jnum(i)=jnum(i-1)+djnum(i)*dtime;
    end

    Inafnum(j,:)=gNaf.*mfnum.^3.*hfnum.*jnum.*(V(j)-Ena)*10^9;
end

% Nas

V=-10:10:30;
% Analytical
for j=1:length(V)
    ms0=1.0/(1.0+exp((-80+20.35+shiftnas)/-4.45));
    hs0=1.0/(1.0+exp((-80+18.00+shiftnas)/4.5));
    msinf=1.0/(1.0+exp((V(j)+20.35+shiftnas)/-4.45));
    taums=1.50*exp(-(0.0595)^2*(V(j)+20.35)^2)+0.15;
    taums=taums*Q10m^((22-T+274)/10);
    hsinf=1.0/(1.0+exp((V(j)+18.00+shiftnas)/4.5));
    tauhs=4.95*exp(-(0.0335)^2*(V(j)+20.00)^2)+0.75;
    tauhs=tauhs*Q10h^((22-T+274)/10);

    t=0:0.005:10-0.005;

    ms=(ms0-msinf)*exp(-t/taums)+msinf;
    hs=(hs0-hsinf)*exp(-t/tauhs)+hsinf;

    Inas(j,:)=gNas.*ms.^3.*hs.*(V(j)-Ena)*10^9;

    dtime=t(2);
    msnum(1)=ms0;
    hsnum(1)=hs0;
    dmsnum(1)=0;
    dhsnum(1)=0;
% Numerical
    for i=2:length(t)
        dmsnum(i)=(msinf-msnum(i-1))/taums;
        dhsnum(i)=(hsinf-hsnum(i-1))/tauhs;
        msnum(i)=msnum(i-1)+dmsnum(i)*dtime;
        hsnum(i)=hsnum(i-1)+dhsnum(i)*dtime;
    end

    Inasnum(j,:)=gNas.*msnum.^3.*hsnum.*(V(j)-Ena)*10^9;
end

% Kd

V=-40:10:40;
% Analytical
for j=1:length(V)
    Q10n=1.40;
    n0=1.0/(1.0+exp((-80+14.62+shiftk)/-18.38));
    ninf=1.0/(1.0+exp((V(j)+14.62+shiftk)/-18.38));
    an=(0.001265*(V(j)+14.273))/(1.0-exp((V(j)+14.273)/-10.0));
    bn=0.125*exp((V(j)+55.0)/-2.5);
    taun=(1/(an+bn))+1;
    taun=taun*Q10n^((22-T+274)/10);

    t=0:0.005:350-0.005;

    n=(n0-ninf)*exp(-t/taun)+ninf;

    Ikd(j,:)=gKd.*n.*(V(j)-Ek)*10^9;

    dtime=t(2);
    nnum(1)=n0;
    dnnum(1)=0;
% Numerical
    for i=2:length(t)
        dnnum(i)=(ninf-nnum(i-1))/taun;
        nnum(i)=nnum(i-1)+dnnum(i)*dtime;
    end

    Ikdnum(j,:)=gKd.*nnum.*(V(j)-Ek)*10^9;
end

% Ktrans

V=-40:10:30;
% Analytical 
for j=1:length(V)
    Q10kt=1.93;
    p0=1.0/(1.0+exp((-80+28.0+shiftk)/-28.0));
    q0=1.0/(1.0+exp((-80+58.0+shiftk)/7.0));
    pinf=1.0/(1.0+exp((V(j)+28.0+shiftk)/-28.0));
    taup=5.0*exp(-(0.022)^2*(V(j)+65.0)^2)+2.5;
    taup=taup*Q10kt^((22-T+274)/10);
    qinf=1.0/(1.0+exp((V(j)+58.0+shiftk)/7.0));
    tauq=100.0*exp(-(0.035)^2*(V(j)+30.00)^2)+10.5;
    tauq=tauq*Q10kt^((22-T+274)/10);

    t=0:0.005:600-0.005;

    p=(p0-pinf)*exp(-t/taup)+pinf;
    q=(q0-qinf)*exp(-t/tauq)+qinf;

    Ika(j,:)=gKa.*p.^3.*q.*(V(j)-Ek)*10^9;

    dtime=t(2);
    pnum(1)=p0;
    qnum(1)=q0;
    dpnum(1)=0;
    dqnum(1)=0;
% Numerical
    for i=2:length(t)
        dpnum(i)=(pinf-pnum(i-1))/taup;
        dqnum(i)=(qinf-qnum(i-1))/tauq;
        pnum(i)=pnum(i-1)+dpnum(i)*dtime;
        qnum(i)=qnum(i-1)+dqnum(i)*dtime;
    end

    Ikanum(j,:)=gKa.*pnum.^3.*qnum.*(V(j)-Ek)*10^9;
end
% Analytical
for j=1:length(V)
    x0=1.0/(1.0+exp((-80+39.59+shiftk)/-14.68));
    y0=1.0/(1.0+exp((-80+48.0+shiftk)/7.0));
    xinf=1.0/(1.0+exp((V(j)+39.59+shiftk)/-14.68));
    taux=5.0*exp(-(0.022)^2*(V(j)+65.0)^2)+2.5;
    taux=taux*Q10kt^((22-T+274)/10);
    yinf=1.0/(1.0+exp((V(j)+48.0+shiftk)/7.0));
    tauy=7500.0;
    tauy=tauy*Q10kt^((22-T+274)/10);

    t=0:0.005:600-0.005;

    x=(x0-xinf)*exp(-t/taux)+xinf;
    y=(y0-yinf)*exp(-t/tauy)+yinf;

    Ikds(j,:)=gKds.*x.^3.*y.*(V(j)-Ek)*10^9;

    dtime=t(2);
    xnum(1)=x0;
    ynum(1)=y0;
    dxnum(1)=0;
    dynum(1)=0;
% Numerical
    for i=2:length(t)
        dxnum(i)=(xinf-xnum(i-1))/taux;
        dynum(i)=(yinf-ynum(i-1))/tauy;
        xnum(i)=xnum(i-1)+dxnum(i)*dtime;
        ynum(i)=ynum(i-1)+dynum(i)*dtime;
    end

    Ikdsnum(j,:)=gKds.*xnum.^3.*ynum.*(V(j)-Ek)*10^9;
end

Iktrans=Ika+Ikds;
Iktransnum=Ikanum+Ikdsnum;

% Can

V=-20:10:20;
% Analytical
for j=1:length(V)
    Q10can=4.30;
    dn0=1.0/(1.0+exp((-80+20.0+shiftca)/-4.5));
    fn10=1.0/(1.0+exp((-80+20.0+shiftca)/25.0));
    rn0=0.2/(1.0+exp((-80+5.0+shiftca)/-10.0));
    fn20=rn0+1.0/(1.0+exp((-80+40.0+shiftca)/10.0));
    dninf=1.0/(1.0+exp((V(j)+20.0+shiftca)/-4.5));
    taudn=3.25*exp(-(0.042)^2*(V(j)+31.0)^2)+0.395;
    taudn=taudn*Q10can^((22-T+274)/10);
    fn1inf=1.0/(1.0+exp((V(j)+20.0+shiftca)/25.0));
    taufn1=33.5*exp(-(0.0395)^2*(V(j)+30.0)^2)+5.0;
    taufn1=taufn1*Q10can^((22-T+274)/10);
    rn=0.2/(1.0+exp((V(j)+5.0+shiftca)/-10.0));
    fn2inf=rn+1.0/(1.0+exp((V(j)+40.0+shiftca)/10.0));
    taufn2=225.0*exp(-(0.0275)^2*(V(j)+40.0)^2)+75.0;
    taufn2=taufn2*Q10can^((22-T+274)/10);

    t=0:0.005:200-0.005;

    dn=(dn0-dninf)*exp(-t/taudn)+dninf;
    fn1=(fn10-fn1inf)*exp(-t/taufn1)+fn1inf;
    fn2=(fn20-fn2inf)*exp(-t/taufn2)+fn2inf;

    Ican(j,:)=gCan.*dn.*(0.55*fn1+0.45*fn2).*(V(j)-Eca)*10^9;

    dtime=t(2);
    ddnnum(1)=dn0;
    fn1num(1)=fn10;
    fn2num(1)=fn20;
    ddnnum(1)=0;
    dfn1num(1)=0;
    dfn2num(1)=0;
% Numerical
    for i=2:length(t)
        ddnnum(i)=(dninf-ddnnum(i-1))/taudn;
        dfn1num(i)=(fn1inf-fn1num(i-1))/taufn1;
        dfn2num(i)=(fn2inf-fn2num(i-1))/taufn2;
        ddnnum(i)=ddnnum(i-1)+ddnnum(i)*dtime;
        fn1num(i)=fn1num(i-1)+dfn1num(i)*dtime;
        fn2num(i)=fn2num(i-1)+dfn2num(i)*dtime;
    end

    Icannum(j,:)=gCan.*ddnnum.*(0.55*fn1num+0.45*fn2num).*(V(j)-Eca)*10^9;
end

% save('33CMatlabClamp/ConstClamps','Ican','Icat','Ikd','Iktrans','Inaf','Inas')