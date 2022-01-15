clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
n=0:127;
s1=hamming(128)'.*exp(2*pi*1i*(0.05*n+0.4*n.^3/(128*128*3)));
s2=hanning(128)'.*exp(2*pi*1i*(0.12*n+0.4*n.^3/(128*128*3)));
perc=0.4;

%IF_O(1,:)=0.05+0.45*3*n.^2/(128*128*3);
%IF_O(2,:)=0.5-0.45*3*n.^2/(128*128*3);
%IF_O=IF_O.';
s = [(s1.') (s2.') ];
theta = [-10,10]*pi/180;

n_sources=2;
N_sensors=2;
N_C=2;
s_orig=s;
A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  % mixing matrix A
X = A*s.';
X_O=X;

%IF_O(:,3)=cccc*t.^2/2+20/2;
%IF_O(:,4)=-cccc*t.^2/2+115/2;



%IF_O(:,3)=90*t.^2/2+15;
WN=64;
wind_step=32;

%Sig=Sig.*([1:128 128:-1:1]);
num=2;
NS=2;
% HADTFD BASED
iiii=0;
jjjj=0;
N_S=8;
X=X_O;

for ii=1:N_sensors
    p=[];
    for i=1:4
        pp = 32*(i-1)+ randperm(32-N_S-1,1);
        p1=pp:1:pp+N_S;
        p=[ p p1];
    end
    X(ii,p)=0;
    AV{ii}=find(X(ii,:)~=0);
    NAV{ii}=p;
end

[~,X_IT_new,X_IT_old] = Multi_Sensor_FASTEST_IF_Recover_ICCD(X,N_sensors,63, n_sources, 3,30,0,0,1,128,AV);
figure;
I2=HTFD_new1(X_O(2,:),3,8,64);
imagesc(I2);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('a','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

I1=HTFD_new1(X(2,:),3,8,64);
figure;
imagesc(I1)
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('b','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

I=HTFD_new1(X_IT_new(2,:),3,8,64);
figure;
imagesc(I)
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
title('c','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',20);

