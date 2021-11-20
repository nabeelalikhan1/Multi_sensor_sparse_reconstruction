close all;
clear all;
fs=200;
N_S=2500;
addpath('D:\tfsa_5-5\windows\win64_bin');
P=0;
Q=0;

M=3;

%s=exp(2*pi*1i*(0.45*n-0*0.2*n.^2/(2*128)-0.5*n.^3/(128*128*3)));%+exp(2*pi*1i*(0.3*n-0*0.2*n.^2/(2*128)-0.3*n.^3/(128*128*3)));
%s=exp(2*pi*1i*(0.5*n-0.5*n.^2/(2*128)));%+0.5*exp(2*pi*1i*(0.25*n-0.2*n.^2/(2*128)));%+exp(2*pi*1i*(0.3*n-0*0.2*n.^2/(2*128)-0.3*n.^3/(128*128*3)));
T1 = zeros(1,N_S);
T2 = zeros(1,N_S);
T3 = zeros(1,N_S);
T4 = zeros(1,N_S);
Ruho=2;
clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
addpath('E:\tfsa_5-5\windows\win64_bin');

t = 0:1/SampFreq:1-1/SampFreq;

n=0:127;
s1=exp(2*pi*1i*(0.05*n+0.35*n.^3/(128*128*3)));
s2=1*exp(2*pi*1i*(0.5*n-0.45*n.^3/(128*128*3)));
%s2=1*exp(2*pi*1i*(0.15*n+0.35*n.^3/(128*128*3)));

perc=0.4;

s = [(s1.') (s2.') ];
theta = [-10,10]*pi/180;

n_sources=2;
N_sensors=2;
N_C=2;
s_orig=s;
A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  % mixing matrix A
X = A*s.';
X_O=X;
S_level=32;
%indexs=randi(128,N_sensors,64);

num=2;

N_S=2500;
%N_S=20;
llll=0;
Ruho=2;
SSS=0;
for S_L=32:16:96
 %   for S_L=0:32:96

    SSS=SSS+1;
    llll=0;
for SNR=-10:2:0
    llll=llll+1;
    P=0;
    Q=0;
    for i=1:N_S
        Xn=awgn(X,SNR,'measured');
        
        sorig = s;
        
       % A=ones(1,NS);
        X = A*s.';
        X1=X;
       
 
        ii = mod(i,2);
        if  ii > 0
            Mul = 1;
            c(i) = 1;
            P = P + 1;
        else
            Mul = 0;
            c(i) = 0;
            Q = Q + 1;
        end
        %----------------------------------------------------
        %----------------------------------------------------
        
        sigma = 10^(-SNR/20);
        StdN = sqrt(sigma);
        % generate noise
        
        NoisePower_Uncetainity   = 1/Ruho + (Ruho - (1/Ruho)).*rand;
        Noise_UNC = NoisePower_Uncetainity*ones(N_sensors,1);
        StdN_NU   = StdN.*NoisePower_Uncetainity;
        
        
        w = StdN_NU.*(randn(N_sensors,128) + 1j*(randn(N_sensors,128)))./sqrt(2); % noise
        X=Mul*X+w;
         for i4=1:N_sensors
            indexs{i4}=randperm(128,S_L);
            if S_L>0
            X(i4,indexs{i4})=0;
            end
            iiii{i4}=find(X(i4,:)~=0);
         end
        
         
      %    [f,sig_den] = Multi_Sensor_FASTEST_IF_Recover(X,N_sensors,63, 2, 3,30,0,0,2,128,iiii);
                    [f,sig_den] = Multi_Sensor_FASTEST_IF_Recover(X,N_sensors,63, 4, 3,30,0.15,0,2,128,iiii);
%[f,sig_den,~] = Multi_Sensor_FASTEST_IF_Recover_ICCD(X,N_sensors,63, n_sources, 3,30,0,0,2,128,iiii);

  AA=cov(sig_den.');
 
        
        %  T1(i)=abs(det(AA))/(abs(prod(diag(AA))));
        T1(i)=1-det(AA)/prod(diag(AA));
        %T1(i)=sum(abs(AA(:)))/sum(diag(abs(AA)));
        T1(i)=real(T1(i));
        AA1=cov(X.');
 
        T2(i)=1-det(AA1)/prod(diag(AA1));
                T2(i)=real(T2(i));

               
    end
    [PF1, PD1] = roc1(T1,c,Q,P);
    %trapz(PD1,PF1)
    
    [PF2, PD2] = roc1(T2,c,Q,P);
    
    for iii=1:length(PF1)
        if PF1(iii+1)>=0.05
            PDD1(SSS,llll)= PD1(iii);
            break;
        end
    end
    for iii=1:length(PF2)
        if PF2(iii+1)>=0.05
            PDD2(SSS,llll)= PD2(iii);
            break;
        end
    end
end
end
figure;
SNR=-10:2:0;
plot(SNR,PDD1(1,:),'o-');
hold on;
plot(SNR,PDD2(1,:),'rx-');

xlabel('Signal to noise ratio (dB)');
ylabel('Detection Probability');
legend('The Proposed method', 'Generalized Coherence');


figure;
plot(SNR,PDD1(2,:),'o-');
hold on;
plot(SNR,PDD2(2,:),'rx-');

xlabel('Signal to noise ratio (dB)');
ylabel('Detection Probability');
legend('The Proposed method', 'Generalized Coherence');



figure;
plot(SNR,PDD1(3,:),'o-');
hold on;
plot(SNR,PDD2(3,:),'rx-');

xlabel('Signal to noise ratio (dB)');
ylabel('Detection Probability');
legend('The Proposed method', 'Generalized Coherence');
% 
 figure;
 plot(SNR,PDD1(4,:),'o-');
 hold on;
 plot(SNR,PDD2(4,:),'rx-');
% 
 xlabel('Signal to noise ratio (dB)');
 ylabel('Detection Probability');
 legend('The Proposed method', 'Generalized Coherence');
% 
% 
 figure;
 plot(SNR,PDD1(5,:),'o-');
 hold on;
 plot(SNR,PDD2(5,:),'rx-');
% 
 xlabel('Signal to noise ratio (dB)');
 ylabel('Detection Probability');
 legend('The Proposed method', 'Generalized Coherence');
% 
% 
% figure;
% plot(SNR,PDD1(6,:),'o-');
% hold on;
% plot(SNR,PDD2(6,:),'rx-');
% 
% xlabel('Signal to noise ratio (dB)');
% ylabel('Detection Probability');
% legend('The Proposed method', 'Generalized Coherence');
