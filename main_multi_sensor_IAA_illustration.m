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
s2=1*exp(2*pi*1i*(0.15*n+0.35*n.^3/(128*128*3)));
%s1=hamming(128)'.*exp(2*pi*1i*(0.05*n+0.35*n.^3/(128*128*3)));
%s2=1.*exp(2*pi*1i*(0.5*n-0.45*n.^3/(128*128*3)));
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
         
         %indexs=randi(128,N_sensors,64);
             for i=1:N_sensors
                indexs{i}=randperm(128,96);
                 X(i,indexs{i})=0;
                 iiii{i}=find(X(i,:)~=0);
             end
           
num=2;
NS=2;
%IF_O=2*IF_O/length(IF_O);
% HADTFD BASED

%for snr=-10:2:10
    %    Sig=awgn(SigO,30,'measured');
 %      Sig=awgn(SigO,30,'measured');
%p = randperm(length(Sig));

%p=p(1:44);
%  p=[11:32  65:85  100:111];
% 
% % p=[17:32  85:90  100:106];
 % p=[12:23  66:75  100:109];
% % 
% p=[12:30 40:55 66:75  100:109];
% 
% p=[10:30 41:54 66:80  95:109];
% 
% % p=[20:35 60:75 95:115];
% 
% Sig(p)=0;
% iii=find(Sig~=0);
% % ORIGINAL
% delta=5;
% alpha = 5;
% mean(abs(SigO-Sig))

[Y,~] = FAST_IF_Recover(X(1,:),63, num, 2,100,0,0,iiii{1});
[fidexmult,X_IT] = Multi_Sensor_FASTEST_IF_Recover(X,N_sensors,63, 2, 3,30,0,0,2,128,iiii);
mean(mean(abs(X-X_O)))
mean(mean(abs(X_O-X_IT)))
figure;
plot(real(X_O(2,:)));
hold on;
plot(real(X_IT(2,:)),'r:');



%xlabel('Samples')
%ylabel('Amplitude')
%legend('The Original Signal','Sparsely Sampled Signal','The Reconstructed Signal (Single Iteration)','The Reconstructed Signal (3 Iterations)');

  I1=HTFD_new1(X_IT(2,:),3,8,64);
  I2=HTFD_new1(X_O(2,:),3,8,64);
 I3=HTFD_new1(X(2,:),3,8,64);
%  I4=HTFD_new1(Sig_out,3,8,64);
% 
 figure;
 imagesc(I3)

figure; 
 imagesc(I2)
 figure;
 imagesc(I1)
% figure;
% imagesc(I4);
