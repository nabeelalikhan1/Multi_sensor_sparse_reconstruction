clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
n=0:127;
s1=hamming(128)'.*exp(2*pi*1i*(0.05*n+0.4*n.^3/(128*128*3)));
s2=hanning(128)'.*exp(2*pi*1i*(0.12*n+0.4*n.^3/(128*128*3)));
perc=0.4;
%s1=hamming(128)'.*exp(2*pi*1i*(0.05*n+0*0.4*n.^3/(128*128*3)));
%s2=hanning(128)'.*exp(2*pi*1i*(0.12*n+0*0.4*n.^3/(128*128*3)));
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
NS=100;
% HADTFD BASED
iiii=0;
jjjj=0;
for snr=40:10:40
    jjjj=jjjj+1;
    iiii=0;
    for N_S=4:4:16
        iiii=iiii+1;
        clear AV;
        for k1=1:NS
            X=X_O;
            if snr~=40
            X=awgn(X,snr,'measured');
            end
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
            
            for kkkkk=0:3
                
                % ORIGINAL
                delta=5;
                alpha = 5;
                if kkkkk==0   %ICCD 
                    [~,X_IT_new,X_IT_old] = Multi_Sensor_FASTEST_IF_Recover_ICCD(X,N_sensors,63, n_sources, 3,30,0,0,1,128,AV);
                    MSEE=[0 0];
                    for k22=1:2
                        if k22==1
                            XX=X_IT_new;
                        else
                            XX=X_IT_old;
                        end
                        for i22=1:N_sensors
                            for j22=1:N_sensors
                                MSE(j22)=mean(abs(X_O(i22,:)-XX(j22,:)));
                            end
                            [v,ind]=min(MSE);
                            MSEE(k22)=MSEE(k22)+v;
                            XX(ind,:)=1000;
                        end
                    end
                    MSEE=MSEE/N_sensors;
                    mse_JOINT_ICCD(k1)=MSEE(1);
                    mse_IND_ICCD(k1)=MSEE(2);
                    
                elseif kkkkk==1 %the new algorithm
                    MS1E=0;
                    for in=1:N_sensors
                        [ext_sig,findex] = sparse_reconstruction_FAST_IF(X(in,:), n_sources,3,NAV{in},5,61);
                        MS1E=MS1E+mean(abs(ext_sig-X_O(in,:)));
                    end
                    mse_TF(k1)=MS1E/N_sensors;
                elseif kkkkk==2
                    MS1E=0;
                    for in=1:N_sensors
                        [ext_sig] = STFT_RECONSTRUCTION(X(in,:),WN,wind_step,AV{in});
                        MS1E=MS1E+mean(abs(ext_sig-X_O(in,:)));
                    end
                    mse_STFT(k1)=MS1E/N_sensors;
                    
                elseif kkkkk==3
                    MS1E=0;
                    for in=1:N_sensors
                        [ext_sig,findex] = GradRec(X(in,:),NAV{in});
                        MS1E=MS1E+mean(abs(ext_sig-X_O(in,:)));
                    end
                    mse_GD(k1)=MS1E/N_sensors;
                
                else
                    
                    
                end
                
                
                
                
            end
            
            
        end
        
        mse_J_ICCD(jjjj,iiii)=mean(mse_JOINT_ICCD);
        mse_ICCD(jjjj,iiii)=mean(mse_IND_ICCD);
        mse_ST(jjjj,iiii)=mean(mse_STFT);
        mse_TFF(jjjj,iiii)=mean(mse_TF);
        mse_G(jjjj,iiii)=mean(mse_GD);
        
    end
    mse_J_ICCD(jjjj,:)
    mse_ICCD(jjjj,:)
    mse_ST(jjjj,:)
    
    mse_TFF(jjjj,:)
end
mse_J_ICCD=20*log10(mse_J_ICCD);
mse_ICCD=20*log10(mse_ICCD);
mse_ST=20*log10(mse_ST);
mse_TFF=20*log10(mse_TFF);
mse_G=20*log10(mse_G);


figure;

plot(4:4:16,mse_J_ICCD(1,:),'k','linewidth',3);
hold on;
plot(4:4:16,mse_ICCD(1,:),'r','linewidth',3);
hold on;
plot(4:4:16,mse_ST(1,:),'b','linewidth',3);
hold on;
plot(4:4:16,mse_TFF(1,:),'g','linewidth',3);
hold on;

plot(4:4:16,mse_G(1,:),'g','linewidth',3);

xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('e')
legend('Joint Reconstruction','Channel wise reconstruction','Modified OMP','TF filtering','Gradient Descent');

% MSE FIXED AMPLITUDE
%
%
%
%
%
