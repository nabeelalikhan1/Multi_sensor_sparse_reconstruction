clc
clear all;
close all
SampFreq = 256/2;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;


Sig1 = 1*exp(1i*(1*pi*(30*t.^3))+1i*(2*pi*(0*t))); %300t»òÕß150t
Sig2 = 1*exp(1i*(1*pi*(30*t.^3))+1i*(2*pi*(10*t))); %300t»òÕß150t

Sig =1*Sig1 +0.75*Sig2;
Sig=hamming(length(Sig)).'.*Sig;
SigO =Sig;
cccc=30*3;
IF_O(:,1)=cccc*t.^2/2;
IF_O(:,2)=-cccc*t.^2/2+100/2;
%IF_O(:,3)=cccc*t.^2/2+20/2;
%IF_O(:,4)=-cccc*t.^2/2+115/2;



%IF_O(:,3)=90*t.^2/2+15;
WN=64;
wind_step=32;

%Sig=Sig.*([1:128 128:-1:1]);
num=2;
NS=100;
IF_O=2*IF_O/length(IF_O);
% HADTFD BASED
iiii=0;
jjjj=0;
for snr=0:10:30
 jjjj=jjjj+1;
 iiii=0;
    for N_S=8:4:20 
         iiii=iiii+1;

    for k1=1:NS
        Sig=SigO;
        Sig=awgn(Sig,snr,'measured');
        p=[];
        for i=1:4
        pp = 32*(i-1)+ randperm(32-N_S-1,1);
        p1=pp:1:pp+N_S;
        p=[ p p1];
        end
        Sig(p)=0;
        [NA]=find(Sig~=0);
        for kkkkk=0:4
            
            % ORIGINAL
            delta=5;
            alpha = 5;
            if kkkkk==0   %ADTFD+VITERBI
                             [ext_sig,findex] = sparse_reconstruction_FAST_IF(Sig, num+2,11,0.4,p,5,61);
               
            elseif kkkkk==1 %the new algorithm
                %[ext_sig,findex] = FAST_IF_ICCD_Sparse(Sig,length(Sig)/2-1, num, 2,100,0,0,NA,5-2);
                [ext_sig,findex] = FAST_IF_Recover(Sig,length(Sig)/(2)-1, num+2, 2,100,0.4,0,NA);

              %  [~,findex] = FAST_IF_ICCD_Sparse(ext_sig,length(Sig)/2-1, num, 2,100,0,0,NA,5);
             %   [ext_sig,~,~] = ICCD_sparse(Sig,1,findex,1,5,0,NA);
               % ext_sig=sum(ext_sig);
               % ext_sig(NA)=Sig(NA);
            elseif kkkkk==2
               [ext_sig] = STFT_RECONSTRUCTION(Sig,WN,wind_step,NA);
            elseif kkkkk==3
                                [ext_sig,findex] = FAST_IF_Recover_iterative(Sig,length(Sig)/(2)-1, num+2, 2,100,0.4,0,NA);

            else
            y=GradRec(Sig,p);
            end
            if kkkkk==0
                                mse_FAST_IF_TF_FILTER(k1)=mean(abs(ext_sig-SigO));

            elseif kkkkk==1
                                mse_FAST_IF_IAA(k1)=mean(abs(ext_sig-SigO));

            elseif kkkkk==2
                mse_STFT(k1)=mean(abs(ext_sig-SigO));
            elseif kkkkk==3
                mse_FAST_IF_IAA_iterative(k1)=mean(abs(ext_sig-SigO));
                
            else
                mse_GD(k1)=mean(abs(y-SigO));
            end
            
            
            
        end
        
        
    end
    mse_IAA(jjjj,iiii)=mean(mse_FAST_IF_IAA);
    mse_IAA_it(jjjj,iiii)=mean(mse_FAST_IF_IAA_iterative);
    
    mse_TF(jjjj,iiii)=mean(mse_FAST_IF_TF_FILTER);
    mse_ST(jjjj,iiii)=mean(mse_STFT);
    mse_GDD(jjjj,iiii)=mean(mse_GD);
    
    end
    mse_IAA_it(jjjj,:)

 mse_IAA(jjjj,:)
  mse_ST(jjjj,:)
 mse_TF(jjjj,:)
end

figure;
plot(8:4:20,mse_IAA_it(1,:),'k','linewidth',3);
hold on;

plot(8:4:20,mse_IAA(1,:),'g','linewidth',3);
hold on;
plot(8:4:20,mse_ST(1,:),'r','linewidth',3);
hold on;
plot(8:4:20,mse_TF(1,:),'b','linewidth',3);
hold on;
plot(8:4:20,mse_GDD(1,:),'y','linewidth',3);

xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('a')
legend('IRNC (3-iterations)','IRNC (1-iteration)','Modified OMP','TF filtering','Gradient Descent');
set(gca,'fontsize', 36)

figure;
plot(8:4:20,mse_IAA_it(2,:),'k','linewidth',3);
hold on;

plot(8:4:20,mse_IAA(2,:),'g','linewidth',3);
hold on;
plot(8:4:20,mse_ST(2,:),'r','linewidth',3);
hold on;
plot(8:4:20,mse_TF(2,:),'b','linewidth',3);
hold on;
plot(8:4:20,mse_GDD(2,:),'y','linewidth',3);

xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('b')

legend('IRNC (3-iterations)','IRNC (1-iteration)','Modified OMP','TF filtering','Gradient Descent');
set(gca,'fontsize', 36)

figure;
plot(8:4:20,mse_IAA_it(3,:),'k','linewidth',3);
hold on;

plot(8:4:20,mse_IAA(3,:),'g','linewidth',3);
hold on;
plot(8:4:20,mse_ST(3,:),'r','linewidth',3);
hold on;
plot(8:4:20,mse_TF(3,:),'b','linewidth',3);
hold on;
plot(8:4:20,mse_GDD(3,:),'y','linewidth',3);

xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('c')
legend('IRNC (3-iterations)','IRNC (1-iteration)','Modified OMP','TF filtering','Gradient Descent');
set(gca,'fontsize', 36)




figure;
plot(8:4:20,mse_IAA_it(4,:),'k','linewidth',3);
hold on;

plot(8:4:20,mse_IAA(4,:),'g','linewidth',3);
hold on;
plot(8:4:20,mse_ST(4,:),'r','linewidth',3);
hold on;
plot(8:4:20,mse_TF(4,:),'b','linewidth',3);
hold on;
plot(8:4:20,mse_GDD(4,:),'y','linewidth',3);

xlabel('Number of missing samples in each gap')
ylabel('Mean absolute errror')
title('d')
legend('IRNC (3-iterations)','IRNC (1-iteration)','Modified OMP','TF filtering','Gradient Descent');

set(gca,'fontsize', 36)

% MSE FIXED AMPLITUDE
%
%
%
%
%
%>> mse_ICDD
%
%mse_ICDD =
%
    %1.1191    1.1168    1.1535    1.2485    1.3135    1.4049
    %0.3427    0.3444    0.3408    0.3543    0.3887    0.5537
    %0.1127    0.1183    0.1372    0.1493    0.1940    0.2143
    %0.0478    0.0528    0.0764    0.0817    0.1367    0.1755
%
%>> mse_TF
%
%mse_TF =
%
    %1.2394    1.2470    1.2562    1.2871    1.2855    1.2994
    %0.4069    0.4203    0.4406    0.4661    0.4907    0.5593
    %0.1656    0.1865    0.2168    0.2571    0.3102    0.3714
    %0.0966    0.1245    0.1600    0.1938    0.2743    0.3327
%
%>> mse_ST
%
%mse_ST =
%
    %1.1890    1.2427    1.2568    1.2818    1.3107    1.2964
    %0.4153    0.4706    0.5292    0.5857    0.6262    0.6767
    %0.1666    0.2268    0.2982    0.3734    0.4272    0.4920
    %0.0955    0.1556    0.2359    0.3163    0.3770    0.4475

