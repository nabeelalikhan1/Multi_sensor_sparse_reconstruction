function [extr_Sig,iniIFfit] = ICCD_sparse_multi(Sig,SampFreq,iniIFset,orderamp,alpha,iii,H)
% Intrinsic Chirp Component Decomposition(ICCD)
% the code is only for complex-valued data analysis
%%%%%%%%%%%%%%%%%%%%%%%  input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sig??measured signal,a row vector
% SampFreq: sampling frequency
% iniIFset: the instantaneous frequency series
% orderIF: the order of the Fourier model used for fitting the instantaneous frequency of the signal
% orderamp??Fourier order for characterizing signal amplitudes
% alpha??Tikhonov regularization parameter for ICCD.
% iii Samples that have been correctly sampled
%%%%%%%%%%%%%%%%%%%%%%%  output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extr_Sig: the reconstructed components
% ampmatrix: the estimated amplitudes
% iniIFfit??the fitted instantaneous frequencies (IFs)
% H is mixing matrix.
%if (isreal(Sig))
%Sig = hilbert(Sig);
%end
[NS,~]=size(Sig);
[multin,N] = size(iniIFset);%multin denotes the number of the components??N is the number of the samples
% when multin is larger than one, the algorithm performs the joint-optimization scheme.


dt = [0:N-1]/SampFreq; %time

phase = zeros(multin,N);
iniIFfit = zeros(multin,N);
for i = 1:multin
 %[temp1,temp2] = coef_ovefour(iniIFset(i,:),SampFreq,orderIF);  % fitting the input IFs with a Fourier model
 %iniIFfit(i,:) = temp1;%fitted IF
 %phase(i,:) = temp2; % fitted phase
 
 iniIFfit(i,:) = iniIFset(i,:);%fitted IF
 phase(i,:) = filter(1,[1 -1],iniIFset(i,:)); % fitted phase
 
 
end
%%%%%%%%%%%%%%%%%%%%Construction of the Fourier matrix%%%%%%%%%%%%%%%%%%%%%%%
f0 = SampFreq/2/N;%%Base Frequency
l_t = 2*orderamp + 1;%%the number of the Fourier coefficients
tmatrix = zeros(N,l_t);
tmatrix(:,1) = ones(N,1);
for j = 2:l_t
        tmatrix(:,j) = cos(2*pi*f0*(j-1)*dt); 
        if j >(l_t+1)/2
            tmatrix(:,j) = sin(2*pi*f0*(j-((l_t+1)/2))*dt);
        end
end

l_kernel = multin*l_t;%
KMAT=[];
SIGG=[];
for k=1:NS
    kmatrix = zeros(N,l_kernel);%

for i = 1:multin
   C = exp(sqrt(-1)*2*pi*phase(i,:));%
   kmatrix(:,((i-1)*l_t+1):i*l_t) = H(k,i)* spdiags(C.', 0, N, N)*tmatrix;
end

kmatrix1=kmatrix(iii{k},:);
    KMAT=[KMAT;kmatrix1];
Sigg=Sig(k,iii{k});
Sigg=Sigg.';
SIGG=[SIGG;Sigg ];
end

%Sig = Sig(:);
%% regularized least-squares solution


%Imatrix1 = speye(size(kmatrix1,2));%Identity matrix
%theta =(alpha*Imatrix1 + kmatrix1'*kmatrix1)\(kmatrix1'*Sig(iii));%coefficient vector

Imatrix1 = speye(size(KMAT,2));%Identity matrix
theta =( KMAT'*KMAT)\(KMAT'*SIGG);%coefficient vector
%% extract each component
extr_Sig = zeros(NS,multin,N);
for k=1:NS
    kmatrix = zeros(N,l_kernel);%

for i = 1:multin
   C = exp(sqrt(-1)*2*pi*phase(i,:));%
   kmatrix(:,((i-1)*l_t+1):i*l_t) = H(k,i)* spdiags(C.', 0, N, N)*tmatrix;
    extr_Sig(k,i,:) = kmatrix(:,(i-1)*l_t+1:i*l_t)*theta((i-1)*l_t+1:i*l_t);
end

end
end
    


