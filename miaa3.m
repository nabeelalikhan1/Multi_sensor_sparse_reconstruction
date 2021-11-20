function ym = miaa3(data, dataIX, K, iterMax, A_m)
%Input:
%   data: column vector
%   dataIX: the index of data,row or column vector
%   K: Number of sampling frequency
%   iterMax:Maximum number of iteration
%Output:
%   alpha:estimated amplitude spectrum, column vector
%   alpha_0:initialization
%   difference:norm-2 difference between each iteration

    L = length(data);
    A = exp(2j*pi*dataIX(:)*(1:K)/K);

    alpha  = A'*data/L;
    energy = conj(alpha(:,1)).*alpha(:,1);
    P = diag(energy);

    for iterNum = 1:iterMax-1
        R = A*P*A';
        R_inv = inv(R);
        for k = 1:K
            alpha(k)= A(:,k)'*R_inv*data/(A(:,k)'*R_inv*A(:,k));
        end
        energy = conj(alpha).*alpha;
        energy(energy<10^-13)=10^-13;
        P=diag(energy);
    end

    coef = [];
    R = A*P*A';
    R_inv = inv(R);
    for k = 1:K
        temp = A(:,k)'*R_inv*data;
        alpha(k)= temp/(A(:,k)'*R_inv*A(:,k));
        coef = [coef temp];
    end
    energy = conj(alpha).*alpha;
                   
    ym = A_m*(energy.*coef(:));
%    dataIAA=alpha(:,iterNum);
    % fprintf('iterNum=%d\n',iterNum)
