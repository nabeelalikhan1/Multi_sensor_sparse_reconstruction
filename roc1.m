function [PF, PD] = roc1(Y,Z,Q,P)
FP = 0;
TP = 0;

% Sorting input Detection Vector
[Tx_sort,index] = sort(Y,'descend');
c_sort = Z(index);

Tx_prev = -inf;
PD = [];
PF = [];
i = 1;
while i <= length(Y) 
    
    if Tx_sort(i) ~= -inf
        
        PD(i) = TP/P;
        PF(i) = FP/Q;
        
%         R(:,i) = [PD(i) PF(i)];
        Tx_prev = Tx_sort(i);
    else
         
        PD(i) = TP/P;
        PF(i) = FP/Q;
        
%         R(:,i) = [PD(i) PF(i)];
        Tx_prev = Tx_sort(i);
    end
    
    if c_sort(i) == 1
        
        TP = TP + 1;
        
    else
        
        FP = FP + 1;
        
    end
    % Last induction of points
    PD1 = TP/P;
    PF1 = FP/Q;
    i = i + 1;
    
end
PD = [PD PD1];
PF = [PF PF1];