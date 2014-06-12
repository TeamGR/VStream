function Ypred = kNNClassify(Xtr, Ytr, k, Xte)
%
% function Y_pred = kNNClassify(Xtr, Ytr,  k, X)
%
% INPUT PARAMETERS
%   Xtr training input
%   Ytr training output 
%   k number of neighbours
%   Xte test input
% 
% OUTPUT PARAMETERS
%   Ypred estimated test output
%
% EXAMPLE
%   Ypred=kNNClassify(Xtr,Ytr,5,Xte);

    n = size(Xtr,1);
    m = size(Xte,1);
    
    
    if k > n
        k = n;
    end
    
    Ypred = zeros(m,1);

    DistMat = SquareDist(Xtr, Xte);
    
    
    
    for j= 1:m
        SortdDistMat = DistMat(:,j);
        
        [~, I] = sort(SortdDistMat);
        idx = I(1:k);
        val = sum(Ytr(idx))/k;

        Ypred(j) = sign(val);
    end
end

