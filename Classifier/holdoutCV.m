function [l, s, Vm, Vs, Tm, Ts] = holdoutCV(algorithm, X, Y, kernel, perc, nrip, intRegPar, intKerPar)
%[l, s, Vm, Vs, Tm, Ts] = holdoutCV(algorithm, X, Y, kernel, perc, nrip, intRegPar, intKerPar)
% algorithm =   'knn' k-Nearest neighbors,
%               'rls' regularized least squares,
%               'krls' kernel regularized least squares,
%               'lr' logistic regression,
%               'klr' kernel logistic regression, 
% X: the training examples
% Y: the training labels
% kernel: the kernel function (see help KernelMatrix).
% perc: fraction of the dataset to be used for validation
% nrip: number of repetitions of the test for each couple of parameters
% intRegPar: list of regularization parameters 
% intKerPar: list of kernel parameters 
%
% Output:
% l, s: the couple of regularization and kernel parameter that minimize the median of the
%       validation error
% Vm, Vs: median and variance of the validation error for each couple of parameters
% Tm, Ts: median and variance of the error computed on the training set for each couple
%       of parameters
%
% Example of code for kNN:
% intRegPar = 1:-0.1:0.1;
% [Xtr, Ytr] = MixGauss([[0;0],[1;1]],[0.5,0.25],100);
% [l, s, Vm, Vs, Tm, Ts] = holdoutCV('knn', Xtr, Ytr,[], 0.5, 5, intRegPar, []);
%
% Example of code for Linear LR:
% intRegPar = 1:-0.1:0.1;
% [Xtr, Ytr] = MixGauss([[0;0],[1;1]],[0.5,0.25],100);
% [l, s, Vm, Vs, Tm, Ts] = holdoutCV('lr', Xtr, Ytr,[], 0.5, 5, intRegPar, []);
%
% Example of code for Linear Regularized Least Squares:
% intRegPar = 1:-0.1:0.1;
% [Xtr, Ytr] = MixGauss([[0;0],[1;1]],[0.5,0.25],100);
% [l, s, Vm, Vs, Tm, Ts] = holdoutCV('rls', Xtr, Ytr,[], 0.5, 5, intRegPar, []);
%
% Example of code for Kernel LR:
% intRegPar = 1:-0.1:0.1;
% intKerPar = 10:-0.5:1;
% [Xtr, Ytr] = MixGauss([[0;0],[1;1]],[0.5,0.25],100);
% [l, s, Vm, Vs, Tm, Ts] = holdoutCV('klr', Xtr, Ytr,'gaussian', 0.5, 5, intRegPar, intKerPar);
%
% Example of code for Kernel RLS:
% intRegPar = 1:-0.1:0.1;
% intKerPar = 10:-0.5:1;
% [Xtr, Ytr] = MixGauss([[0;0],[1;1]],[0.5,0.25],100);
% [l, s, Vm, Vs, Tm, Ts] = holdoutCV('krls', Xtr, Ytr,'gaussian', 0.5, 5, intRegPar, intKerPar);
    if isequal(algorithm,'knn') || isequal(algorithm,'rls') || isequal(algorithm,'lr') || isequal(kernel,'linear') || isequal(kernel,[])
        intKerPar = [0];
    end

    nKerPar = numel(intKerPar);
    nRegPar = numel(intRegPar);
 
    
    n = size(X,1);
    ntr = ceil(n*(1-perc));
    
    tmn = zeros(nRegPar, nKerPar, nrip);
    vmn = zeros(nRegPar, nKerPar, nrip);
    
    for rip = 1:nrip
        I = randperm(n);
        Xtr = X(I(1:ntr),:);
        Ytr = Y(I(1:ntr),:);
        Xvl = X(I(ntr+1:end),:);
        Yvl = Y(I(ntr+1:end),:);
        
        
        is = 0;
        for s=intKerPar
            is = is + 1;
            il = 0;  
            
            if isequal(algorithm,'klr');
                c = zeros(ntr,1);
            elseif isequal(algorithm,'lr');
                c = zeros(size(Xtr,2),1);
            end
            
            for l=intRegPar
                il = il + 1;
         
                if isequal(algorithm, 'rls')
                    w = regularizedLSTrain(Xtr, Ytr, l);
                    tmn(il, is, rip) =  calcErr(regularizedLSTest(w, Xtr), Ytr);               
                    vmn(il, is, rip)  = calcErr(regularizedLSTest(w, Xvl), Yvl);
                elseif isequal(algorithm, 'krls')
                    w = regularizedKernLSTrain(Xtr, Ytr, kernel, s, l);
                    tmn(il, is, rip) =  calcErr(regularizedKernLSTest(w, Xtr, kernel, s, Xtr), Ytr);               
                    vmn(il, is, rip)  = calcErr(regularizedKernLSTest(w, Xtr, kernel, s, Xvl), Yvl);                
                elseif isequal(algorithm, 'klr')
                    c = kernLRTrain1(Xtr, Ytr, kernel, s, l, c);
                    tmn(il, is, rip) =  calcErrLR(kernLRTest(c, Xtr, kernel, s, Xtr), Ytr);               
                    vmn(il, is, rip)  = calcErrLR(kernLRTest(c, Xtr, kernel, s, Xvl), Yvl);
                elseif isequal(algorithm, 'lr')
                    c = linearLRTrain1(Xtr, Ytr, l, c);
                    tmn(il, is, rip) =  calcErrLR(linearLRTest(c, Xtr), Ytr);               
                    vmn(il, is, rip)  = calcErrLR(linearLRTest(c, Xvl), Yvl);   
                else  
                    tmn(il, is, rip) =  calcErr(kNNClassify(Xtr, Ytr, l, Xtr),Ytr);
                    vmn(il, is, rip)  = calcErr(kNNClassify(Xtr, Ytr, l, Xvl),Yvl);
                end
                
                str = sprintf('rip\tRegPar\tKerPar\tvalErr\ttrErr\n%d\t%f\t%f\t%f\t%f\n',rip, l, s, vmn(il,is,rip), tmn(il,is,rip));
                disp(str);
            end
        end
    end
    
    Tm = median(tmn,3);
    Ts = std(tmn,0,3);
    Vm = median(vmn,3);
    Vs = std(vmn,0,3);
    
    [row, col] = find(Vm <= min(min(Vm)));
    
    l = intRegPar(row(1));
    s = intKerPar(col(1));
end

function err = calcErrLR(T, Y)
    err = mean(log(1+exp(-Y.*T)));
end

function err = calcErr(T, Y)
    err = mean(sign(T)~=sign(Y));
end

function c = kernLRTrain1(Xtr, Ytr, kernel, kerpar, lambda, c)
    iter = 10000;
    epsilon = min(1e-6, 0.1*lambda);
    [n,  ~] = size(Xtr);
    K = KernelMatrix(Xtr, Xtr, kernel, kerpar);
    L = eigs(K, 1)/n + lambda;
    gamma = 1/L;
    j = 0;
    fold = 0;
    f = inf;
    while(j < iter && abs(f - fold) >= epsilon)
        fold = f;
        j = j + 1;
        p = exp(Ytr.*(K*c));
        c = c - gamma*(-(1/n)*(Ytr./(1+p)) + 2*lambda*c);
        f = sum(log(1 + exp(-Ytr.*(K*c))))/n + lambda*c'*c;
        %disp(['iter:', num2str(j),'  err:', num2str(abs(f - fold))]);
    end
end


function w = linearLRTrain1(Xtr, Ytr, lambda, c)    
    epsilon = min(1e-6, 0.1*lambda);
    iter = 10000;
    [n, ~] = size(Xtr);
    w = c;
    L = eigs(Xtr*Xtr', 1)/n + lambda;
    gamma = 1/L;
    
    j=0;
    fold = 0;
    f = inf;
    while j < iter && abs(f - fold) >= epsilon
        fold = f;
        j = j + 1;
        p = exp(-Ytr.* (Xtr*w));
        w = w - gamma*(-(1/n)*Xtr'*(Ytr.*p./(1+p)) + 2*lambda*w);
        f = sum(log(1 + exp(-Ytr.*(Xtr*w))))/n + lambda*w'*w;
        %disp(['iter:', num2str(j),'  err:', num2str(abs(f - fold))]);
    end
end
