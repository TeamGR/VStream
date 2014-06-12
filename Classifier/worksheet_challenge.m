% MNIST lab warm-up
% Author: raffaello.camoriano@iit.it

%% Loading
clear all;
close all;

% Load dataset
load('Challenge_Train.mat');

%% Model selection

numRep = 20;
testErrVecOverall = [];
trainSplit = 0.8;
n = numel(Ytr);
n_train = floor(trainSplit * n);
n_val = n - n_train;
intRegPar = [ 0.00001 0.0001 0.001 0.01 0.1 1 10 100 1000 10000 ];
intKerPar = [ 0.00001 0.0001 0.001 0.01 0.1 1 10 100 1000 10000 ];
intKpar = 1:50;
knnValErrMedian = [];
kernel = 'gaussian';
algo = {'kNN' ; 'RLS' ; 'KRLS'};

for j = 1:size(algo,1)
    valErrVec = [];

    switch algo{j}
        
       case 'kNN'
           % HP selection
           
           for k = intKpar
               knnValErrVec = [];
               for i  = 1 : numRep
                   [Xtr1, Ytr1, Xvs1, Yvs1] = randomSplitDataset(Xtr, Ytr, n_train, n_val);
                   Ypred = kNNClassify(Xtr1, Ytr1, k, Xvs1);

                   % Misclassifications
                   ind = find((sign(Ypred)~=sign(Yvs1)));
                   missNum = numel(ind);
                   knnValErrVec = [knnValErrVec ; missNum/numel(Yvs1)];
               end
               knnValErrMedian = [knnValErrMedian ; median(knnValErrVec)];
           end
           
           [Vm_kNN_best, I] = min(knnValErrMedian);
           k_best = intKpar(I);
           
       case 'RLS'
           % HP selection
           [l_RLS, s_RLS, Vm_RLS, Vs_RLS, Tm_RLS, Ts_RLS] = holdoutCV('rls', Xtr, Ytr, kernel, 1 - trainSplit, numRep, intRegPar, intKerPar);
           Vm_RLS_best = min(Vm_RLS);

       case 'KRLS'
          % HP selection
           [l_KRLS, s_KRLS, Vm_KRLS, Vs_KRLS, Tm_KRLS, Ts_KRLS] = holdoutCV('krls', Xtr, Ytr, kernel, 1 - trainSplit, numRep, intRegPar, intKerPar);
           Vm_KRLS_best = min(min(Vm_KRLS));
    end
end