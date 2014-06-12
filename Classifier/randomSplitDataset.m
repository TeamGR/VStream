function [Xtr, Ytr, Xts, Yts] = randomSplitDataset(X, Y, n_train, n_test)
    n = size(X,1);
    I = randperm(n);
    I1 = I(1:n_train);
    I2 = I((n_train+1):(n_train+n_test));
    
    Xtr = X(I1,:);
    Ytr = Y(I1);
    
    Xts = X(I2,:);
    Yts = Y(I2);
end