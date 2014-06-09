% Simple classifier for 1-L signatures, n --> 1
% Author: Raffaello Camoriano

load signature_gabor_car.mat
sig_car = signatures;
clear signatures;

load signature_gabor_plane.mat
sig_pln = signatures;
clear signatures;

X = [ sig_car ; sig_pln ];
Y = [ ones(size(sig_car,1),1) ; -ones(size(sig_pln,1),1) ];

% Randomly split the dataset between training and testing

n_train = 2;
n_test = 48;
[Xtr, Ytr, Xts, Yts] = randomSplitDataset(X, Y, n_train, n_test);

% Apply 1-NN classification
k = 1;
Ypred = kNNClassify(Xtr, Ytr, k, Xts);
ind = find((sign(Ypred) ~= sign(Yts)));
[Yts  Ypred ]
missNum = numel(ind)
