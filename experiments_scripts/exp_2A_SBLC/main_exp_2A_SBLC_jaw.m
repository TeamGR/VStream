%% Experiment description
%
% Number of layers of the net: 2
% Employed dataset: SBLC
%
% Transformations:
%   Independent scaling, translation and rotation on a uniform
%   background
%
% Binary classification, 1-35.jpg (car) vs 2696-2730.jpg (plane)

close all;
clear all;

%% Load images and infer parameters


% Uniform BG images
n_uniform = 35;
Seq = cell(n_uniform,2);
imgPath = 'project\SBLC\images';
imgType = '*.jpg'; % change based on image type
images  = dir([imgPath '\' imgType]);
for jdx = 1:2
for idx = 1:n_uniform
    if jdx == 1
        Seq{idx, jdx} = double(rgb2gray(imread([imgPath '\' images(idx).name])));
    else
        Seq{idx, jdx} = double(rgb2gray(imread([imgPath '\' images(2695+idx).name])));
    end
end    
end
uniform_images = Seq;


n_images = 70;

%% Load templates of layer 1

%T = load('gabor_filters.mat');
T = load('pascal_filters_5.mat');

gabors = T.templates;

%% Load templates of layer 2

%T = load('templatesL2.mat');
%T = load('templatesL2_gabor.mat');
T = load('templatesL2-5_L1-5_range0.5.mat');

templatesL2_hist = T.templatesL2_hist;
%templatesL2_moms = T.templatesL2_moms;

%% Init data structure for responses

% set parameters for histogram computation at C1 & C2 layers
n_splits = T.n_splits; % the image is divided in a grid of n_splits x n_splits regions

%%%% NOTE: the # of bins should be made customizable between layers!
n_binsL1 = T.n_bins; % bars of the histograms at L1
n_binsL2 = T.n_bins; % bars of the histogram at L2
n_bins = T.n_bins; % bars of the histograms at L1
range = T.range; % range of the histogram

rangeL2 = [-0.001 0.001];

% init data structures which will contain the output signatures

S1uniform = cell(n_uniform,2);
C1uniform = cell(n_uniform,2);
S2uniform = cell(n_uniform,2);
C2uniform = cell(n_uniform,2);


%% S1 responses 

% output of dotproductL1_giulia is of format
% cell(n_scales, 1)
% responses{idx_scale} = zeros(n_templates, n_ori, sizeY-(taps(idx_scale)-1), sizeX-(taps(idx_scale)-1))

for idx_class = 1:2
for idx_uniform=1:n_uniform
    S1uniform{idx_uniform, idx_class} = dotproductL1_giulia(uniform_images{idx_uniform, idx_class}, gabors, n_splits);
end
end

%% C1 responses
 
% output of poolingL1_giulia is of format
% zeros(n_binsL1, n_reg, n_templatesL1, 2)

for idx_class = 1:2
for idx_uniform=1:n_uniform
   histograms = poolingL1_giulia(S1uniform{idx_uniform, idx_class}, n_splits, n_binsL1, range,  'histogram');
   signature = histograms(:, :, :, 1);
   signature = signature(:);
   C1uniform{idx_uniform, idx_class} = signature;
end
end

clear S1*

%% S2 responses 

% output of dotproductL2_giulia is of format
% zeros(signature_length, n_oriL2, n_scalesL2, n_templatesL2)
% signature_length = n_bins * n_reg * n_templatesL1

for idx_class = 1:2
for idx_uniform=1:n_uniform
   S2uniform{idx_uniform, idx_class} = dotproductL2_giulia(C1uniform{idx_uniform, idx_class}, templatesL2_hist);
end
end


%clear C1*
%% C2 responses

% output of poolingL2_giulia is of format
% zeros(n_bins, n_templatesL2, 2)

C2tot = cell(2, 1);

for idx_class = 1:2
for idx_uniform = 1:n_uniform
    histograms = poolingL2_giulia(S2uniform{idx_uniform, idx_class}, n_binsL2, rangeL2,  'histogram');
    signature = histograms(:, :, 1);
    signature = signature(:); 
    C2uniform{idx_uniform, idx_class} = signature';
    C2tot{idx_class} = [ C2tot{idx_class} ; C2uniform{idx_uniform, idx_class} ];
end
end

clear S2*

%% Response visualization and quantitative comaprison

for idx_class = 1:2
    % Create separate images for each class of objects
    figure( idx_class )
    m = mean(C2tot{idx_class});
    sd = std(C2tot{idx_class});
    f = [ m+3*sd , flipdim(m-3*sd,2)]; 
    fill([1:size(C2tot{idx_class},2) , size(C2tot{idx_class},2):-1:1] , f, [7 7 7]/8)
    hold on;
    plot(1:size(C2tot{idx_class},2) , m , 'b' , 'LineWidth',1);
    max(sd)
end

% Compute inter-class signature similarity

% Class 1
mSim1 = 0;
for i = 1:size(C2tot{1},1)
    for j = i+1:size(C2tot{1},1)

        mSim1 = mSim1 + similarity(C2tot{1}(i,:),C2tot{1}(j,:));
    end
end
mSim1 = mSim1*2/(size(C2tot{1},1)^2-size(C2tot{1},1))

% Class 2
mSim2 = 0;
for i = 1:size(C2tot{2},1)
    for j = i+1:size(C2tot{2},1)

        mSim2 = mSim2 + similarity(C2tot{2}(i,:),C2tot{2}(j,:));
    end
end
mSim2 = mSim2*2/(size(C2tot{2},1)^2-size(C2tot{2},1))

% Compute similarity between signatures of separate classes 

mSimCross = 0;
for i = 1:size(C2tot{1},1)
    for j = 1:size(C2tot{2},1)

        mSimCross = mSimCross + similarity(C2tot{1}(i,:),C2tot{2}(j,:));
    end
end
mSimCross = mSimCross/(size(C2tot{1,:},1)*size(C2tot{2,:},1))

%% Binary classifier

Y1 = ones(size(C2tot{1},1),1);
Y2 = -ones(size(C2tot{2},1),1);

% Randomly split the dataset between training and testing

empErrM = [];
empErrV = [];

for i = 2:2:50
    n_train = i;
    n_test = 70-i;

    % Run classifier some times
    numClassRuns = 100;
    missErr = [];
    for iClass = 1:numClassRuns

        [Xtr1, Ytr1, Xts1, Yts1] = randomSplitDataset(C2tot{1}, Y1, ceil(n_train/2), floor(n_test/2));
        [Xtr2, Ytr2, Xts2, Yts2] = randomSplitDataset(C2tot{2}, Y2, floor(n_train/2), floor(n_test/2));

        Xtr = [Xtr1 ; Xtr2];
        Ytr = [Ytr1 ; Ytr2];
        Xts = [Xts1 ; Xts2];
        Yts = [Yts1 ; Yts2];

        % Apply 1-NN classification
        k = 1;
        Ypred = kNNClassify(Xtr, Ytr, k, Xts);
        ind = find((sign(Ypred) ~= sign(Yts)));
        %[Yts  Ypred ]
        missErr = [ missErr (numel(ind)/n_test) ];
    end
    empErrM = [empErrM mean(missErr)]
    empErrV = [empErrV var(missErr,1)]
end

figure
plot(empErrM); hold on;
errorbar(empErrM , sqrt(empErrV),'rx');
hold off;