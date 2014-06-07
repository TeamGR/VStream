addpath((genpath('.')));

%% Simple layer

%% Import templates

%T = load('gabor_filters.mat');
T = load('pascal_filters.mat');

n_scales = length(T.templates);
n_templates = size(T.templates{1},1);
n_ori = size(T.templates{1},2);

taps = zeros(n_scales,1);
for idx_scale=1:n_scales
    taps(idx_scale) = size(T.templates{idx_scale},3);
end

templates = T.templates;

for idx_scale=1:n_scales
    
    min(min(min(min(templates{idx_scale}))))
    max(max(max(max(templates{idx_scale}))))
       
end

%% Import, normalize and zero-center the input image

image_path = fullfile('images', '000019.jpg');

inputImg = imread(image_path,'jpg');
if size(inputImg,3)==3
    inputImg = double(rgb2gray(inputImg));
end
inputImg  = inputImg - mean2(inputImg);
inputImg = inputImg / norm(inputImg);

%% Apply transformations to the input image

inputImg = imrotate(inputImg , 30 , 'bilinear' , 'crop');
[inSizeY, inSizeX] = size(inputImg);

%% Filter all the transformed images with all the transformed templates

S1responses = cell(n_scales, 1);

% 1) Loop on the scales
for idx_scale=1:n_scales
    
    S1responses{idx_scale,1} = zeros(n_templates, n_ori, inSizeY-(taps(idx_scale)-1), inSizeX-(taps(idx_scale)-1));

    for idx_template=1:n_templates
        for idx_ori=1:n_ori                     
            S1responses{idx_scale}(idx_template, idx_ori, :, :) = conv2(inputImg,squeeze(templates{idx_scale,1}(idx_template, idx_ori, :, :)),'valid');
        end
    end
    
    % should be already normalized if filter and image are: (v-mean(v))/norm(v-mean(v)) * (u-mean(u))/norm(u-mean(u))
    %min(min(min(min(S1responses{idx_scale}))))
    %max(max(max(max(S1responses{idx_scale}))))
    
end

%% Complex layer - Pooling

n_splits = 4;

n_bins = 10;
range = [-1 1];

C1responses = pooling_giulia(S1responses, n_splits, n_bins, range,  'histogram');

% 
% % Side sizes of the pooling regions
% regSideLenX = floor(inSizeX/poolingSplitNum);
% regSideLenY = floor(inSizeY/poolingSplitNum);
% regArea = regSideLenX * regSideLenY;
% 
% ma = max(max(cell2mat(filteredImg(1,1,1))));
% mi = min(min(cell2mat(filteredImg(1,1,1))));
% 
% step = (ma-mi)/numBars;
% x_hist = mi:step:ma;
% L1hist = cell(n_templates, poolRegsNum, 2);
% 
% % Templates loop
% for j = 1:n_templates
%     
%     % Reshape filtered image in column array form
%     reshFiltI = reshape(cell2mat(filteredImg(1,j,:)),numel(cell2mat(filteredImg(1,j,:))),1);
% 
%     % Pooling regions loop
%     for k = 1:poolRegsNum
%         [L1hist{j,k,1} L1hist{j,k,2}] = hist(reshFiltI((k-1)*regArea + 1 : k*regArea), x_hist);
%     end
% end
% Plot a histogram relative to a single template

bar(C1responses{8,1,2},C1responses{8,1,1});

% Generate pooling area 1 signature by concatenating histograms generated
% by the filtering with all templates
signature = horzcat(C1responses{:,1,1});
plot(signature);
