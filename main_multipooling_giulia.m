addpath((genpath('.')));

%% Simple layer

%% Import templates

T = load('gabor_filters.mat');
%T = load('pascal_filters.mat');

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
inputImg = imresize(inputImg, [100 100]);
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

n_splits = 10;

n_bins = 10;
range = [-1 1];

C1responses_hist = pooling_giulia(S1responses, n_splits, n_bins, range,  'histogram'); % out_hist = zeros(n_bins, n_reg, n_templates, 2)
C1responses_moms = pooling_giulia(S1responses, n_splits, n_bins, range,  'moments'); % out_moms = zeros(2, n_reg, n_templates)

% Concatenate histograms to obtain a 1-D array of size n_bins*n_reg*n_templates
C1responses_hist_signature = C1responses_hist(:, :, :, 1);
C1responses_hist_signature = C1responses_hist_signature(:);
C1responses_hist_xdomain = C1responses_hist(:, :, :, 2);
C1responses_hist_xdomain = C1responses_hist_xdomain(:);

% Concatenate moments to obtain a 1-D array of size 3*n_reg*n_templates
C1responses_moms_signature = C1responses_moms(:);

% Plot the histogram signature
figure
plot(C1responses_hist_signature);
title('C1 signature histogram')

% Plot the moments signature
figure
plot(C1responses_moms_signature);
title('C1 signature 3 moments')
