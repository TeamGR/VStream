%% Experiment description
%
% Number of layers of the net: 2
% Employed dataset: SUFR
% Transformations:
%   Independent scaling, translation and rotation of faces on a uniform
%   background
%
%

close all;
clear all;

%% Load images and infer parameters

% Rotated images
n_rotations = 7;
Seq = cell(n_rotations);
imgPath = 'SUFRData\image_files\uniform_bg\roll';
imgType = '*.jpg'; % change based on image type
images  = dir([imgPath '\' imgType]);
for idx = 1:n_rotations
    Seq{idx} = double(rgb2gray(imread([imgPath '\' images(idx).name])));
end
rotated_images = Seq;

% Translated Images
n_translations = 5;
Seq = cell(n_translations);
imgPath = 'SUFRData\image_files\uniform_bg\translation';
imgType = '*.jpg'; % change based on image type
images  = dir([imgPath '\' imgType]);
for idx = 1:n_translations
    Seq{idx} = double(rgb2gray(imread([imgPath '\' images(idx).name])));
end
translated_images = Seq;

% Scaled images
n_scales = 5;
Seq = cell(n_scales);
imgPath = 'SUFRData\image_files\uniform_bg\scaling';
imgType = '*.jpg'; % change based on image type
images  = dir([imgPath '\' imgType]);
for idx = 1:n_scales
    Seq{idx} = double(rgb2gray(imread([imgPath '\' images(idx).name])));
end
scaled_images = Seq;

%% Load templates of layer 1

T = load('gabor_filters.mat');
%T = load('pascal_filters.mat');

gabors = T.templates;

%% Load templates of layer 2

T = load('templatesL2.mat');

templatesL2_hist = T.templatesL2_hist;
templatesL2_moms = T.templatesL2_moms;


%% Init data structure for responses

% set parameters for histogram computation at C1 layer
n_splits = T.n_splits; % the image is divided in a grid of n_splits x n_splits regions
n_bins = T.n_bins; % bars of the histogram
range = [-1 1]; % range of the histogram

% init data structures which will contain the output signatures
S1transl = cell(n_translations);
C1transl = cell(n_translations);
S2transl = cell(n_translations);
C2transl = cell(n_translations);
S1rot = cell(n_rotations);
C1rot = cell(n_rotations);
S2rot = cell(n_rotations);
C2rot = cell(n_rotations);
S1scale = cell(n_scales);
C1scale = cell(n_scales);
S2scale = cell(n_scales);
C2scale = cell(n_scales);

%% S1 responses 

for idx_rot=1:n_rotations
    S1rot{idx_rot} = convolve_giulia(rotated_images{idx_rot}, gabors);
end

for idx_transl=1:n_translations
    S1transl{idx_transl} = convolve_giulia(translated_images{idx_transl}, gabors);
end

for idx_scale=1:n_scales
    S1scale{idx_scale} = convolve_giulia(scaled_images{idx_scale}, gabors);
end

%% C1 responses
 
for idx_rot=1:n_rotations
    histograms = poolingL1_giulia(S1rot{idx_rot}, n_splits, n_bins, range,  'histogram');
    signature = histograms(:, :, :, 1);
    normalized_signature = signature - mean(signature);
    normalized_signature = normalized_signature / norm(normalized_signature);
    C1rot{idx_rot} = normalized_signature;
end

for ix_transl=1:n_translations
    histograms = poolingL1_giulia(S1transl{idx_transl}, n_splits, n_bins, range,  'histogram');
    signature = histograms(:, :, :, 1);
    signature = signature(:);
    normalized_signature = signature - mean(signature);
    normalized_signature = normalized_signature / norm(normalized_signature);
    C1transl{idx_transl} = normalized_signature;
end

for idx_scale=1:n_scales
    histograms = poolingL1_giulia(S1scale{idx_scale}, n_splits, n_bins, range,  'histogram');
    signature = histograms(:, :, :, 1);
    normalized_signature = signature - mean(signature);
    normalized_signature = normalized_signature / norm(normalized_signature);
    C1scale{idx_scale} = normalized_signature;
end

%% S2 responses 

% output of dot_product_giulia is of format
% zeros(signature_length, n_ori, n_scales, n_templates);

for idx_rot=1:n_rotations
    S2rot{idx_rot} = dotproduct_giulia(C1rot{idx_rot}, templatesL2_hist);
end

for ix_transl=1:n_translations
    S2transl{idx_transl} = dotproduct_giulia(C1transl{idx_transl}, templatesL2_hist);
end

for idx_scale=1:n_scales
    S2scale{idx_scale} = dotproduct_giulia(C1scale{idx_scale}, templatesL2_hist);
end

%% C2 responses

for ix_transl=1:n_translations
        C2transl{idx_transl} = poolingL2_giulia(S2transl{idx_transl}, n_splits, n_bins, range,  'histogram');
end

for idx_rot=1:n_rotations
    C2rot{idx_rot} = poolingL2_giulia(S2rot{idx_rot}, n_splits, n_bins, range,  'histogram');
end

for idx_scale=1:n_scales
    C2scale{idx_scale} = poolingL2_giulia(S2scale{idx_scale}, n_splits, n_bins, range,  'histogram');
end