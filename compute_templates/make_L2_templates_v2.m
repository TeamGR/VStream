addpath((genpath('.')));
clear all;
close all;

%% Simple layer

% Import templates

% Load cooked L1 templates

T = load('gabor_filters.mat');

n_scales = length(T.templates);
n_templates = size(T.templates{1},1);
n_ori = size(T.templates{1},2);

taps = zeros(n_scales,1);
for idx_scale=1:n_scales
    taps(idx_scale) = size(T.templates{idx_scale},3);
end

gabors = T.templates;

% Load raw L2 templates

load('pascal_templates_L2_raw.mat');

% Format of templatesL2_raw:
% n_scales_raw X 1 cell array
%   Each cell contains:
%   n_templates X n_ori 2D raw templates matrices

n_scales_raw = length(templatesL2_raw);
n_templates_raw = size(templatesL2_raw{1},1);
n_ori_raw = size(templatesL2_raw{1},2);

% Pooling parameters
n_splits = 2;
n_bins = 20;
range = [-1 1];

%% Filter the raw L2 templates with the L1 templates, then pool and concatenate to obtain L2 templates

templatesL2_hist = cell(n_templates_raw , n_ori_raw , n_scales_raw);
templatesL2_moms = cell(n_templates_raw , n_ori_raw , n_scales_raw);

% Outer loop on the raw templates



% Inner loops
S1responses = cell(n_templates, n_ori, n_scales);
for idx_scales_raw=1:n_scales_raw            % Scales
%[inSizeY , inSizeX] = size(squeeze(templatesL2_raw{idx_scales_raw,1}(1,1,:,:)));
    for idx_templates_raw=1:n_templates_raw      % Templates
        for idx_ori_raw=1:n_ori_raw                 % Orientations           
            S1responses{idx_templates_raw , idx_ori_raw , idx_scales_raw} =...
            dotproductL1_giulia(squeeze(templatesL2_raw{idx_scales_raw , 1}...
            ( idx_templates_raw , idx_ori_raw , : , : )) , gabors , n_splits);
        end
    end

    idx_templates_raw
    idx_ori_raw
    idx_scales_raw
end

% Complex layer - Pooling

C1responses_hist = cell(n_templates, n_ori, n_scales);
for idx_scales_raw=1:n_scales_raw            % Scales
%S1responses{idx_scale,1} = zeros(n_templates, n_ori, inSizeY - (taps(idx_scale)-1) , inSizeX - (taps(idx_scale)-1));
    for idx_templates_raw=1:n_templates_raw      % Templates
        for idx_ori_raw=1:n_ori_raw                 % Orientations           
 
            histograms = poolingL1_giulia(S1responses{idx_templates_raw , idx_ori_raw , idx_scales_raw}, n_splits, n_bins, range,  'histogram');

            % Concatenate histograms to obtain a 1-D cell array of size n_bins*n_reg*n_templates
            signature = histograms(:, :, :, 1);
            signature = signature(:); 
            C1responses_hist{idx_templates_raw, idx_ori_raw, idx_scales_raw} = signature;
        end
    end

    idx_templates_raw
    idx_ori_raw
    idx_scales_raw
end


save ('templatesL2.mat', 'templatesL2_hist' , 'n_splits' , 'n_bins' );