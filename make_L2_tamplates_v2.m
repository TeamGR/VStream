addpath((genpath('.')));

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

templates = T.templates;


% Load raw L2 templates

T2 = load('pascal_filters.mat');
%templatesL2_raw = T2.templates;

%----------- DEBUG
templatesL2_raw = cell(3,1);
for i=1:3
    templatesL2_raw{i} = ones(1,1,100,100);
end
%-----------------

% Format of templatesL2_raw:
% n_scales_raw X 1 cell array
%   Each cell contains:
%   n_templates X n_ori 2D raw templates matrices

n_scales_raw = length(templatesL2_raw);
n_templates_raw = size(templatesL2_raw{1},1);
n_ori_raw = size(templatesL2_raw{1},2);

%% Filter the raw L2 templates with the L1 templates, then pool and concatenate to obtain L2 templates

templatesL2_hist = cell(n_templates_raw , n_ori_raw , n_scales_raw);
templatesL2_moms = cell(n_templates_raw , n_ori_raw , n_scales_raw);

% Outer loop on the raw templates

for idx_templates_raw = 1:n_templates_raw
for idx_ori_raw = 1:n_ori_raw
for idx_scales_raw = 1:n_scales_raw

S1responses = cell(n_scales, 1);
[inSizeY , inSizeX] = size(squeeze(templatesL2_raw{idx_scales_raw,1}(1,1,:,:)));

% Inner loops
for idx_scale=1:n_scales            % Scales
S1responses{idx_scale,1} = zeros(n_templates, n_ori, inSizeY - (taps(idx_scale)-1) , inSizeX - (taps(idx_scale)-1));

    for idx_template=1:n_templates      % Templates
        for idx_ori=1:n_ori                 % Orientations                     
            S1responses{idx_scale,1}(idx_template, idx_ori, :, :) = conv2(squeeze(templatesL2_raw{idx_scales_raw , 1}( idx_templates_raw , idx_ori_raw , : , : )) , squeeze(templates{idx_scale,1}(idx_template, idx_ori, :, :)),'valid');
        end
    end
end

% Complex layer - Pooling

n_splits = 10;

n_bins = 10;
range = [-1 1];

C1responses_hist = pooling_giulia(S1responses, n_splits, n_bins, range,  'histogram'); % out_hist = zeros(n_bins, n_reg, n_templates, 2)
C1responses_moms = pooling_giulia(S1responses, n_splits, n_bins, range,  'moments'); % out_moms = zeros(2, n_reg, n_templates)
    
% Concatenate histograms to obtain a 1-D array of size n_bins*n_reg*n_templates
templatesL2_hist{idx_templates_raw , idx_ori_raw , idx_scales_raw} = C1responses_hist(:, :, :, 1);
templatesL2_hist{idx_templates_raw , idx_ori_raw , idx_scales_raw} = templatesL2_hist{idx_templates_raw , idx_ori_raw , idx_scales_raw}(:);
% C1responses_hist_xdomain = C1responses_hist(:, :, :, 2);
% C1responses_hist_xdomain = C1responses_hist_xdomain(:);

% Concatenate moments to obtain a 1-D array of size 3*n_reg*n_templates
templatesL2_moms{idx_templates_raw , idx_ori_raw , idx_scales_raw} = C1responses_moms(:);
    
    
end
end
end

