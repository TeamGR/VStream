%% Load images and infer parameters

load_prefixname = 'pascal';

input = load([load_prefixname '_images.mat'], 'images');
images = input.images;
input = load([load_prefixname '_translated_images.mat'], 'translated_images');
translated_images = input.images;
input = load([load_prefixname '_rotated_images.mat'], 'rotated_images');
rotated_images = input.images;
input = load([load_prefixname '_scaled_images.mat'], 'scaled_images');
scaled_images = input.images;

n_images = length(images);
n_xtranslations = size(translated_images, 2);
n_ytranslations = size(translated_images, 3);
n_rotations = size(rotated_images,2);
n_scales = size(scaled_images,2);

%% Load templates

T = load('gabor_filters.mat');
%T = load('pascal_filters.mat');

templates = T.templates;

%% Compute S1 responses

% set parameters for histogram computation at C1 layer
n_splits = 10; % the image is divided in a grid of n_splits x n_splits regions
n_bins = 10; % bars of the histogram
range = [-1 1]; % range of the histogram

% init data structures
S1 = cell{n_images, 1};
C1 = cell{n_images, 1};
S2 = cell{n_images, 1};
C2 = cell{n_images, 1};
S1transl = cell{n_images, n_xtranslations, n_ytranslations};
C1transl = cell{n_images, n_xtranslations, n_ytranslations};
S2transl = cell{n_images, n_xtranslations, n_ytranslations};
C2transl = cell{n_images, n_xtranslations, n_ytranslations};
S1rot = cell{n_images, n_rotations};
C1rot = cell{n_images, n_rotations};
S2rot = cell{n_images, n_rotations};
C2rot = cell{n_images, n_rotations};
S1scale = cell{n_images, n_scales};
C1scale = cell{n_images, n_scales};
S2scale = cell{n_images, n_scales};
C2scale = cell{n_images, n_scales};

for idx_image=1:n_images
    S1{idx_image} = convolve_giulia(images{idx_image}, templates);
end

for ix_transl=1:n_xtranslations
    for iy_transl=1:n_ytranslations
        S1transl{idx_image, ix_transl, iy_transl} = convolve_giulia(translated_images{idx_image, ix_transl, iy_transl}, templates);
    end
end

for idx_rot=1:n_rotations
    S1rot{idx_image, idx_rot} = convolve_giulia(rotated_images{idx_image, idx_rot}, templates);
end

for idx_scale=1:n_scales
    S1scale{idx_image, idx_scale} = convolve_giulia(scaled_images{idx_image, idx_scale}, templates);
end

%% Compute C1 responses
 
for idx_image=1:n_images
   histograms = pooling_giulia(S1{idx_image}, n_splits, n_bins, range,  'histogram');
   signature = histograms(:, :, :, 1);
   C1{idx_image} = signature(:);
end

for ix_transl=1:n_xtranslations
    for iy_transl=1:n_ytranslations
        histograms = pooling_giulia(S1transl{idx_image, ix_transl, iy_transl}, n_splits, n_bins, range,  'histogram');
        signature = histograms(:, :, :, 1);
        C1transl{idx_image, ix_transl, iy_transl} = signature(:);
    end
end

for idx_rot=1:n_rotations
    histograms = pooling_giulia(S1rot{idx_image, idx_rot}, n_splits, n_bins, range,  'histogram');
    signature = histograms(:, :, :, 1);
    C1rot{idx_image, idx_rot} =  signature(:);
end

for idx_scale=1:n_scales
    histograms = pooling_giulia(S1scale{idx_image, idx_scale}, n_splits, n_bins, range,  'histogram');
    signature = histograms(:, :, :, 1);
    C1scale{idx_image, idx_scale} = signature(:);
end



C1responses_hist = pooling_giulia(S1responses, n_splits, n_bins, range,  'histogram'); % out_hist = zeros(n_bins, n_reg, n_templates, 2)
% C1responses_moms = pooling_giulia(S1responses, n_splits, n_bins, range,  'moments'); % out_moms = zeros(2, n_reg, n_templates)

% Concatenate histograms to obtain a 1-D array of size n_bins*n_reg*n_templates
C1responses_hist_signature = C1responses_hist(:, :, :, 1);
C1responses_hist_signature = C1responses_hist_signature(:);
% C1responses_hist_xdomain = C1responses_hist(:, :, :, 2);
% C1responses_hist_xdomain = C1responses_hist_xdomain(:);

% Concatenate moments to obtain a 1-D array of size 3*n_reg*n_templates
% C1responses_moms_signature = C1responses_moms(:);

% Plot the histogram signature
figure
plot(C1responses_hist_signature, '.');
title('C1 signature histogram')

% Plot the moments signature
figure
plot(C1responses_moms_signature, '.');
title('C1 signature 3 moments')
