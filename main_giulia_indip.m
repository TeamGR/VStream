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

%% S1 responses 

for idx_image=1:n_images
    S1{idx_image} = convolve_giulia(images{idx_image}, gabors);
end

for ix_transl=1:n_xtranslations
    for iy_transl=1:n_ytranslations
        S1transl{idx_image, ix_transl, iy_transl} = convolve_giulia(translated_images{idx_image, ix_transl, iy_transl}, gabors);
    end
end

for idx_rot=1:n_rotations
    S1rot{idx_image, idx_rot} = convolve_giulia(rotated_images{idx_image, idx_rot}, gabors);
end

for idx_scale=1:n_scales
    S1scale{idx_image, idx_scale} = convolve_giulia(scaled_images{idx_image, idx_scale}, gabors);
end

%% C1 responses
 
for idx_image=1:n_images
   histograms = poolingL1_giulia(S1{idx_image}, n_splits, n_bins, range,  'histogram');
   signature = histograms(:, :, :, 1);
   signature = signature(:);
   normalized_signature = signature - mean(signature);
   normalized_signature = normalized_signature / norm(normalized_signature);
   C1{idx_image} = normalized_signature;
end

for ix_transl=1:n_xtranslations
    for iy_transl=1:n_ytranslations
        histograms = poolingL1_giulia(S1transl{idx_image, ix_transl, iy_transl}, n_splits, n_bins, range,  'histogram');
        signature = histograms(:, :, :, 1);
        signature = signature(:);
        normalized_signature = signature - mean(signature);
        normalized_signature = normalized_signature / norm(normalized_signature);
        C1transl{idx_image, ix_transl, iy_transl} = normalized_signature;
    end
end

for idx_rot=1:n_rotations
    histograms = poolingL1_giulia(S1rot{idx_image, idx_rot}, n_splits, n_bins, range,  'histogram');
    signature = histograms(:, :, :, 1);
    normalized_signature = signature - mean(signature);
    normalized_signature = normalized_signature / norm(normalized_signature);
    C1rot{idx_image, idx_rot} = normalized_signature;
end

for idx_scale=1:n_scales
    histograms = poolingL1_giulia(S1scale{idx_image, idx_scale}, n_splits, n_bins, range,  'histogram');
    signature = histograms(:, :, :, 1);
    normalized_signature = signature - mean(signature);
    normalized_signature = normalized_signature / norm(normalized_signature);
    C1scale{idx_image, idx_scale} = normalized_signature;
end

%% S2 responses 

% output of dot_product_giulia is of format
% zeros(signature_length, n_ori, n_scales, n_templates);

for idx_image=1:n_images
   S2{idx_image} = dotproduct_giulia(C1{idx_image}, templatesL2_hist);
end

for ix_transl=1:n_xtranslations
    for iy_transl=1:n_ytranslations
        S2transl{idx_image, ix_transl, iy_transl} = dotproduct_giulia(C1transl{idx_image, ix_transl, iy_transl}, templatesL2_hist);
    end
end

for idx_rot=1:n_rotations
    S2rot{idx_image, idx_rot} = dotproduct_giulia(C1rot{idx_image, idx_rot}, templatesL2_hist);
end

for idx_scale=1:n_scales
    S2scale{idx_image, idx_scale} = dotproduct_giulia(C1scale{idx_image, idx_scale}, templatesL2_hist);
end

%% C2 responses

for idx_image=1:n_images
   C2{idx_image} = poolingL2_giulia(S2{idx_image}, n_splits, n_bins, range,  'histogram');
end

for ix_transl=1:n_xtranslations
    for iy_transl=1:n_ytranslations
        C2transl{idx_image, ix_transl, iy_transl} = poolingL2_giulia(S2transl{idx_image, ix_transl, iy_transl}, n_splits, n_bins, range,  'histogram');
    end
end

for idx_rot=1:n_rotations
    C2rot{idx_image, idx_rot} = poolingL2_giulia(S2rot{idx_image, idx_rot}, n_splits, n_bins, range,  'histogram');
end

for idx_scale=1:n_scales
    C2scale{idx_image, idx_scale} = poolingL2_giulia(S2scale{idx_image, idx_scale}, n_splits, n_bins, range,  'histogram');
end