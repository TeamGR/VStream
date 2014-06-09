addpath((genpath('.')));

%% Set transformations for input images

xtranslations = [-10 10];
ytranslations = [-10 10];

rotations = pi/3:pi/3:2*pi;

scales = [2 4];

n_xtranslations = length(xtranslations);
n_ytranslations = length(ytranslations);
n_translations = n_xtranslations*n_ytranslations;
n_rotations = length(rotations);
n_scales = length(scales);

%% Set registry file with image paths

% set registry path
registry_path = fullfile(pwd, 'pascal_registry_train.txt');
registry_fid = fopen(registry_path,'r');
if (registry_fid==-1)
    error('Error! Please provide a valid path for registry file.');
end

% set extension, images folder and number of different images
ext = 'jpg';
images_path = 'D:\IIT\CODICI\pipeline_overfeat\PASCAL2007\VOCdevkit\VOC2007\JPEGImages';
n_images = 5;

% set prefix name for saving transformed images
save_prefixname = 'pascal';

% count lines in registry file and check if there are enough image paths
if isunix
    [~,n] = system(['wc -l < ' registry_path]);
    n = str2num(n);
elseif ispc
    [~,n] = system(['find /v /c "&*fake&*" "' registry_path '"']);
    last_space = strfind(n,' ');
    n = str2num(n((last_space(end)+1):(end-1)));
end
if n<n_images
    error('Error! Not enough image paths in registry file.');
end

%% Read input images and apply transformations

% init data structures for images and responses
images = cell(n_images, 1);
translated_images = cell(n_images, n_xtranslations, n_ytranslations);
rotated_images = cell(n_images, n_rotations);
scaled_images = cell(n_images, n_scales);

for idx_image=1:n_images
    
    % read image
    image_path = fullfile(images_path, [fgetl(registry_fid) '.' ext]);
    I = imread(image_path, ext);
    % convert to grayscale if necessary
    if (size(I,3) == 3)
        I = rgb2gray(I);
    end
    % convert to double and assign
    images{idx_image} = double(I);
  
    % apply translations
    for ix_transl = 1:n_xtranslations
        for iy_transl = 1:n_ytranslations
            translated_images{idx_image, ix_transl, iy_transl} = imtranslate(images{idx_image},[xtranslations(ix_transl) ytranslations(iy_transl)], 'OutputView', 'same', 'FillValues', mean2(images{idx_image}));
        end
    end
    
    % apply rotations
    for idx_rot=1:n_rotations
        rotated_images{idx_image, idx_rot} = imrotate(images{idx_image}, rotations(idx_rot)*180/pi, 'bilinear' , 'crop');
    end
    
    % apply scales
    for idx_scale=1:n_scales
        scaled_images{idx_image, idx_scale} = imresize(images{idx_image}, size(images{idx_image})*scales(idx_scale));
    end
   
end

save([save_prefixname '_images.mat'], 'images');
save([save_prefixname '_translated_images.mat'], 'translated_images');
save([save_prefixname '_rotated_images.mat'], 'rotated_images');
save([save_prefixname '_scaled_images.mat'], 'scaled_images');