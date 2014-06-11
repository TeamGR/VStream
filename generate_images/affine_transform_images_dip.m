%addpath((genpath('.')));

%% Set transformations for input images

% xtranslations = [-10 10];
% ytranslations = [-10 10];

rotations = pi/3:pi/3:2*pi;

scales = [1 1.5 2];

% n_xtranslations = length(xtranslations);
% n_ytranslations = length(ytranslations);
% n_translations = n_xtranslations*n_ytranslations;

n_rotations = length(rotations);
n_scales = length(scales);

%% Set registry file with image paths

% set registry path
registry_path = fullfile(pwd, 'pascal_registry_L2.txt');
registry_fid = fopen(registry_path,'r');
if (registry_fid==-1)
    error('Error! Please provide a valid path for registry file.');
end

% set extension, images folder and number of different images
ext = 'jpg';
images_path = '..\PASCAL2007\VOCdevkit\VOC2007\JPEGImages';
n_images = 20;

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

% init data structure for input and output images
images = cell(n_images, 1);

% Format of templatesL2_raw:
% n_scales_raw X 1 cell array
%   Each cell contains:
%   n_templates X n_ori 2D raw templates matrices

templatesL2_raw = cell(n_scales,1);
% for i=1:3
%     templatesL2{i} = ones(1,1,100,100);
% end

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
  
     % Crop image
     images{idx_image} = imcrop(images{idx_image}, [ 0 0 ceil(256*sqrt(2)) ceil(256*sqrt(2)) ]);
    
    
    % apply rotations
    for idx_rot=1:n_rotations

        tmp = imrotate(images{idx_image}, rotations(idx_rot)*180/pi, 'bilinear' , 'crop');

        % apply scaling    
        for idx_scale=1:n_scales

            tmp = imresize( tmp, size(images{idx_image})*scales(idx_scale));

            % Crop image
            templatesL2_raw{idx_scale}(idx_image, idx_rot , : , : ) = imcrop(tmp, [ floor((size(tmp,1) - 256 ) /2) floor((size(tmp,2)-256)/2) 256 256 ]);

            % Centering and normalization
            templatesL2_raw{idx_scale}(idx_image, idx_rot) = templatesL2_raw{idx_scale}(idx_image, idx_rot) - mean2(templatesL2_raw{idx_scale}(idx_image, idx_rot));
            templatesL2_raw{idx_scale}(idx_image, idx_rot) = templatesL2_raw{idx_scale}(idx_image, idx_rot) / norm(templatesL2_raw{idx_scale}(idx_image, idx_rot));
        end
    end

    % AT THE END normalize original image
%     images{idx_image} = images{idx_image} - mean2(images{idx_image});
%     images{idx_image} = images{idx_image} / norm(images{idx_image});
end

save('pascal_templates_L2.mat', 'templatesL2_raw');