addpath(genpath('.'));

%% Random patches

% parameters

n_ori = 8;
res_ori = 2*pi / n_ori;
ori = (0:(n_ori-1))*res_ori;

n_scales = 3;
taps = [11 23 47];

n_transformations = n_ori * n_scales;

n_templates = 10;

n_filters = n_templates * n_transformations;

% source images

registry_path = fullfile(pwd, 'pascal_registry.txt');
ext = '.jpg';
images_path = 'D:\IIT\CODICI\pipeline_overfeat\PASCAL2007\VOCdevkit\VOC2007\JPEGImages';

registry_fid = fopen(registry_path,'r');
if (registry_fid==-1)
    error('Error! Please provide a valid path for registry file.');
end

% count lines in registry file and check if there are enough images
if isunix
    [~,n_images] = system(['wc -l < ' registry_path]);
    n_images = str2num(n_images);
elseif ispc
    [~,n_images] = system(['find /v /c "&*fake&*" "' registry_path '"']);
    last_space = strfind(n_images,' ');
    n_images = str2num(n_images((last_space(end)+1):(end-1)));
end
if n_images<n_templates
    error('Error! Not enough image paths in registry file.');
end

% patch extraction and transformations

templates = cell(n_scales,1);

for idx_scale=1:n_scales
    templates{idx_scale,1} = zeros(n_templates, n_ori, taps(idx_scale), taps(idx_scale));
end

for idx_template=1:n_templates
    
    % import image
    template_path = fgetl(registry_fid);
    template_path = fullfile(images_path, [template_path ext]);
    I = imread(template_path);
    if (size(I,3) == 3)
        I = rgb2gray(I);
    end
    sx_I = size(I,2); % x-dim along horizontal axis
    sy_I = size(I,1); % y-dim along vertical axis
    
    % extract random x and y coordinates (patch's center)
    x_crop = unidrnd(sx_I-taps(idx_scale));
    y_crop = unidrnd(sy_I-taps(idx_scale));
    
    % compute x and y ranges for the crop
    x_min = x_crop - (taps(1)-1)/2;
    x_max = x_crop + (taps(1)-1)/2;
    y_min = y_crop - (taps(1)-1)/2;
    y_max = y_crop + (taps(1)-1)/2;
    
    % rotate the image around the patch's center and pick the crop 
    figure(idx_template)
    for idx_ori=1:n_ori 
        
        rotated_image = rotateAround(I, y_crop, x_crop, ori(idx_ori)*180/pi, 'bilinear');

        rotated_template = rotated_image(y_min:y_max,x_min:x_max);

        for idx_scale=1:n_scales
            
            rotated_scaled_template = imresize(squeeze(templates{1,1}(idx_template, idx_ori, :, :)),[taps(idx_scale) taps(idx_scale)]);
            templates{idx_scale,1}(idx_template, idx_ori, :, :) = (rotated_scaled_template-mean2(rotated_scaled_template))/norm(rotated_scaled_template-mean2(rotated_scaled_template));  
            subplot(n_scales, n_ori, (idx_scale-1)*n_ori+idx_ori);
            imagesc(squeeze(templates{idx_scale,1}(idx_template, idx_ori, :, :)));
            colormap(gray);
        end 
   
    end

end

fclose(registry_fid);

save('pascal_filters.mat', 'templates');

%% Gabors 

% parameters

n_ori = 8;
res_ori = pi / n_ori;
ori = (0:(n_ori-1))*res_ori;

n_scales = 3;
taps = [11 23 47];

n_transformations = n_ori * n_scales;

n_templates = 2; % even and odd parts

n_filters = n_templates * n_transformations;

% filters computation

templates = cell(n_scales,1);

f0 = 1/4;
B = f0/3;
for idx_scale=1:n_scales
    
    templates{idx_scale,1} = zeros(n_templates, n_ori, taps(idx_scale), taps(idx_scale));
    
    F = design_gabor_filt(f0,taps(idx_scale),B,0);
    [E,O] = compose_gabor_filters(F);
   
    figure
    for idx_ori=1:n_ori
            
        templates{idx_scale,1}(1, idx_ori, :, :) = (E(:,:,idx_ori)-mean2(E(:,:,idx_ori)))/norm(O(:,:,idx_ori)-mean2(O(:,:,idx_ori)));
        templates{idx_scale,1}(2, idx_ori, :, :) = (O(:,:,idx_ori)-mean2(O(:,:,idx_ori)))/norm(O(:,:,idx_ori)-mean2(O(:,:,idx_ori)));
        
        subplot(2, n_ori, idx_ori);
        imagesc(squeeze(templates{idx_scale,1}(1, idx_ori, :, :)));
        subplot(2, n_ori, idx_ori+8);
        imagesc(squeeze(templates{idx_scale,1}(2, idx_ori, :, :)));
    end
    
    f0 = f0/2;
    B = f0/3;
end

save('gabor_filters.mat', 'templates');

%% Filters loading

templates_gabor = load('gabor_filters.mat');
templates_pascal = load('pascal_filters.mat');