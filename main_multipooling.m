clear all;
close all;




%% Random patches

%% Parameters

n_ori = 8;
res_ori = 2*pi / n_ori;
ori = (0:(n_ori-1))*res_ori;

n_scales = 1;
scale_factor = 2;
sx_scale0 = 22;
sy_scale0 = 22;
sx_scales = sx_scale0 * ones(1, n_scales);
sy_scales = sy_scale0 * ones(1, n_scales);

n_transformations = n_ori * n_scales;

n_templates = 10;

n_filters = n_templates * n_transformations;

%% Source images

registry_path = fullfile(pwd, 'pascal_registry.txt');
ext = '.jpg';
images_path = '\..\PASCAL2007\VOCdevkit\VOC2007\JPEGImages';

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

%% Patch extraction and transformations

templates = cell(n_scales,1);

for idx_scale=1:n_scales
    
    templates{idx_scale,1} = zeros(n_templates, n_ori, sx_scales(idx_scale), sy_scales(idx_scale));
    
    for idx_template=1:n_templates
        
        % generation of the base template
        
        % 1) import image
        template_path = fgetl(registry_fid);
        template_path = fullfile(images_path, [template_path ext]);
        I = imread(template_path);
        if (size(I,3) == 3)
           I = rgb2gray(I);
        end
        sx_I = size(I,2); % x-dim along horizontal axis
        sy_I = size(I,1); % y-dim along vertical axis
        
        % 2) extract random x and y coordinates (patch's left-upper corner)
        x_crop = unidrnd(sx_I-sx_scales(idx_scale));
        y_crop = unidrnd(sy_I-sy_scales(idx_scale));
        
        % 3) crop the image at extracted x and y
        basic_template = double(I(y_crop:(y_crop+sy_scales(idx_scale)-1), x_crop:(x_crop+sx_scales(idx_scale)-1)));
        
        % Normalize and zero-mean the template
        basic_template  = basic_template - mean(mean(basic_template));
        basic_template = basic_template ./ norm(basic_template, 1);

        %         figure(1)
%         subplot(10, 10, idx_template);
%         imshow(basic_template, [0 255]);
        
        % generations of its transformations
        
        % 1) loop on the orientations
        for idx_dir=1:n_ori
            
            templates{idx_scale,1}(idx_template, idx_dir, :, :) = imrotate(basic_template, ori(idx_dir)*180/pi, 'nearest', 'crop');
        
%             figure(idx_template+1)
%             subplot(2, 4, idx_dir);
%             imshow(squeeze(templates{idx_scale,1}(idx_template, idx_dir, :, :)), []);
        end
   
    end
end

fclose(registry_fid);

%% Gabors 



%% Parameters


%% Import, normalize and zero-center the input image
inputImg = double(rgb2gray(imread('lena.jpg','jpg')));
inputImg  = inputImg - mean(mean(inputImg));
inputImg = inputImg ./ norm(inputImg, 1);
inputImg = imrotate(inputImg , 30 , 'bilinear' , 'crop');
[inSizeX inSizeY] = size(inputImg);

%% Simple layer

% Initialize filtered images array
filteredImg = cell(n_scales, n_templates, n_ori);

for idx_scale=1:n_scales

    % 1) Loop on the templates
    for idx_template=1:n_templates
        
        % 2) Loop on the orientations
        for idx_dir=1:n_ori
                            
            filteredImg{idx_scale, idx_template, idx_dir} = conv2(inputImg,squeeze(templates{idx_scale,1}(idx_template, idx_dir, :, :)),'valid');
        end
    end
end

figure;
imshow(cell2mat(filteredImg(1,1,1)) , []);

%% Complex layer

% Image pooling
poolingSplitNum = 1;    % Number of splits along axes for determining pooling regions
poolRegsNum = poolingSplitNum^2; 
numBars = 20;   % Number of bars of the histogram

[ L1hist ] = pooling( filteredImg , n_templates , poolingSplitNum , numBars , 'histogram');

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

bar(L1hist{8,1,2},L1hist{8,1,1});

% Generate pooling area 1 signature by concatenating histograms generated
% by the filtering with all templates
signature = horzcat(L1hist{:,1,1});
plot(signature);
    

%% Response visualization and quantitative comaprison