%% HW module

n_dirs = 8;
res_dir = 2*pi / n_dirs;
dirs = (0:(n_dirs-1))*res_dir;

n_transl_per_dir = 3; % symmetric
half_transl_per_dir = (n_transl_per_dir-1)/2 ;
res_transl_per_dir = 3; % [pixels]
transls_per_dir = [(-half_transl_per_dir:half_transl_per_dir)*res_transl_per_dir; zeros(n_transl_per_dir,1)];

n_transls = n_dirs * n_transl_per_dir;
transls = zeros(2, n_transls);
for idx_dir=1:n_dirs
    
    theta = dirs(idx_dir);
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    transls(:, ((idx_dir-1)*n_transls_per_dir+1):(idx_dir*n_transls_per_dir)) = R*transls_per_dir;

end

n_scales = 1;
scale_factor = 2;
sx_scale0 = 100;
sy_scale0 = 100;
sx_scales = sx_scale0 * ones(1, n_scales);
sy_scales = sy_scale0 * ones(1, n_scales);

n_transformations = n_transls * n_dirs * n_scales;

n_templates = 100;
registry_path = fullfile(pwd, 'pascal_registry.txt');
ext = '.jpg';
images_path = 'D:\IIT\CODICI\pipeline_overfeat\PASCAL2007\VOCdevkit\VOC2007\JPEGImages';

n_RFs = n_templates * n_transformations;

%% Templates construction

% Open registry file
registry_fid = fopen(registry_path,'r');
if (registry_fid==-1)
    error('Error! Please provide a valid path for registry file.');
end

% Count lines in registry file and check if there are enough images
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

templates = cell(n_scales,1);

for idx_scale=1:n_scales
    
    templates{idx_scale,1} = zeros(n_templates, n_dirs, n_transls, sx_scales(idx_scale), sy_scales(idx_scale));
    
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
        basic_template = I(y_crop:(y_crop+sy_scales(idx_scale)-1), x_crop:(x_crop+sx_scales(idx_scale)-1));
        
        figure(1)
        subplot(10, 10, idx_template);
        imshow(basic_template, [0 255]);
        
        % generations of its transformations
        
        % 1) loop on the orientations
        for idx_dir=1:n_dirs
            
            rotated_template = imrotate(basic_template, rotations(idx_dir)*180/pi, 'nearest', 'crop');
            
            % 2) loop on the translations
            for idx_translation=1:n_transls
                
                templates{idx_scale,1}(idx_template, idx_dir, idx_translation, :, :) = rotated_template;
            end
            
        end
   
    end
end

fclose(registry_fid);

%% Simple layer

%% Complex layer

%% 2Response visualization and quantitative comaprison