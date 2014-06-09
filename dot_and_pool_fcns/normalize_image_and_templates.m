function [norm_images, norm_templates] = normalize_image_and_templates(image, templates, n_splits)

n_scales = length(templates);
n_templates = size(templates{1},1);
n_ori = size(templates{1},2);
taps = zeros(n_scales,1);
for idx_scale=1:n_scales
    taps(idx_scale) = size(templates{idx_scale},3);
end

sizeX = size(image,2);
sizeY = size(image,1);

norm_images = cell(idx_scales,1);
norm_templates = cell(idx_scales,1);

reg_w = sizeX/n_splits;
reg_h = sizeY/n_splits;
reg_x = ceil( (1:sizeX)/ reg_w);
reg_y = ceil( (1:sizeY)/ reg_h);

for idx_scale=1:n_scales
    
    norm_images{idx_scale} = zeros(sizeY,sizeX);
    norm_templates{idx_scale} = zeros(n_templates, n_ori, taps(idx_scale), taps(idx_scale));
    
    thickness = (taps(idx_scale)-1)/2;

    x_borders_prev = reg_x((thickness+1):end)-reg_x(1:(end-thickness));
    y_borders_prev = reg_y((thickness+1):end)-reg_y(1:(end-thickness));
    x_borders_post = [zeros(1,thickness) -x_borders_prev];
    y_borders_post = [zeros(1,thickness) -y_borders_prev];
    x_borders_prev = [x_borders_prev zeros(1,thickness)];
    y_borders_prev = [y_borders_prev zeros(1,thickness)];
    y_borders_prev = y_borders_prev*n_splits;
    y_borders_post = y_borders_post*n_splits;
    x_borders = x_borders_prev + x_borders_post;
    y_borders = y_borders_prev + y_borders_post;

    [grid_x, grid_y] = meshgrid(reg_x, reg_y);
    [grid_borders_x, grid_borders_y] = meshgrid(x_borders, y_borders);
    grid_borders_crosses = (grid_borders_x+grid_borders_y) .* ((grid_borders_x~=0) & (grid_borders_y~=0));
    
    grid_x=grid_x';
    grid_y=grid_y';
    grid_y = grid_y(:)';
    grid_x = grid_x(:)';
    
    grid_borders_x = grid_borders_x';
    grid_borders_y = grid_borders_y';
    grid_borders_crosses = grid_borders_crosses';
    grid_borders_x = grid_borders_x(:)';
    grid_borders_y = grid_borders_y(:)';
    grid_borders_crosses = grid_borders_crosses(:)';
    
    reg_number = zeros(4,sizeX*sizeY);
    reg_number(1,:) = (grid_y-1)*n_splits + grid_x;
    reg_number(2,:) = reg_number_1 + grid_borders_x;
    reg_number(3,:) = reg_number_1 + grid_borders_y;
    reg_number(4,:) = reg_number_1 + grid_borders_crosses;
    
    n_reg = n_splits^2; % should be equal to the last reg_index
    
    % For each region, store its pixels
    px_in_reg = cell(1,n_reg);
    for idx_reg=1:n_reg
        for idx=1:4
            indices = find(reg_number(idx,:) == idx_reg);
            if isempty(indices),
                continue;
            end
            px_in_reg{idx_reg} = [px_in_reg{idx_reg} indices];
        end
        px_in_reg{idx_reg} = unique(px_in_reg{idx_reg});
    end
    
    image_unrolled = image';
    image_unrolled = image_unrolled(:)';
    
    reg_mean = zeros(n_reg,1);
    reg_norm = zeros(n_reg,1);
    
    for idx_reg=1:n_reg
        reg_mean(idx_reg) = mean(image_unrolled(px_in_reg{idx_reg}));
        reg_norm(idx_reg) = norm(image_unrolled(px_in_reg{idx_reg})-reg_mean(idx_reg));
    end
    
    for ix=1:sizeX
        for iy=1:sizeY
            idx_px = (iy-1)*sizeX + ix;
            its_region = region_number(1,idx_px);
            norm_images{idx_scale}(iy,ix) = (image(iy,ix) - reg_mean(its_region)) / reg_norm(its_region);
        end
    end
    
    templ_mean = zeros(n_templates, n_ori);
    templ_norm = zeros(n_templates, n_ori);

    for idx_template=1:n_templates
        for idx_ori=1:n_ori
            xpad = regW - taps(idx_scale) + thickness;
            ypad = regH - taps(idx_scale) + thickness;
            padded_template = [templates{idx_scale}(idx_template, idx_ori, :, :) zeros(taps(idx_scale),xpad); zeros(ypad, taps(idx_scale)+xpad)];
            
            templ_mean(idx_scale, idx_templ, idx_ori) = mean(padded_template);
            templ_norm(idx_scale, idx_templ, idx_ori) = norm(padded_template - templ_mean(idx_scale, idx_templ, idx_ori));
            norm_templates{idx_scale}(idx_template, idx_ori, :, :) = (templates{idx_scale}(idx_template, idx_ori, :, :) - templ_mean(idx_scale, idx_templ, idx_ori)) / templ_norm(idx_scale, idx_templ, idx_ori);
        end
    end
    
end

end
