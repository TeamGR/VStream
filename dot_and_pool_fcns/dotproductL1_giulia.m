function responses = dotproductL1_giulia(image, templates, n_splits)

n_scales = length(templates);
n_templates = size(templates{1},1);
n_ori = size(templates{1},2);

taps = zeros(n_scales,1);
for idx_scale=1:n_scales
    taps(idx_scale) = size(templates{idx_scale},3);
end

sizeX = size(image,2);
sizeY = size(image,1);

[norm_images, norm_templates] = normalize_image_and_templates(image, templates, n_splits);

responses = cell(n_scales, 1);
for idx_scale=1:n_scales
    responses{idx_scale,1} = zeros(n_templates, n_ori, sizeY-(taps(idx_scale)-1), sizeX-(taps(idx_scale)-1));
    for idx_template=1:n_templates
        for idx_ori=1:n_ori                     
            responses{idx_scale}(idx_template, idx_ori, :, :) = conv2(norm_images{idx_scale},squeeze(norm_templates{idx_scale}(idx_template, idx_ori, :, :)),'valid');
        end
    end
end

end