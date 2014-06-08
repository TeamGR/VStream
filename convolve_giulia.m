function responses = convolve_giulia(image, templates)

n_scales = length(templates);
n_templates = size(templates{1},1);
n_ori = size(templates{1},2);

taps = zeros(n_scales,1);
for idx_scale=1:n_scales
    taps(idx_scale) = size(templates{idx_scale},3);
end

responses = cell(n_scales, 1);
for idx_scale=1:n_scales
    responses{idx_scale,1} = zeros(n_templates, n_ori, inSizeY-(taps(idx_scale)-1), inSizeX-(taps(idx_scale)-1));
    for idx_template=1:n_templates
        for idx_ori=1:n_ori                     
            responses{idx_scale}(idx_template, idx_ori, :, :) = conv2(image,squeeze(templates{idx_scale,1}(idx_template, idx_ori, :, :)),'valid');
        end
    end
end

end