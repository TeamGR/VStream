function out = dotproduct_giulia(response, templates)

signature_length = length(response); % should be n_bins*n_reg*n_templates

n_templates = size(templates,1);
n_ori = size(templates,2);
n_scales = size(templates,3);

out = zeros(signature_length, n_ori, n_scales, n_templates);

for idx_template=1:n_templates
    for idx_ori=1:n_ori
        for idx_scale=1:n_scales  
            out(:, idx_ori, idx_scale, idx_template) = response.*templates{idx_template, idx_ori, idx_scale};
        end
    end
end

end