function out = poolingL2_giulia(responses, n_splits, n_bins, range, policy)

signature_length = size(responses,1);
n_ori = size(responses,2);
n_scales = size(responses,3);
n_templates = size(responses,4);

% Compute the histogram or the moments
switch policy
    case 'histogram'
        
        min_hist = range(1);
        max_hist = range(2);
        x_hist = linspace(min_hist,max_hist,n_bins);
        
        out_hist = zeros(n_bins, n_templates, 2);
        
        for idx_template=1:n_templates 
            domain_hist = responses(:, :, :, idx_template);
            [out_hist(:, idx_template, 1), out_hist(:, idx_template, 2)] = hist(domain_hist(:), x_hist);
        end
        
        out = out_hist;
         
    case 'moments'
        
        out_moms = zeros(3, n_templates);
        
        for idx_template=1:n_templates
            domain_hist = responses(:, :, :, idx_template);
            out_moms(1, idx_template) = mean(domain_hist(:));
            out_moms(2, idx_template) = var(domain_hist(:));
            out_moms(3, idx_template) = max(domain_hist(:));
        end
     
        out = out_moms;
end

end