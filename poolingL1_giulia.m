function out = poolingL1_giulia(responses, n_splits, n_bins, range, policy)

%   responses: reponses from previous layer of type cell(n_scales,1) 
%               with each elem of format zeros(n_templates*n_ori*sizeY*sizeX)
%   n_bins : bumber of bars of the histogram
%   range : range of the histogram
%   policy : pooling policy
%               'histogram': provides the concatenated histogram signature
%               'moments': provides the first 3 moments of the distribution
%               (mean, variance and max value)
%   out: either    out_hist = zeros(n_bins, n_reg, n_templates, 2)
%               or out_moms = zeros(3, n_reg, n_templates)

% Extract info from responses
n_scales = length(responses);
n_templates = size(responses{1},1);
n_ori = size(responses{1},2);

% Pick the size of the responses at the greatest scale
sizeX = size(responses{n_scales},4);
sizeY = size(responses{n_scales},3);

% Crop borders so that all responses are as small as those at the greatest scale
for idx_scale=1:(n_scales-1)
    
    border_cut = (size(responses{idx_scale},4)-sizeX)/2;
    responses{idx_scale}(:, :, 1:border_cut, :) = [];
    responses{idx_scale}(:, :, (end-border_cut+1):end, :) = [];
    responses{idx_scale}(:, :, :, 1:border_cut) = [];
    responses{idx_scale}(:, :, :, (end-border_cut+1):end) = [];
end

% For each pixel, store its region index
reg_w = sizeX/n_splits;
reg_h = sizeY/n_splits;
reg_x = ceil( (1:sizeX)/ reg_w);
reg_y = ceil( (1:sizeY)/ reg_h);
[grid_x, grid_y] = meshgrid(reg_x, reg_y);
grid_x=grid_x';
grid_y=grid_y';
grid = [grid_y(:)'; grid_x(:)'];
reg_number = (grid(1,:)-1)*n_splits + grid(2,:);
n_reg = n_splits^2; % should be equal to the last reg_index

% For each region, store its pixels
px_in_reg = cell(1,n_reg);
for idx_reg=1:n_reg
    indices = find(reg_number == idx_reg);
    if isempty(indices),
        continue;
    end
    px_in_reg{idx_reg} = indices;
end

% Unroll the response matrices
responses_unrolled = zeros(n_templates, sizeX*sizeY, n_ori, n_scales);
for idx_template=1:n_templates
    for idx_scale=1:n_scales
        for idx_ori=1:n_ori        
             resp_matrix = squeeze(responses{idx_scale}(idx_template, idx_ori, :, :))'; 
             responses_unrolled(idx_template, :, idx_ori, idx_scale) = resp_matrix(:); 
        end
    end
end

% Group along scale and orientation and separate the regions using 'px_in_reg'
domain_hist = cell(n_templates, n_reg);
for idx_template=1:n_templates
    for idx_reg=1:n_reg
         domain_single_hist = squeeze(responses_unrolled(idx_template, px_in_reg{idx_reg}, :, :));
         domain_hist{idx_template, idx_reg} = domain_single_hist(:);
    end
end
 
% Compute the histogram or the moments
switch policy
    case 'histogram'
        
        min_hist = range(1);
        max_hist = range(2);
        x_hist = linspace(min_hist,max_hist,n_bins);
        
        out_hist = zeros(n_bins, n_reg, n_templates, 2);
        
        for idx_template=1:n_templates
            for idx_reg = 1:n_reg
                [out_hist(:, idx_reg, idx_template, 1), out_hist(:, idx_reg, idx_template, 2)] = hist(domain_hist{idx_template, idx_reg}, x_hist);
            end
        end
        
        out = out_hist;
         
    case 'moments'
        
        out_moms = zeros(3, n_reg, n_templates);
        
        for idx_template=1:n_templates
            for idx_reg = 1:n_reg
                out_moms(1, idx_reg, idx_template) = mean(domain_hist{idx_template, idx_reg});
                out_moms(2, idx_reg, idx_template) = var(domain_hist{idx_template, idx_reg});
                out_moms(3, idx_reg, idx_template) = max(domain_hist{idx_template, idx_reg});
            end
        end
     
        out = out_moms;
end

end

