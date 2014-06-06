function F = design_gabor_filt(f0,taps,B,show)

%
% function design_gabor_filt(f0,taps,B,show)
%
% f0: filter peak frequency, in pixels (default 1/4)
% taps: filter size, in pixels (default 11)
% B: filter bandwidth, in pixels (default f0/3 - 1 octave)
% show: display filters and errors (default 0)
%


if (nargin<1)
  f0 = 1/32;
end
if (nargin<4)
  show = 0;
end
if (nargin<3)
  B = f0/3;
end
if (nargin<2)
  taps = 95;
end

n_orient = 8;


%%%%%%%%%%%%%%%%%%%%
% Check Parameters %
%%%%%%%%%%%%%%%%%%%%


% spatial standard deviation gabor
sigma = sqrt(2*log(2)) ./ (2*pi*B);

% check if number of taps sufficient (larger than 4 times spatial std)
ok = (4.*sigma <= taps);

if (~ok)
   fprintf('More taps required!\n');
end

% check if highest frequency (f0+bandwith) does not exceed Nyquist frequency
ok = (B + f0) <= 1/2;

if (~ok)
   fprintf('Bandwidth too large, exceeding Nyquist frequency!\n'); 
end

% check if bandwith is sufficiently large to cover all orientations
ok = ( n_orient >= pi.*f0./B );

if (~ok)
   fprintf('Bandwidth too small (orientation coverage)!\n');
   fprintf('We recommend at least %d orientations!\n',ceil(pi.*f0./B));
end

% check if bandwith is sufficiently large to cover all scales
ok = (B >= 1/3*f0);

if (~ok)
   fprintf('Bandwidth too small (scale-space coverage)!\n'); 
end


%%%%%%%%%%%%%%%%%%%%%%
% Construct 1D masks %
%%%%%%%%%%%%%%%%%%%%%%


x = -floor(taps/2):floor(taps/2);
sigma = sqrt(2*log(2)) ./ (2*pi*B);
g = 1/sqrt(2*pi*sigma)*exp(-x.^2 ./ (2.*sigma.^2));
sq = sqrt(2)/2;
sq_plus = sqrt(2+sqrt(2))/2;
sq_min = sqrt(2-sqrt(2))/2;

F = zeros(9,taps); % 1D masks (9)

F(1,:) = g;
F(2,:) = g.*cos(2*pi*f0*x);
F(3,:) = g.*sin(2*pi*f0*x);
F(4,:) = g.*cos(2*pi*f0*x*sq);
F(5,:) = g.*sin(2*pi*f0*x*sq);
F(6,:) = g.*cos(2*pi*f0*x*sq_plus);
F(7,:) = g.*sin(2*pi*f0*x*sq_plus);
F(8,:) = g.*cos(2*pi*f0*x*sq_min);
F(9,:) = g.*sin(2*pi*f0*x*sq_min);

% F = F ./ max(abs(F(:))); % largest absolute value = 1

F_noDC = zeros(size(F));
F_noDC_balE = zeros(size(F));


%%%%%%%%%%%%%%%%%%%%%%%
% Remove DC component %
%%%%%%%%%%%%%%%%%%%%%%%

  
acc = 1e-13;

options = optimset('LargeScale','off','Tolcon',acc,'Tolfun',acc, ...
                   'TolX',acc,'Display','final');

if show
  options = optimset(options,'Display','iter');
end

half = ceil(taps/2);

% Horizontal and vertical (0 and 90)

Xt = [ F(1,1:half) F(2,1:half) F(3,1:half) ];

X0 = Xt;

X = fmincon(@opt_gabor_match_HV,X0,[],[],[],[],[],[], ...
            @opt_gabor_DC_HV,options,Xt);

tmp = X(1:half);
F_noDC(1,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(half+1:2*half);
F_noDC(2,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(2*half+1:3*half);
F_noDC(3,:) = [ tmp -tmp(end-1:-1:1) ];

% Diagonal (45 and 135)

Xt = [ F(4,1:half) F(5,1:half) ];

X0 = Xt;

X = fmincon(@opt_gabor_match_D,X0,[],[],[],[],[],[], ...
            @opt_gabor_DC_D,options,Xt);

tmp = X(1:half);
F_noDC(4,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(half+1:2*half);
F_noDC(5,:) = [ tmp -tmp(end-1:-1:1) ];

% 'In-betweens' (22.5 67.5 112.5 157.5)

Xt = [ F(6,1:half) F(7,1:half) F(8,1:half) F(9,1:half) ];

X0 = Xt;

X = fmincon(@opt_gabor_match_IB,X0,[],[],[],[],[],[], ...
            @opt_gabor_DC_IB,options,Xt);

tmp = X(1:half);
F_noDC(6,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(half+1:2*half);
F_noDC(7,:) = [ tmp -tmp(end-1:-1:1) ];
tmp = X(2*half+1:3*half);
F_noDC(8,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(3*half+1:4*half);
F_noDC(9,:) = [ tmp -tmp(end-1:-1:1) ];

% Normalize (max value = 1)

F_noDC = F_noDC ./ max(abs(F_noDC(:)));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove DC component AND balance energy even/odd %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
acc = 1e-13;

options = optimset('LargeScale','off','Tolcon',acc,'Tolfun',acc, ...
                   'TolX',acc,'Display','final');

if show
  options = optimset(options,'Display','iter');
end

half = ceil(taps/2);

% Horizontal and vertical (0 and 90)

Xt = [ F(1,1:half) F(2,1:half) F(3,1:half) ];

X0 = Xt;

X = fmincon(@opt_gabor_match_HV,X0,[],[],[],[],[],[], ...
            @opt_gabor_DC_E_HV,options,Xt);

tmp = X(1:half);
F_noDC_balE(1,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(half+1:2*half);
F_noDC_balE(2,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(2*half+1:3*half);
F_noDC_balE(3,:) = [ tmp -tmp(end-1:-1:1) ];

% Diagonal (45 and 135)

Xt = [ F(4,1:half) F(5,1:half) ];

X0 = Xt;

X = fmincon(@opt_gabor_match_D,X0,[],[],[],[],[],[], ...
            @opt_gabor_DC_E_D,options,Xt);

tmp = X(1:half);
F_noDC_balE(4,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(half+1:2*half);
F_noDC_balE(5,:) = [ tmp -tmp(end-1:-1:1) ];

% 'In-betweens' (22.5 67.5 112.5 157.5)

Xt = [ F(6,1:half) F(7,1:half) F(8,1:half) F(9,1:half) ];

X0 = Xt;

X = fmincon(@opt_gabor_match_IB,X0,[],[],[],[],[],[], ...
            @opt_gabor_DC_E_IB,options,Xt);

tmp = X(1:half);
F_noDC_balE(6,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(half+1:2*half);
F_noDC_balE(7,:) = [ tmp -tmp(end-1:-1:1) ];
tmp = X(2*half+1:3*half);
F_noDC_balE(8,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(3*half+1:4*half);
F_noDC_balE(9,:) = [ tmp -tmp(end-1:-1:1) ];

% Normalize (max value = 1)

% F_noDC_balE = F_noDC_balE ./ max(abs(F_noDC_balE(:)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Round to integer precision % (-255:255)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


F_int_noDC = integer_rounding(F_noDC,taps);
F_int_noDC_balE = integer_rounding(F_noDC_balE,taps);
  

if show


  %%%%%%%%%%%%%%%%%%%
  % Display Results %
  %%%%%%%%%%%%%%%%%%%


  % Construct 2D masks
  
  [E,O] = compose_gabor_filters(F);
  [E_noDC,O_noDC] = compose_gabor_filters(F_noDC);

  % Compute some errors
  
  MSE_E = mean(reshape((E - E_noDC).^2,[taps^2 8]));
  MSE_O = mean(reshape((O - O_noDC).^2,[taps^2 8]));

  DC_E = squeeze(sum(sum(E)));
  DC_O = squeeze(sum(sum(O)));
  
  DC_E_noDC = squeeze(sum(sum(E_noDC)));
  DC_O_noDC = squeeze(sum(sum(O_noDC)));
  
  % Show

  cl = [ -1 1 ];
  
  angles = 0:pi/n_orient:pi-(pi/n_orient)  

  for a = 1:n_orient

    subplot(2,2,1);
    plgr(E(:,:,a));
    title('E');
    colormap gray
    set(gca,'clim',cl);
    colorbar

    subplot(2,2,2);
    plgr(O(:,:,a));
    title('O');
    colormap gray
    set(gca,'clim',cl);
    colorbar

    subplot(2,2,3);
    plgr(E_noDC(:,:,a));
    title('E_noDC');
    colormap gray
    set(gca,'clim',cl);
    colorbar

    subplot(2,2,4);
    plgr(O_noDC(:,:,a));
    title('o_noDC');
    colormap gray
    set(gca,'clim',cl);
    colorbar

    fprintf('angle: %3.1f\n',(a*180/pi));
    fprintf('MSE E : %2.1e  O : %2.1e\n',MSE_E(a),MSE_O(a));
    fprintf('DC  E : %2.1e  O : %2.1e  E_noDC : %2.1e  O_noDC : %2.1e\n', ...
            DC_E(a),DC_O(a),DC_E_noDC(a),DC_O_noDC(a));
    
    pause;
    
    fprintf('\n\n');
    
  end
  
end


%%%%%%%%%%%%%%%%%%
% Save Filter(s) %
%%%%%%%%%%%%%%%%%%


type = 'g';

basename = [ 'Gt' int2str(taps) 'B' num2str(B,'%1.4f') ...
             'f' num2str(f0,'%1.3f')];

% energy unbalanced and dc remaining

eval(['save ' basename 'DE.mat F taps type f0 B;']);

% energy unbalanced and dc removed

F = F_noDC;
eval(['save ' basename 'E.mat F taps type f0 B;']);

% energy balanced and dc removed

F = F_noDC_balE;
eval(['save '  basename '.mat F taps type f0 B;']);

% integer, energy unbalanced and dc removed

F = F_int_noDC;
eval(['save ' basename 'IE.mat F taps type f0 B;']);

% integer, energy balanced and dc removed

F = F_int_noDC_balE;
eval(['save ' basename 'I.mat F taps type f0 B;']);














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F_int = integer_rounding(F,taps)
  
% Round the filter to integer precision and enforce zero DC (no energy
% balancing)

half = ceil(taps/2);

F_int = 255.*F;
F_intr = round(F_int);

% only 1D filters 2,4 and 6 must sum to zero; the others are odd and filter
% number 8 is always used in combination with 6 or an odd filter

check = [ 2 4 6 ];

for f = check
  
  DC = sum(F_intr(f,:));
  
  switch abs(DC)
    
   case 0
    
    % nothing to do
    
   case 1
    
    % Add or subtract one from the center value
    
    F_intr(f,half) = F_intr(f,half) - DC;
    
   case 2
    
    % Add or subtract one symmetrically from the value with the largest
    % rounding error
    
    [dummy,I] = max(abs(F_intr(f,1:half-1) - F_int(f,1:half-1)));
    
    F_intr(f,I) = F_intr(f,I) - sign(DC);
    F_intr(f,taps-I+1) = F_intr(f,taps-I+1) - sign(DC);

   case 3
    
    % Apply both modifications from above

    F_intr(f,half) = F_intr(f,half) - sign(DC);

    [dummy,I] = max(abs(F_intr(f,1:half-1) - F_int(f,1:half-1)));
    
    F_intr(f,I) = F_intr(f,I) - sign(DC);
    F_intr(f,taps-I+1) = F_intr(f,taps-I+1) - sign(DC);
    
   case 4
    
    % Add or subtract one symmetrically from the values with the largest
    % rounding errors

    [dummy,I] = sort(abs(F_intr(f,1:half-1) - F_int(f,1:half-1)), ...
                     'descend');
    
    for i = 1:2
      F_intr(f,I(i)) = F_intr(f,I(i)) - sign(DC);
      F_intr(f,taps-I(i)+1) = F_intr(f,taps-I(i)+1) - sign(DC);
    end
    
   otherwise
    
    fprintf('\n\nROUNDING PROBLEM !!!\n');
    
  end
  
end

F_int = F_intr;
