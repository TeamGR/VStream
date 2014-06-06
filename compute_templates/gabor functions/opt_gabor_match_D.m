function f = opt_gabor_match_D(X,Xt)
  
  
half = length(X)./2;
taps = half + half-1;


%%%%%%%%%%%%%%%%%%%%%%
% Extract 1D filters % (enforce symmetry or anti-symmetry)
%%%%%%%%%%%%%%%%%%%%%%


% Current estimates

F_noDC = zeros(9,taps); % only 2 are used

tmp = X(1:half);
F_noDC(4,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(half+1:2*half);
F_noDC(5,:) = [ tmp -tmp(end-1:-1:1) ];

% Theoretical

F = zeros(9,taps); % only 2 are used

tmp = Xt(1:half);
F(4,:) = [ tmp tmp(end-1:-1:1) ];
tmp = Xt(half+1:2*half);
F(5,:) = [ tmp -tmp(end-1:-1:1) ];


%%%%%%%%%%%%%%%%%%%%%%
% Compose 2D filters %
%%%%%%%%%%%%%%%%%%%%%%


% Current estimates

G_noDC = zeros(taps,taps,4);

G_noDC(:,:,1) = F_noDC(4,:).'*F_noDC(4,:) - F_noDC(5,:).'*F_noDC(5,:);
G_noDC(:,:,2) = F_noDC(5,:).'*F_noDC(4,:) + F_noDC(4,:).'*F_noDC(5,:);

% OPIMIZATION OF THE SYMMETRICAL FILTERS IS NOT NECESSARY AND CONFUSES
% THE OPTIMIZATION ALGORITHM !!!

%G_noDC(:,:,3) = F_noDC(4,:).'*F_noDC(4,:) + F_noDC(5,:).'*F_noDC(5,:);
%G_noDC(:,:,4) = F_noDC(5,:).'*F_noDC(4,:) - F_noDC(4,:).'*F_noDC(5,:);

% Theoretical

G = zeros(taps,taps,4);

G(:,:,1) = F(4,:).'*F(4,:) - F(5,:).'*F(5,:); % pi/4 even
G(:,:,2) = F(5,:).'*F(4,:) + F(4,:).'*F(5,:); % pi/4 odd
%G(:,:,3) = F(4,:).'*F(4,:) + F(5,:).'*F(5,:); % 3pi/4 even
%G(:,:,4) = F(5,:).'*F(4,:) - F(4,:).'*F(5,:); % 3pi/4 odd


%%%%%%%%%%%%%%%%%%%%%%
% Mean Squared Error %
%%%%%%%%%%%%%%%%%%%%%%


f = mean((G_noDC(:) - G(:)).^2);
