function f =  opt_gabor_match_HV(X,Xt)
  
  
half = length(X)./3;
taps = half + half-1;


%%%%%%%%%%%%%%%%%%%%%%
% Extract 1D filters % (enforce symmetry or anti-symmetry)
%%%%%%%%%%%%%%%%%%%%%%


% Current estimates

F_noDC = zeros(9,taps); % only 3 are used

tmp = X(1:half);
F_noDC(1,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(half+1:2*half);
F_noDC(2,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(2*half+1:3*half);
F_noDC(3,:) = [ tmp -tmp(end-1:-1:1) ];

% Theoretical

F = zeros(9,taps); % only 3 are used

tmp = Xt(1:half);
F(1,:) = [ tmp tmp(end-1:-1:1) ];
tmp = Xt(half+1:2*half);
F(2,:) = [ tmp tmp(end-1:-1:1) ];
tmp = Xt(2*half+1:3*half);
F(3,:) = [ tmp -tmp(end-1:-1:1) ];


%%%%%%%%%%%%%%%%%%%%%%
% Compose 2D filters %
%%%%%%%%%%%%%%%%%%%%%%


% Current estimates

G_noDC = zeros(taps,taps,4);

G_noDC(:,:,1) = F_noDC(1,:).'*F_noDC(2,:); % 0 even
G_noDC(:,:,2) = F_noDC(1,:).'*F_noDC(3,:); % 0 odd
G_noDC(:,:,3) = F_noDC(2,:).'*F_noDC(1,:); % pi/2 even
G_noDC(:,:,4) = F_noDC(3,:).'*F_noDC(1,:); % pi/2 odd

% Theoretical

G = zeros(taps,taps,4);

G(:,:,1) = F(1,:).'*F(2,:);
G(:,:,2) = F(1,:).'*F(3,:);
G(:,:,3) = F(2,:).'*F(1,:);
G(:,:,4) = F(3,:).'*F(1,:);


%%%%%%%%%%%%%%%%%%%%%%
% Mean Squared Error %
%%%%%%%%%%%%%%%%%%%%%%


f = mean((G_noDC(:) - G(:)).^2);
