function f = opt_gabor_match_IB(X,Xt)
  
  
half = length(X)./4;
taps = half + half-1;


%%%%%%%%%%%%%%%%%%%%%%
% Extract 1D filters % (enforce symmetry or anti-symmetry)
%%%%%%%%%%%%%%%%%%%%%%


% Current estimates

F_noDC = zeros(9,taps); % only 4 are used

tmp = X(1:half);
F_noDC(6,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(half+1:2*half);
F_noDC(7,:) = [ tmp -tmp(end-1:-1:1) ];
tmp = X(2*half+1:3*half);
F_noDC(8,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(3*half+1:4*half);
F_noDC(9,:) = [ tmp -tmp(end-1:-1:1) ];

% Theoretical

F = zeros(9,taps); % only 4 are used

tmp = Xt(1:half);
F(6,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(half+1:2*half);
F(7,:) = [ tmp -tmp(end-1:-1:1) ];
tmp = X(2*half+1:3*half);
F(8,:) = [ tmp tmp(end-1:-1:1) ];
tmp = X(3*half+1:4*half);
F(9,:) = [ tmp -tmp(end-1:-1:1) ];


%%%%%%%%%%%%%%%%%%%%%%
% Compose 2D filters %
%%%%%%%%%%%%%%%%%%%%%%


% Current estimates

G_noDC = zeros(taps,taps,8);

G_noDC(:,:,1) = F_noDC(8,:).' * F_noDC(6,:) - (F_noDC(9,:).' * F_noDC(7,:));
G_noDC(:,:,2) = F_noDC(9,:).' * F_noDC(6,:) + (F_noDC(8,:).' * F_noDC(7,:));
G_noDC(:,:,3) = F_noDC(6,:).' * F_noDC(8,:) - (F_noDC(7,:).' * F_noDC(9,:));
G_noDC(:,:,4) = F_noDC(7,:).' * F_noDC(8,:) + (F_noDC(6,:).' * F_noDC(9,:));
G_noDC(:,:,5) = F_noDC(6,:).' * F_noDC(8,:) + (F_noDC(7,:).' * F_noDC(9,:));
G_noDC(:,:,6) = F_noDC(7,:).' * F_noDC(8,:) - (F_noDC(6,:).' * F_noDC(9,:));
G_noDC(:,:,7) = F_noDC(8,:).' * F_noDC(6,:) + (F_noDC(9,:).' * F_noDC(7,:));
G_noDC(:,:,8) = F_noDC(9,:).' * F_noDC(6,:) - (F_noDC(8,:).' * F_noDC(7,:));

% Theoretical

G = zeros(taps,taps,8);

G(:,:,1) = F(8,:).' * F(6,:) - (F(9,:).' * F(7,:)); % pi/8 even
G(:,:,2) = F(9,:).' * F(6,:) + (F(8,:).' * F(7,:)); % pi/8 odd
G(:,:,3) = F(6,:).' * F(8,:) - (F(7,:).' * F(9,:)); % 3pi/8 even
G(:,:,4) = F(7,:).' * F(8,:) + (F(6,:).' * F(9,:)); % 3pi/8 odd
G(:,:,5) = F(6,:).' * F(8,:) + (F(7,:).' * F(9,:)); % 5pi/8 even
G(:,:,6) = F(7,:).' * F(8,:) - (F(6,:).' * F(9,:)); % 5pi/8 odd
G(:,:,7) = F(8,:).' * F(6,:) + (F(9,:).' * F(7,:)); % 7pi/8 even
G(:,:,8) = F(9,:).' * F(6,:) - (F(8,:).' * F(7,:)); % 7pi/8 odd


%%%%%%%%%%%%%%%%%%%%%%
% Mean Squared Error %
%%%%%%%%%%%%%%%%%%%%%%


f = mean((G_noDC(:) - G(:)).^2);
