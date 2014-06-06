function [c,ceq] = opt_gabor_DC_E_IB(X,Xt)
  
  
c = 0;

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


%%%%%%%%%%%%%%%%
% DC component %
%%%%%%%%%%%%%%%%


tmp = sum(sum(G_noDC,1),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy balancing (even energy - odd energy) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tmp2 = sum(sum(G_noDC.^2,1),2);
tmp2 = [ tmp2(1) - tmp2(2)
         tmp2(3) - tmp2(4) 
         tmp2(5) - tmp2(6)
         tmp2(7) - tmp2(8) ];


%%%%%%%%%%%%%%%%%%%%%%%
% Equality constraint %
%%%%%%%%%%%%%%%%%%%%%%%


ceq = [ tmp(:) ; tmp2(:) ];
