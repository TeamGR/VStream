function [c,ceq] = opt_gabor_DC_E_D(X,Xt)
  
  
c = 0;

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


%%%%%%%%%%%%%%%%%%%%%%
% Compose 2D filters %
%%%%%%%%%%%%%%%%%%%%%%


% Current estimates

G_noDC = zeros(taps,taps,4);

G_noDC(:,:,1) = F_noDC(4,:).'*F_noDC(4,:) - F_noDC(5,:).'*F_noDC(5,:);
G_noDC(:,:,2) = F_noDC(5,:).'*F_noDC(4,:) + F_noDC(4,:).'*F_noDC(5,:);

% OPTIMIZATION OF THE SYMMETRICAL FILTERS IS NOT NECESSARY AND CONFUSES
% THE OPTIMIZATION ALGORITHM !!!

%G_noDC(:,:,3) = F_noDC(4,:).'*F_noDC(4,:) + F_noDC(5,:).'*F_noDC(5,:);
%G_noDC(:,:,4) = F_noDC(5,:).'*F_noDC(4,:) - F_noDC(4,:).'*F_noDC(5,:);


%%%%%%%%%%%%%%%%
% DC component %
%%%%%%%%%%%%%%%%


tmp = sum(sum(G_noDC,1),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy balancing (even energy - odd energy) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tmp2 = sum(sum(G_noDC.^2,1),2);
tmp2 = [ tmp2(1) - tmp2(2)
         tmp2(3) - tmp2(4) ];


%%%%%%%%%%%%%%%%%%%%%%%
% Equality constraint %
%%%%%%%%%%%%%%%%%%%%%%%


ceq = [ tmp(:) ; tmp2(:) ];
