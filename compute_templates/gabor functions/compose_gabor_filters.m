function [E,O] = compose_gabor_filters(F)


taps = size(F,2);
E = zeros(taps,taps,8);
O = zeros(taps,taps,8);


% 0

E(:,:,1) = F(1,:).'*F(2,:);
O(:,:,1) = F(1,:).'*F(3,:);

% pi/8

E(:,:,2) = F(8,:).' * F(6,:) - (F(9,:).' * F(7,:));
O(:,:,2) = F(9,:).' * F(6,:) + (F(8,:).' * F(7,:));

% pi/4

E(:,:,3) = F(4,:).'*F(4,:) - F(5,:).'*F(5,:);
O(:,:,3) = F(5,:).'*F(4,:) + F(4,:).'*F(5,:);

% 3pi/8

E(:,:,4) = F(6,:).' * F(8,:) - (F(7,:).' * F(9,:));
O(:,:,4) = F(7,:).' * F(8,:) + (F(6,:).' * F(9,:));

% pi/2

E(:,:,5) = F(2,:).'*F(1,:);
O(:,:,5) = F(3,:).'*F(1,:);

% 5pi/8

E(:,:,6) = F(6,:).' * F(8,:) + (F(7,:).' * F(9,:));
O(:,:,6) = F(7,:).' * F(8,:) - (F(6,:).' * F(9,:));

% 3pi/4

E(:,:,7) = F(4,:).'*F(4,:) + F(5,:).'*F(5,:);
O(:,:,7) = F(5,:).'*F(4,:) - F(4,:).'*F(5,:);

% 7pi/8

E(:,:,8) = F(8,:).' * F(6,:) + (F(9,:).' * F(7,:));
O(:,:,8) = F(9,:).' * F(6,:) - (F(8,:).' * F(7,:));
