function plot_filter(E,O)
% close all
Real = E;
Imag = O;
% Real = G{1,1};
% Imag = G{1,2};

n_orient = size(Real,3);

Ecomp(1,:)=squeeze(sum(sum(Real.^2)));
Ecomp(2,:)=squeeze(sum(sum(Imag.^2)));
DCcomp(1,:)=squeeze(sum(sum(Real)));
DCcomp(2,:)=squeeze(sum(sum(Imag)));

figure

for j=0:n_orient-1
% j=0;
%         Real2(:,:,j+1)=expand(Real(:,:,j+1));  
%         Real2=Real;
%         [sy,sx,ori]=size(Real2);
%         plot(Real2(:,:,j+1)')
        subplot(n_orient/2,2,j+1)
        mesh(Real(:,:,j+1)), axis square, colormap(gray)
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        xlabel(['sum =' num2str(Ecomp(1,j+1)) 'sum =' num2str(DCcomp(1,j+1))]);
end
text(-120,-160,'Filter real part')
text(-100,25,'Phase')
text(-220,-60,'Orient')

figure

for j=0:n_orient-1
% j=0;
%         Imag2(:,:,j+1)=expand(Imag(:,:,j+1));  
%         Imag2=Imag;
%         [sy,sx,ori]=size(Imag2);
%         plot(Imag2(:,:,j+1)')        
        subplot(n_orient/2,2,j+1)
        mesh(Imag(:,:,j+1)), axis square, colormap(gray)
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        xlabel(['sum=' num2str(Ecomp(2,j+1)) 'sum=' num2str(DCcomp(2,j+1))]);
end
text(-120,-160,'Filter imaginary part')
text(-100,25,'Phase')
text(-220,-60,'Orient')