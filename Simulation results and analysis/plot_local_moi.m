%==========================================================================
% FUNCTION plot_local_moi(nb)
% This function calculates the number of virions per cell that attach to
% the six nearest neighbors of the infected parent cell using a hexagonal
% lattice. The total number of viral particles is scaled to match the input
% burst size, 'nb'.
%
% Color maps for plots are generated using customcolormap.m
% Víctor Martínez-Cagigal (2023). Custom Colormap
% (https://www.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap),
% MATLAB Central File Exchange.
%==========================================================================

function plot_local_moi(nb)

K = 15;
for k = 1:K
  fname = ['MOI_data_v',num2str(k),'.mat'];
  load(fname);
  for n0 = 1:length(p0);
    for n1 = 1:length(p1);
      if(p1(n1)>=p0(n0));
        temp = sort(S{n0,n1}*nb/N,'descend');
        s(n0,n1,k) = mean(temp(1:6));
      end
    end
  end
end

s0 = mean(s,3);
 
[p0,p1] = meshgrid(p0,p1);
s0(p0<p1) = nan;
surf(s0); view(2); shading flat; axis equal;
mycolormap = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
colormap(mycolormap);
axis off; colorbar;