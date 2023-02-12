%==========================================================================
% FUNCTION plot_results(p_i,f,nb)
% This function calculates the expected number of infected cells using
% cellular MOI results from simulations and parameters that determine
% virion infectious potential. The input 'p_i' denotes the segment packaging
% probability (value between 0 and 1). The input 'f' denotes the fraction
% of virions that are infection-competent (value between 0 and 1). The
% input 'nb' denotes the burst size.
%
% Color maps for plots are generated using customcolormap.m
% Víctor Martínez-Cagigal (2023). Custom Colormap
% (https://www.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap),
% MATLAB Central File Exchange.
%==========================================================================

function plot_results(p_i,f,nb)

N = 0; Nd = 0;
K = 15;
for k = 1:K
  k
  fname = ['MOI_data_v',num2str(k),'.mat'];
  [n] = infection_prob(p_i,f,nb,fname);
  N = N + n;
end
n = N/K;

load('MOI_data_v1.mat')
[p0,p1] = meshgrid(p0,p1);
n(p0<p1) = nan; nd(p0<p1) = nan; nu(p0<p1) = nan;
surf(n); view(2); shading flat; axis equal;
mycolormap = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
colormap(mycolormap);
axis off; colorbar;

save infection_data.mat n p0 p1 p_i f nb