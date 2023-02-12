%==========================================================================
% FUNCTION [S,N_disperse] = cell_moi(P,a,N);
% This function calculates cellular MOI values on a ~25x25 hexagonal
% lattice of cells given simulation results from the function
% 'lattice_spread'. 
% 
% INPUTS: 'P' gives the positions of bound virions; 'a' is the
% center-to-center distance of cells in the lattice; 'N' is the number of
% virions in the simulation.
%
% OUTPUTS: 'S' gives the virion count for naive cells in the lattice, with
% each element representing a cell (the parent cell is excluded).
% 'N_disperse' gives the number of virions that attach outside the bounds
% of the hexagonal lattice.
%==========================================================================

function [S,N_disperse] = cell_moi(P,a,N);

close all;

C = [];
CN = 12;        % Define the size of the hexagonal lattice: ~(2*CN+1) x (2*CN+1)
for Cy = -CN:CN
  if(mod(Cy,2)==0)
    temp(:,1) = a*[-CN:CN];
  else
    temp(:,1) = a*[-(CN+0.5):(CN+0.5)];
  end
  temp(:,2) = a*sqrt(3)/2*Cy*ones(size(temp(:,1)));
  C = [C ; temp];
  clear temp;
end

theta = [30:60:330]*pi/180; s = a/sqrt(3);
for k = 1:length(C);
  H(:,1) = cos(theta')*s + C(k,1); H(:,2) = sin(theta')*s + C(k,2);
  in = inpolygon(P(:,1),P(:,2),H(:,1),H(:,2));
  S(k) = sum(in);
end

% Determine the number of particles bound outside the lattice sites:
N_disperse = size(P,1)-sum(sum(S));

% Remove data for virions that remain attached to the parent cell:
i = (length(S)+1)/2;
S(i) = [];