%==========================================================================
% FUNCTION main
% This function simulates virion diffusion, attachment, and infection for a
% range of cis and trans binding probabilities. Results are saved to a
% .mnat file (e.g. 'MOI_data_v1.mat','MOI_data_v2.mat', etc. for each
% simulation).
%==========================================================================

function main

N = 10000;          % Number of viral particles simulated;
T = 8*60*60;        % Total time period simulated [s]

p0 = logspace(-2,0,30);     % Cis binding probabilities;
p1 = logspace(-2,0,30);     % Trans binding probabilities;

K = 15;             % Run the simulations K times (larger values reduce stochastic noise)
for k = 1:K        
  for n0 = 1:length(p0)
    for n1 = 1:length(p1);
      n1
      if(p1(n1)>=p0(n0));
        p2 = p1(n1);
        [S{n0,n1},N_disperse(n0,n1)] = lattice_spread(p0(n0),p1(n1),p2,N,T);
      end
    end
  end
  sname = ['MOI_data_v',num2str(k),'.mat'];
  save(sname,'S','N_disperse','N','p0','p1','T');
end