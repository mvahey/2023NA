%==========================================================================
% [N] = infection_prob(p_i,f,fname);
% This function determines the expected vlue for the number of infected
% cells ('N'), based on the results of simulations saved in the file specified by
% the variable 'fname', the segment packaging probabilities ('p_i'), and the
% fraction of virions that are infection-competent ('f'). The viral burst
% size is given by the input 'nb'.
%==========================================================================

function [N] = infection_prob(p_i,f,nb,fname);

s = p_i*ones(1,8);    % Segment packaging probabilities

load(fname)
np = N;   % Number of particles in the simulation

for n0 = 1:length(p0)
  for n1 = 1:length(p1);
    if(p1(n1)>=p0(n0));
      N = round(S{n0,n1}*nb/np);
      for q = 1:length(N);
        P0(q) = 0;
        for k = 1:N(q)
          pk = binopdf(k,N(q),f);   % Probability that, of N(q) particles, k are 'competent'
          P0(q) = P0(q) + pk*prod(1 - ((1-s).^(k)));
        end
      end
      Ni(n0,n1) = sum(P0);
    end
  end
end

% Approximate the contributions of the dispersed virions:
Nd = (N_disperse*nb/np)*prod(s)*f;

N = Ni + Nd;