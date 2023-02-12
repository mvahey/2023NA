%==========================================================================
% FUNCTION [S,Nd] = lattice_spread(p0,p1,p2,N,T)
% This function determines the cellular MOI for cells in a hexagonal
% lattice. 
%
% INPUTS: 'p0', 'p1', 'p2' indicate binding probability for virions encountering the
% infected cell (p0), nearest neighbors (p1), or distant cells (p2). Input
% 'N' is the number of particles to simulate. 'T' indicates the total
% duraction of the simulation (in seconds).
%
% OUTPUTS: Each element in the output 'S' gives the virion count for a
% naive cell within the lattice.  The output 'Nd' gives the number of
% virions that attach outside the bounds of the hexagonal lattice.
%==========================================================================

function [S,Nd] = lattice_spread(p0,p1,p2,N,T)

a = 30;             % Distance between centers of neighboring cells [um]
s = a/sqrt(3);      % Length of hexagonal cell sides [um]
d = 3;              % Diffusive step size [um]
dt = 1;             % Time step [s]

% Define a central hexagon for the infected cell:
C = [0,0 ; a,0 ; a/2,sqrt(3)*a/2 ; -a/2,sqrt(3)*a/2 ; -a,0 ; -a/2,-sqrt(3)*a/2 ; a/2,-sqrt(3)*a/2];
theta = [30:60:330]*pi/180;
H0(:,1) = cos(theta)*s + C(1,1);
H0(:,2) = sin(theta)*s + C(1,2);

% Define hexagon boundaries for neighboring cells:
H1(:,1) = cos(theta')*s + C(2,1); H1(:,2) = sin(theta')*s + C(2,2);
H2(:,1) = cos(theta')*s + C(3,1); H2(:,2) = sin(theta')*s + C(3,2);
H3(:,1) = cos(theta')*s + C(4,1); H3(:,2) = sin(theta')*s + C(4,2);
H4(:,1) = cos(theta')*s + C(5,1); H4(:,2) = sin(theta')*s + C(5,2);
H5(:,1) = cos(theta')*s + C(6,1); H5(:,2) = sin(theta')*s + C(6,2);
H6(:,1) = cos(theta')*s + C(7,1); H6(:,2) = sin(theta')*s + C(7,2);
H = [];
H = [H ; H2(6,:) ; H2(1,:) ; H2(2,:) ; H2(3,:)];
H = [H ; H3(2,:) ; H3(3,:) ; H3(4,:)];
H = [H ; H4(3,:) ; H4(4,:) ; H4(5,:)];
H = [H ; H5(4,:) ; H5(5,:) ; H5(6,:)];
H = [H ; H6(5,:) ; H6(6,:) ; H6(1,:)];
H = [H ; H1(6,:) ; H1(1,:)];

% Initialize virion positions:
P = a*(rand(2*N,2)-1/2);
in = inpolygon(P(:,1),P(:,2),H0(:,1),H0(:,2));
P = P(in,:); P = P(1:N,:); P(:,3) = P(:,1)*0;

P_stuck = [];
for t = 0:dt:T;
    if(t==0)
      dP = 0*P;
    else
      dP = d*randn(size(P));
    end
  P1 = P + dP;
  in0 = inpolygon(P1(:,1),P1(:,2),H0(:,1),H0(:,2));
  in1 = inpolygon(P1(:,1),P1(:,2),H(:,1),H(:,2));
  i0 = find((P1(:,3)<0).*in0);            % Find the particles that contact the infected cell;
  i1 = find((P1(:,3)<0).*(in1-in0));      % Find the particles that contact the neighbor cells;
  i2 = find((P1(:,3)<0).*(1-in1));        % Find the particles that contact far away;
  r0 = rand(size(i0))<p0;
  r1 = rand(size(i1))<p1;
  r2 = rand(size(i2))<p2;
  i_stuck = [i0(r0) ; i1(r1) ; i2(r2)];   % Assign particles as stuck based on the probabilities p0, p1, p2;
  P1(i_stuck,3) = 0;
  i0 = find(P1(:,3)<0); P1(i0,3) = -P1(i0,3); P = P1;
  P_stuck = [P_stuck ; P(i_stuck,:)];
  P(i_stuck,:) = [];

end

[S,Nd] = cell_moi(P_stuck,a,N); 