%% 1D heat conduction Module
% Input: 
%   Geometry Tau
%   Spatial discretization J (or mesh size)
%   heat source pTriplePrime_j
%   Manufactured source p_MMS_j
%   Boundary condition
% Output: 
%   Cell-averaged temperature
function [T_j]=heat_cond_module(J,Tau,T_L,T_R,pTriplePrime_j,p_MMS_j)
  %% Example also as optional param
  if ~exist('J','var')
    J=5;%*2;%*2*2*2*2*2*2*2*2
  end
  if ~exist('Tau','var')
    Tau=10;
  end
  if ~exist('T_L','var')
    T_L=0;
  end
  if ~exist('T_R','var')
    T_R=100;
  end
  if ~exist('pTriplePrime','var')
    pTriplePrime_j=ones(J,1)*0.2; %kappa=1, sigma_f=0.1; phi=2.0
  end
  if ~exist('p_MMS_j','var')
    p_MMS_j=ones(J,1)*2.2;
  end

%% Solver
  delta_z_j=ones(J,1)*Tau/J;
  delta_z=Tau/J; % Assuming uniform grid

  A=zeros(J,J);
  % The way A is set up assumes uniform grid.
  A(1,1)=-3;
  A(1,2)=1;
  for j=2:J-1
    A(j,j-1)=1;
    A(j,j)=-2;
    A(j,j+1)=1;
  end
  A(J,J-1)=1;
  A(J,J)=-3;
  A=A/(delta_z*delta_z);

  f=zeros(J,1);
  for j=1:J
    f(j)=-pTriplePrime_j(j)+p_MMS_j(j);
  end
  f(1)=f(1)-2*T_L/(delta_z_j(1)*delta_z_j(1));
  f(J)=f(J)-2*T_R/(delta_z_j(J)*delta_z_j(J));

  T_j=A\f;

end
