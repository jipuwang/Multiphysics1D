% A constant and quadratic solution generator including
    % Discretized analytical solution
    % Manufactured boundary conditions
    % Manufactured source
function [phi0_j_ana,psi_b1_n,psi_b2_n,Q_MMS_j_n,...
          T_j_ana,T_L,T_R,p_MMS_j]=manufacturer_const_quadratic(Tau,mat,J,N)
%   % input parameters
%   if ~exist('Tau','var')
%     Tau=10;
%   end
%   if ~exist('J','var')
%     J=5*2;%*2%*2*2*2*2*2*2*2*2
%   end
%   if ~exist('N','var')
%     N=16;
%   end

  % Material
  Sig_ss_j=mat.Sig_ss_j;
  nuSig_f_j=mat.nuSig_f_j;
  Sig_t_j=mat.Sig_t_j;
  
  % Everything below is determined by the analytical solution
  % In this case, it's the constant angular flux and quadratic temp.
  
  % Discretized analytical solution
  phi0_j_ana=zeros(J,1);
  for j=1:J
    phi0_j_ana(j)=2.0;
  end
  % Boundary Condition
  psi_b1_n=ones(N,1)*1.0; % n=N/2+1:N % mu>0
  psi_b2_n=ones(N,1)*1.0; % n=1:N/2 % mu<0
  % MMS source
  Q_MMS_j_n=ones(J,N)*0.3; % isotropic external source, angular quan.

  %% For TH MMS solution and problem
  % Discretized analytical solution
  T_j_ana=zeros(J,1);
  h_j=ones(J,1)*Tau/J;
  for j=1:J
    T_j_ana(j)=(j*j+j*(j-1)+(j-1)*(j-1))*h_j(j)*h_j(j)/3;
  end
  % Boundary Condition
  T_L=0;
  T_R=100;
  % MMS source
  p_MMS_j=zeros(J,1);
  for j=1:J
    p_MMS_j(j)=2.2;
  end

end

