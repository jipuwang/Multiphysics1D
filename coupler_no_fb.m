%% Info
% This is a no_feedback_coupler, it does a lot of stuff for now including:
% Pass the source to neutronics module and get the solution phi_j.
% Build the fission heat source.
% Pass the fission heat source and manufactured heat source (or the
% addition of the two) to the heat conduction module and get the solution
% T_j.

% The name of the coupler should reveal the constant_quadratic case.

function [phi0_j,T_j]=coupler_no_fb(Tau,mat,J,N,psi_b1_n,psi_b2_n,Q_MMS_j_n,...
          T_L,T_R,p_MMS_j)
%   % input parameters
%   if ~exist('J','var')
%     Tau=10;
%   end
%   if ~exist('J','var')
%     J=5*2;%*2%*2*2*2*2*2*2*2*2
%   end
%   if ~exist('N','var')
%     N=16;
%   end
%   if ~exist('psi_b1_n','var')
%     psi_b1_n=ones(N,1)*1.0; % the first/negative half is not useful; n=N/2+1:N % mu>0
%   end
%   if ~exist('psi_b2_n','var')
%     psi_b2_n=ones(N,1)*1.0; % the second/positive half is not useful; n=1:N/2 % mu<0
%   end
%   if ~exist('Q_MMS_j_n','var')
%     Q_MMS_j_n=ones(J,N)*0.3; % isotropic external source, angular quan.
%   end
%   if ~exist('T_L','var')
%     T_L=0;
%   end
%   if ~exist('T_R','var')
%     T_R=100;
%   end
%   if ~exist('p_MMS_j','var')
%     % define p_MMS_j
%     p_MMS_j=zeros(J,1);
%     for j=1:J
%       p_MMS_j(j)=2.2;
%     end
%   end

%% Call the MoC module to get the flux
  phi0_j=MoC_module(Tau,mat,J,N,...
           psi_b1_n,psi_b2_n,Q_MMS_j_n);

%% The coupling
  pTriplePrime_j=1*0.1*phi0_j; % kappa=1, sigma_f=0.1;

%% Call the heat conduction module to get the temperature
  T_j=heat_cond_module(Tau,J,T_L,T_R,pTriplePrime_j,p_MMS_j);

end
