%% Info
% This is a no_feedback_coupler, it does a lot of stuff for now including:
% Pass the source to neutronics module and get the solution phi_j.
% Build the fission heat source.
% Pass the fission heat source and manufactured heat source (or the
% addition of the two) to the heat conduction module and get the solution
% T_j.

% The name of the coupler should reveal the constant_quadratic case.

function [phi0_j,T_j]=coupler_no_fb(J,N,Tau,mat,psi_b1_n,psi_b2_n,Q_MMS_j_n,...
          T_L,T_R,p_MMS_j)
  % input parameters
  if ~exist('J','var')
    Tau=10;
  end
  if ~exist('N','var')
    N=16;
  end
  if ~exist('Tau','var')
    Tau=10;
  end
  if ~exist('mat','var')
    % Material
    field1 = 'Sig_ss_j';  value1 = ones(J,1)*0.5;
    field2 = 'nuSig_f_j';  value2 = ones(J,1)*0.2;
    field3 = 'Sig_t_j';  value3 = ones(J,1);
    field4 = 'thermal_cond_k_j'; value4 = ones(J,1);
    field5 = 'Sig_f_j'; value5 = ones(J,1)*0.1;
    mat = struct(field1,value1,field2,value2,field3,value3,... 
      field4,value4,field5,value5);
  end
  if ~exist('psi_b1_n','var')
    psi_b1_n=ones(N,1)*1.0; % the first/negative half is not useful; n=N/2+1:N % mu>0
  end
  if ~exist('psi_b2_n','var')
    psi_b2_n=ones(N,1)*1.0; % the second/positive half is not useful; n=1:N/2 % mu<0
  end
  if ~exist('Q_MMS_j_n','var')
    Q_MMS_j_n=ones(J,N)*0.3; % isotropic external source, angular quan.
  end
  if ~exist('T_L','var')
    T_L=0;
  end
  if ~exist('T_R','var')
    T_R=100;
  end
  if ~exist('p_MMS_j','var')
    % define p_MMS_j
    p_MMS_j=zeros(J,1);
    for j=1:J
      p_MMS_j(j)=2.2;
    end
  end

%% Call the MoC module to get the flux
  phi0_j=MoC_module(J,N,Tau,mat,...
           psi_b1_n,psi_b2_n,Q_MMS_j_n);

%% The coupling
  pTriplePrime_j=1*0.1*phi0_j; % kappa=1, sigma_f=0.1;

%% Call the heat conduction module to get the temperature
  T_j=heat_cond_module(J,Tau,mat,T_L,T_R,pTriplePrime_j,p_MMS_j);

end
