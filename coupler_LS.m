%% Info
% This is a coupler_with_fb, it does the following things:
% 1. Solve for phi0_j with known source, e.g., MMS source.
% 2. Build the fission heat source for TH problem.
% 3. Solve for T_j with the above heat soruce (+ MMS soruce optionally).
% 4. Use T_j to update the xs
% 5. Return to S1 until convergence in both phi0_j and T_j.

function [phi0_j,T_j]=...
  coupler_LS(J,N,Tau,mat,psi_b1_n,psi_b2_n,Q_MMS_j_n,Q_MMS_hat_j_n, ...
          T_L,T_R,p_MMS_j,fbType,phi0_j_ana,T_j_ana)
  % input parameters
  if ~exist('J','var')
    J=10;
  end
  if ~exist('N','var')
    N=16;
  end
  if ~exist('Tau','var')
    Tau=10;
  end
  if ~exist('mat','var')
    % Material
    field1='Sig_t_j';          value1=ones(J,1);
    field2='Sig_ss_j';         value2=ones(J,1)*0.5;
    field3='Sig_gamma_j';      value3=ones(J,1)*0.4;
    field4='Sig_f_j';          value4=ones(J,1)*0.1;
    field5='nuSig_f_j';        value5=ones(J,1)*0.2;
    field6='thermal_cond_k_j'; value6=ones(J,1);
    field7='kappaSig_f_j';     value7=ones(J,1)*0.1; % kappa=1.0;
    mat = struct(field1,value1,field2,value2,field3,value3,... 
      field4,value4,field5,value5,field6,value6,field7,value7);
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
  if ~exist('fbType','var')
    fbType='linear';
  end

  % Start the Picard Iteration
  T_old_j=zeros(J,1);
  phi0_old_j=zeros(J,1);
  isConverged=false;
  % save the capture xs before the correction.
  Sig_gamma_ref_j=mat.Sig_gamma_j;
  while ~isConverged
    %% Call the MoC module to get the flux
    phi0_j=MoC_LS_module(J,N,Tau,mat,...
             psi_b1_n,psi_b2_n,Q_MMS_j_n,Q_MMS_hat_j_n);
%     phi0_j=phi0_j_ana;
    error_phi=norm(phi0_j-phi0_old_j)/sqrt(J);
    phi0_old_j=phi0_j;

    %% The coupling
    % Build the heat source for TH
    pTriplePrime_j=mat.kappaSig_f_j.*phi0_j; % kappa=1, sigma_f=0.1;

    %% Call the heat conduction module to get the temperature
    T_j=heat_cond_module(J,Tau,mat,T_L,T_R,pTriplePrime_j,p_MMS_j);
%     T_j=T_j_ana;
    error_T=norm(T_j-T_old_j)/sqrt(J);
    T_old_j=T_j;

    %% check convergence
    if (error_phi<1e-12 && error_T<1e-12)
      break;
    end
    
    %% update cross section
    switch char(fbType)
      case 'noFeedback'
        break;
      case 'linear'
        % gamma_coeff*(T_j_new-T0) is fb.
        T0=50;
        gamma_coeff=0.004;
%         gamma_coeff=0.000;
        mat.Sig_gamma_j=Sig_gamma_ref_j+gamma_coeff*(T_j-T0); 
      case 'squareRootPlus1'
        T0=50;
        mat.Sig_gamma_j=Sig_gamma_ref_j.*sqrt((T0+1)./(T_j+1));
    end
    mat.Sig_t_j=mat.Sig_ss_j+mat.Sig_gamma_j+mat.Sig_f_j;

  end

end
