% 1D MoC Module
% Input: 
%   Geometry Tau
%   Spatial discretization J (or mesh size)
%   Angular discretization N
%   Material: all cross sections and stuff
%   Boundary condition
%   Distributed source, can be MMS
% Output: 
%   Cell-averaged scalar flux

function [phi0_j]=MoC_LS_module(J,N,Tau,mat,...
       psi_b1_n,psi_b2_n,Q_MMS_j_n,Q_MMS_hat_j_n)

%   Input parameter
  if ~exist('Tau','var')
  Tau=10;
  end
  if ~exist('J','var')
  J=5*2;%*2%*2*2*2*2*2*2*2*2
  end
  if ~exist('N','var')
  N=16;
  end
  if ~exist('mat','var')
  % Material
  field1='Sig_t_j';      value1=ones(J,1);
  field2='Sig_ss_j';     value2=ones(J,1)*0.5;
  field3='Sig_gamma_j';    value3=ones(J,1)*0.4;
  field4='Sig_f_j';      value4=ones(J,1)*0.1;
  field5='nuSig_f_j';    value5=ones(J,1)*0.2;
  field6='thermal_cond_k_j'; value6=ones(J,1);
  field7='kappaSig_f_j';   value7=ones(J,1)*0.1; % kappa=1.0;
  mat = struct(field1,value1,field2,value2,field3,value3,... 
    field4,value4,field5,value5,field6,value6,field7,value7);
  end
  if ~exist('psi_b1_n','var')
    psi_b1_n=ones(N,1)*1.0;
  end
  if ~exist('psi_b2_n','var')
    psi_b2_n=ones(N,1)*1.0;
  end
  if ~exist('Q_MMS_j_n','var')
    Q_MMS_j_n=ones(J,N)*0.3; % removed *2.0 (angular quantity)
  end
  if ~exist('Q_MMS_hat_j_n','var')
    Q_MMS_hat_j_n=ones(J,N)*0.1; % removed *2.0 (angular quantity)
  end
  
  % Material
  Sig_ss_j=mat.Sig_ss_j;
  nuSig_f_j=mat.nuSig_f_j;
  Sig_t_j=mat.Sig_t_j;
%   Sig_ss_j=ones(J,1)*0.5;
%   nuSig_f_j=ones(J,1)*0.2;
%   Sig_t_j=ones(J,1);

  Sig_t_inv_j=1./Sig_t_j;
  
  % Default variables, can be customized. 
  maxIterate=2000;
  epsilon_phi=1e-10;
  delta=1E-13;
  [mu_n,weight_n]=lgwt(N,-1,1); mu_n=flipud(mu_n);
  
  h_j=ones(J,1)*Tau/J;
  % N rays to trace, each angle has only 1 ray, no ray-spacing
  % n for each angle, and j for FSR region index
  segLen_j_n=zeros(J,1);
  for n=1:N
    for j=1:J
      segLen_j_n(j,n)=h_j(j)/abs(mu_n(n));
    end
  end

% From copy and paste
phi0_old_j=ones(1,J)*1.0; % so the 1st dimension is consistently the angle. 
phi0_old_hat_j=ones(1,J)*1.0; % so the 1st dimension is consistently the angle. 

Q_x_j_n=zeros(J,N); % these are actually angular quantities, already have 0.5's in them.
q_j_n=zeros(J,N);
q_sm_j_n=zeros(J,N);
% new quantities for linear source
Q_x_hat_j_n=zeros(J,N);
q_hat_j_n=zeros(J,N);
q_sm_hat_j_n=zeros(J,N);
%% temporary
Q_MMS_hat_j_n=zeros(J,N);
%%
for iIterate=1:maxIterate
  for j=1:J
    for n=1:N
      Q_x_j_n(j,n)=0.5*(Sig_ss_j(j)+nuSig_f_j(j))*phi0_old_j(j)+Q_MMS_j_n(j,n);
      q_j_n(j,n)=Q_x_j_n(j,n);
      q_sm_j_n(j,n)=q_j_n(j,n);

      Q_x_hat_j_n(j,n)=0.5*(Sig_ss_j(j)+nuSig_f_j(j))*phi0_old_hat_j(j)+Q_MMS_hat_j_n(j,n); 
      q_hat_j_n(j,n)=Q_x_hat_j_n(j,n)/(h_j(j)*h_j(j)/12);
      q_sm_hat_j_n(j,n)=q_hat_j_n(j,n)*(mu_n(n));   % NO ABS IS NEEDED!
    end
  end
%   phi_j_old_hat
%   q_n_j
%   q_n_j_hat

  phi0_new_j=zeros(1,J);
  phi0_hat_new_j=zeros(1,J);
  % ray tracing
  for n=N/2+1:N
    psi_in=psi_b1_n(n);
    for j=1:J
      tau_temp=Sig_t_j(j)*segLen_j_n(j,n);
      F1=1-exp(-tau_temp);
      F2=2*(tau_temp-F1)-tau_temp*F1;
      psi_out=psi_in+(q_sm_j_n(j,n)*Sig_t_inv_j(j)-psi_in)*F1...
        +(q_sm_hat_j_n(j,n)*0.5*Sig_t_inv_j(j)*Sig_t_inv_j(j))*F2;
      psi_avg=q_sm_j_n(j,n)*Sig_t_inv_j(j)+(psi_in-psi_out)/tau_temp;
      phi0_new_j(j)=phi0_new_j(j)+weight_n(n)*psi_avg;
      G1=1+tau_temp*0.5-(1+1/tau_temp)*F1;
      G2=2/3*tau_temp-(1+2/tau_temp)*G1;
      psi_hat=psi_in*0.5*segLen_j_n(j,n) ...
        + (q_sm_j_n(j,n)*Sig_t_inv_j(j)-psi_in)*G1*Sig_t_inv_j(j) ...
        + (q_sm_hat_j_n(j,n)*0.5*Sig_t_inv_j(j)*Sig_t_inv_j(j))*segLen_j_n(j,n)*G2; %changed here

      phi0_hat_new_j(j)=phi0_hat_new_j(j)+...
        weight_n(n)*(-h_j(j)*0.5*psi_avg+abs(mu_n(n))*psi_hat);
      psi_in=psi_out;
    end
  end
  
  for n=1:N/2
    psi_in=psi_b2_n(n);
    for j=J:-1:1
      tau_temp=Sig_t_j(j)*segLen_j_n(j,n);
      F1=1-exp(-tau_temp);
      F2=2*(tau_temp-F1)-tau_temp*F1;
      psi_out=psi_in+(q_sm_j_n(j,n)*Sig_t_inv_j(j)-psi_in)*F1 ...
        +(q_sm_hat_j_n(j,n)*0.5*Sig_t_inv_j(j)*Sig_t_inv_j(j))*F2;
      psi_avg=q_sm_j_n(j,n)*Sig_t_inv_j(j)+(psi_in-psi_out)/tau_temp;
      phi0_new_j(j)=phi0_new_j(j)+weight_n(n)*psi_avg;
      G1=1+tau_temp*0.5-(1+1/tau_temp)*F1;
      G2=2/3*tau_temp-(1+2/tau_temp)*G1;
      psi_hat=psi_in*0.5*segLen_j_n(j,n) ...
        + (q_sm_j_n(j,n)*Sig_t_inv_j(j)-psi_in)*G1*Sig_t_inv_j(j) ...
        + (q_sm_hat_j_n(j,n)*0.5*Sig_t_inv_j(j)*Sig_t_inv_j(j))*segLen_j_n(j,n)*G2; %changed here

      phi0_hat_new_j(j)=phi0_hat_new_j(j)+...
        weight_n(n)*(+h_j(j)*0.5*psi_avg-abs(mu_n(n))*psi_hat);
      psi_in=psi_out;
    end
  end
  
  % test for convergence
  error=norm(phi0_new_j-phi0_old_j)
  if error<epsilon_phi
    break;
  end
  phi0_old_j=phi0_new_j;
  phi0_old_hat_j=phi0_hat_new_j;
end
error
% phi_j_old=phi_j_new;
phi0_new_j=phi0_new_j';
display(iIterate);
% figure(19);
% plot(phi0_new_j,'*-')
% openvar('phi_j_new')

  phi0_j=phi0_new_j;
  
end
