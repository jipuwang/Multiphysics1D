% A constant and quadratic solution generator including
    % Discretized analytical solution
    % Manufactured boundary conditions
    % Manufactured source
function [phi0_MMS_j,psi_b1_n,psi_b2_n,Q_MMS_j_n,...
          T_MMS_j,T_L,T_R,q_MMS_j]=manufacturer_sine_sine(J,N,Tau,mat)
  % input parameters
  if ~exist('J','var')
    J=5*2;%*2%*2*2*2*2*2*2*2*2
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
    field6 = 'kappaSig_f_j'; value6 = ones(J,1)*0.1; % kappa=1.0;
    mat = struct(field1,value1,field2,value2,field3,value3,... 
      field4,value4,field5,value5,field6,value6);
  end

  % Material
  Sig_ss_j=mat.Sig_ss_j;
  nuSig_f_j=mat.nuSig_f_j;
  kappaSig_f_j=mat.kappaSig_f_j;
  Sig_f_j=mat.Sig_f_j;
  Sig_t_j=mat.Sig_t_j;
  k_F=mat.thermal_cond_k_j(1);
  
  h=Tau/J;
  [mu_n,weight_n]=lgwt(N,-1,1); mu_n=flipud(mu_n);
  
  %% For MoC MMS solution and problem  
  % Everything below is determined by the analytical solution
  % In this case, it's the constant angular flux and quadratic temp.
  
  % Assumed manufactured solution \psi(x,\mu)=sin(\pi*x/Tau), 0<x<Tau
  % Define function handle for angular flux, phi assumed to be 2*psi
  psi_MMS =@(x) sin(pi*x/Tau);
  psi_MMS_Diff =@(x) pi/Tau*cos(pi*x/Tau);
  
  % Boundary condition and source
  % psi expression evaluated at x=0
  psi_b1_n=psi_MMS(0)*ones(N,1); % n=N/2+1:N % mu>0
  % psi expression evaluated at x=Tau
  psi_b2_n=psi_MMS(Tau)*ones(N,1); % n=1:N/2 % mu<0
  
  Q_MMS_j_n=zeros(J,N); % preallocate memory, avg'ed over tau_(j-1/2) and tau_(j+1/2)
  % MMS source: mu_n * derivative(psi_MMS) ...
  % + Sig_t* psi_MMS - (Sig_ss+nuSig_f)*0.5*phi0_MMS
  
  phi0_MMS_j=zeros(J,1);
  psi_MMS_j=zeros(J,1);
  psi_MMS_Diff_j=zeros(J,1); % This is needed to build MMS source

  for j=1:J
    x_L=(j-1)*h;x_R=j*h;
    psi_MMS_j(j)=1/h*integral(psi_MMS,x_L,x_R);
    phi0_MMS_j(j)=2.0*psi_MMS_j(j);
    psi_MMS_Diff_j(j)= 1/h*integral(psi_MMS_Diff,x_L,x_R);
    for n=1:N
    Q_MMS_j_n(j,n)=mu_n(n)* psi_MMS_Diff_j(j) +Sig_t_j(j)*psi_MMS_j(j) ...
      -(Sig_ss_j(j)+nuSig_f_j(j))*0.5* phi0_MMS_j(j);
    end % n
  end % j

  %% For TH MMS solution and problem
  % Assumed manfuactured solution T(x)=sin(pi*x/Tau), 0<x<Tau
  T_MMS =@(x) sin(pi*x/Tau);
  T_MMS_xx =@(x) -(pi*pi/Tau/Tau)*sin(pi*x/Tau);
  
  % Boundary condition and source
  % Left boundary, T_MMS evaluated at x=0;
  T_L=T_MMS(0);
  % Right boundary, T_MMS evaluated at x=Tau;
  T_R=T_MMS(Tau);
  
  % Discretized MMS solution
  T_MMS_j=zeros(J,1);
  T_MMS_xx_j=zeros(J,1);

  % MMS source
  q_MMS_j=zeros(J,1);
  for j=1:J
    x_L=(j-1)*h;x_R=j*h;
    T_MMS_j(j)=1/h*integral(T_MMS,x_L,x_R);
    T_MMS_xx_j(j)=1/h*integral(T_MMS_xx,x_L,x_R);
    q_MMS_j(j)=k_F*T_MMS_xx_j(j)+kappaSig_f_j(j)*phi0_MMS_j(j);
  end  

end

