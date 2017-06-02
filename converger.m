% Grid refiner is for grid refinement analysis. 
% Right now, it needs to be aware of the geometry and the material. So it
% can call the manufacturer to get the MMS problem and solution
% The geometry and material also need to be passed to the coupler, so the
% coupler can keep passing the info on to the modules, because it's part of
% the problem description. 
% It needs to know the geometry and is responsible for generating the grid
% and pass the grid information to the coupler. 

nGrids=10;%8;

% Geometry
Tau=10; 
% Material
field1 = 'Sig_ss_j';  value1 = ones(J,1)*0.5;
field2 = 'nuSig_f_j';  value2 = ones(J,1)*0.2;
field3 = 'Sig_t_j';  value3 = ones(J,1);
mat = struct(field1,value1,field2,value2,field3,value3);

error_phi_n=ones(nGrids,1);
error_T_n=ones(nGrids,1);
gridMeshSize_n=ones(nGrids,1);
N=16; % angular discretization, fixed not refined. 
for iGrid=1:nGrids
  J=5*2^iGrid;
  gridMeshSize_n(iGrid)=Tau/J;
  
  % call the manufacturer to get MMS problem and solution
  [phi0_j_ana,psi_b1_n,psi_b2_n,Q_MMS_j_n,...
          T_j_ana,T_L,T_R,p_MMS_j]=manufacturer_const_quadratic(Tau,mat,J,N);
  
  % call the coupler to solve the above manufactured problem
  [phi0_j,T_j]=coupler_no_fb(Tau,mat,J,N,psi_b1_n,psi_b2_n,Q_MMS_j_n,...
          T_L,T_R,p_MMS_j);

  % Calculate the error compared to manufactured solution
  error_phi_n(iGrid)=norm(phi0_j-phi0_j_ana,2)/sqrt(J);
  error_T_n(iGrid)=norm(T_j-T_j_ana,2)/sqrt(J);
end

% Calculate the order of accuracy
order_phi=ones(nGrids-1,1);
order_T=ones(nGrids-1,1);
for j=1:nGrids-1
  order_phi(j)=log(error_phi_n(j)/error_phi_n(j+1))/log(gridMeshSize_n(j)/gridMeshSize_n(j+1));
  order_T(j)=log(error_T_n(j)/error_T_n(j+1))/log(gridMeshSize_n(j)/gridMeshSize_n(j+1));
end

% Display the result
error_phi_n
error_T_n
order_phi
order_T
  
% Visualize the results
figure(1)
loglog(gridMeshSize_n,error_phi_n,'-*');
% title('scalar flux error convergence');
xlabel('mesh size in mean free path');
ylabel('2-norm error of scalar flux');
figure(2)
loglog(gridMeshSize_n,error_T_n,'-*');
% title('temperature error convergence');
xlabel('mesh size in mean free path');
ylabel('2-norm error of temperature');

