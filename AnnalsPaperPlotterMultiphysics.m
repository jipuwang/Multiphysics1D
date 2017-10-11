% One case Runner
% Store the order of accuracy, it could have just been a constant
fbTypeS={'noFeedback','linear','squareRootPlus1'};
% fbTypeS=string(A);
mocSrcS={'flat-source','linear-source'};
% mocSrcS=string(A);
assumedSolnS={'const-cubic','sine-sine','sqrtPlus1-quadratic'};
% assumedSolnS=string(A);

%% Plot 1: flat-source, const-cubic, noFeedback
fbType='noFeedback'; 
mocSrc='flat-source'; 
assumedSoln='const-cubic'; 
[order_phi,order_T]=...
  converger(fbType,mocSrc,assumedSoln);

%% Plot 2: flat-source, const-cubic, linear
fbType='linear'; 
mocSrc='flat-source'; 
assumedSoln='const-cubic'; 
[order_phi,order_T]=...
  converger(fbType,mocSrc,assumedSoln);

%% Plot 3: linear-source, sine-sine, noFeedback
fbType='noFeedback'; 
mocSrc='linear-source';
assumedSoln='sine-sine'; 
[order_phi,order_T]=...
  converger(fbType,mocSrc,assumedSoln);

%% Plot 4: linear-source, sine-sine, squareRootPlus1
fbType='squareRootPlus1'; 
mocSrc='linear-source';
assumedSoln='sine-sine'; 
[order_phi,order_T]=...
  converger(fbType,mocSrc,assumedSoln);
%%

