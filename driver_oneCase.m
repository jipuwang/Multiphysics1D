% One case Runner
% Store the order of accuracy, it could have just been a constant
fbTypeS={'noFeedback','linear','squareRootPlus1'};
% fbTypeS=string(A);
mocSrcS={'flat-source','linear-source'};
% mocSrcS=string(A);
assumedSolnS={'const-cubic','sine-sine','sqrtPlus1-quadratic'};
% assumedSolnS=string(A);


fbType=fbTypeS(2); % linear
mocSrc=mocSrcS(2); % linear-source
assumedSoln=assumedSolnS(1); % const-cubic
[order_phi,order_T]=...
  converger(fbType,mocSrc,assumedSoln);

a=1;