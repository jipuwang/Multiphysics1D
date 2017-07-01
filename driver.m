% Driver
% Case configure options
% fbType options: 'noFeedback','linear','squareRootPlus1'
% assumedSoln options: 'const_cubic','sine_sine','sqrtPlus1_quadratic'
% mocSrc options: 'flat_source','linear_source'
delete diary.txt;
diary('diary.txt')
A={'noFeedback','linear','squareRootPlus1'};
fbTypeS=string(A);
A={'flat_source','linear_source'};
mocSrcS=string(A);
A={'const_cubic','sine_sine','sqrtPlus1_quadratic'};
assumedSolnS=string(A);


% Store the order of accuracy, it could have just been a constant
nCombinations=size(fbTypeS,2)*size(mocSrcS,2)*size(assumedSolnS,2);
order_phi_ensemble=ones(nCombinations,1);
order_T_ensemble=ones(nCombinations,1);
iCombination=0;

for i_fbType=1:size(fbTypeS,2)
  for j_mocSrc=1:size(mocSrcS,2)
    for k_assumedSoln=1:size(assumedSolnS,2)
%       close all;
      fbType=fbTypeS(i_fbType);
      mocSrc=mocSrcS(j_mocSrc);
      assumedSoln=assumedSolnS(k_assumedSoln);
      iCombination=iCombination+1;
      [order_phi_ensemble(iCombination),order_T_ensemble(iCombination)]=...
        converger(fbType,mocSrc,assumedSoln);
    end
  end
end

close all;
diary off;

aa=1;
