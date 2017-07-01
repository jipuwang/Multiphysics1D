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
order_phi_ensemble=zeros(nCombinations,1);
order_T_ensemble=zeros(nCombinations,1);
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

order_phi_gold=[100;2;2;100;4;4;2;2;2;2;2;2;2;2;2;2;2;2;];
order_T_gold=[2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;];

regressionTestPasses=true;
for iCombination=1:nCombinations
  if order_phi_gold(iCombination)~=100 %100 is considered an automatic pass
    if abs(order_phi_ensemble(iCombination)-order_phi_gold(iCombination))>0.1
      regressionTestPasses=false;
    end
    if abs(order_T_ensemble(iCombination)-order_T_gold(iCombination))>0.1
      regressionTestPasses=false;
    end
  end
end

if regressionTestPasses==true
  display 'All Tests Passed'
else
  display 'Test Failed'
end
    