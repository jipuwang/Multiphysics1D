% Driver
% Case configure options
% fbType options: 'noFeedback','linear','squareRootPlus1'
% assumedSoln options: 'const_quadratic','sine_sine','sqrtPlus1_quadratic'
% mocSrc options: 'flat_source','linear_source'
delete diary.txt;
diary('diary.txt')
A={'noFeedback','linear','squareRootPlus1'};
fbTypeS=string(A);
A={'flat_source','linear_source'};
mocSrcS=string(A);
A={'const_quadratic','sine_sine','sqrtPlus1_quadratic'};
assumedSolnS=string(A);

for i_fbType=1:size(fbTypeS,2)
  for j_mocSrc=1:size(mocSrcS,2)
    for k_assumedSoln=1:size(assumedSolnS,2)
      clearvars -except i_fbType j_mocSrc k_assumedSoln ...
        fbTypeS mocSrcS assumedSolnS; 
      close all;
      fbType=fbTypeS(i_fbType);
      mocSrc=mocSrcS(j_mocSrc);
      assumedSoln=assumedSolnS(k_assumedSoln);
      converger;
    end
  end
end

close all;
diary off;
