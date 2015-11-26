function holdoutSSE = holdoutSSE(cleanchrom,MaxIter,ICs)

% Using a single chromatogram, uses a holdout analysis to generate a distribution of SSE for a 
% number of models.
%
% Input: 
%   cleanchrom: nx1 cleaned chromatogram, where n is the number of fractions.
%   maxIter: scalar MaxIter, which is the number of iterations to run the holdout analysis.
%   ICs: 1x5 cell, where each element is a vector containing the ICs for a single Gaussian model.
% Output: 
%   holdoutSSE: 1x5 cell, where each element is a vector containing the ICs.
%
% Adapted from Nichollas Scott's Gaus_build_24_1.m.
% Made by Greg Stacey on Nov 25 2015.



%% 0. Initialize

if nargin<1
  error('Must input a chromatogram.')
end

if nargin<2
  MaxIter = 500;
  disp('Setting number of iterations to 500...')
end

[n1,n2] = size(cleanchrom);
if n1~=1 && n2~=1 || (n1==1 && n2==1)
  error('Chromatogram must be a 1-dimensional vector.')
end

% transform chromatogram into row vector if necessary
if n1>n2
  cleanchrom = cleanchrom';
end


%% Run holdout analysis

% Define fit type
ft{1} = fittype('gauss1');
ft{2} = fittype('gauss2');
ft{3} = fittype('gauss3');
ft{4} = fittype('gauss4');
ft{5} = fittype('gauss5');

% Define fit options
fo{1} = fitoptions('method','NonlinearLeastSquares','Robust','off','Lower',[0.1 0.1 0.1]);
fo{2} = fitoptions('method','NonlinearLeastSquares','Robust','off','Lower',[0.1 0.1 0.1 0.1 0.1 0.1]);
fo{3} = fitoptions('method','NonlinearLeastSquares','Robust','off','Lower',[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]);
fo{4} = fitoptions('method','NonlinearLeastSquares','Robust','off','Lower',[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]);
fo{5} = fitoptions('method','NonlinearLeastSquares','Robust','off','Lower',[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]);

x = 1:length(cleanchrom);

% Fitting using crossvalidation, done via leave one out analysis
% where 10% of the curve is removed to check
% THIS BUILDS Output_SSE
holdoutSSE = cell(1,5);
for ii=1:5
  jj=1;
  while jj <= MaxIter
    try %attempt the hold out analysis
      jj=jj+1;
      train = crossvalind('HoldOut',N2, 0.1);
      set(fo{ii},'startpoint',ICs{ii});
      curveFit = fit(x(train),cleanchrom(train),ft{ii},fo{ii});
      yhat = feval(curveFit,x2);
      holdoutSSE{ii}(jj) =  sum((yhat - y2).^2);
    catch
%       %if the crossvalidation process leads to an error repeat till       % I dont think this is
%       %the required successful initrations have been done                 % needed, since jj won't
%       if jj>1 % this ensures  i = 1 will not be turned to 0               % be incremented if the
%         jj=jj-1;                                                          % try block is skipped.
%       end                                                                 %
    end
  end
end

