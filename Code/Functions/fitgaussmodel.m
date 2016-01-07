function [Coef,SSE,adjrsquare,fit_flag] = fitgaussmodel(cleanchrom,model)

% Selects the number of Gaussians that best fits a single chromatogram.
%
% Model selection is done by generating a distribution of SSEs for each model type. The distribution
% of SSEs are generated by fitting each model many times on a subset of the data. This is a hold out
% analysis. These distributions are finally run through Nick's model selection algorithm to pick the
% model type. Each fit requires initial conditions (ICs, also called StartingPoints). These are
% generated by function gaussfitICs.m.
%
% Input:
%   cleanchrom: nx1 cleaned chromatogram, where n is the number of fractions.
%   maxIter: integer, the number of iterations to run the holdout analysis.
% Output:
%   model: integer from 1 to 5, which is the number of Gaussians for the best model.
%   holdoutSSE: 1x5 cell, where each element is a distribution SSE of length MaxIter.
%
% Adapted from Nichollas Scott's Gaus_build_24_1.m.
% Made by Greg Stacey on Nov 25 2015.



%% 0. Initialize

if nargin<1
  error('Must input a chromatogram.')
end

if nargin<2
  disp('    Selecting best model via holdout analysis with 500 iterations...')
  model = choosemodel_holdout(clean_chromatogram,500);
end

[n1,n2] = size(cleanchrom);
if n1~=1 && n2~=1 || (n1==1 && n2==1)
  error('Chromatogram must be a 1-dimensional vector.')
end

% transform chromatogram into column vector if necessary
if n2>n1
  cleanchrom = cleanchrom';
end



%% Fit the best model.

x = (1:length(cleanchrom))';

if model.Ngauss>0
  % A reasonable model was selected, i.e. 1-5 Gaussians
  
  try
    fit_flag = 1;
    [curveFit,goodness] = fit(x,cleanchrom,model.ft,model.fo);
  catch
    fit_flag = 0;
    %Use dummy values
    dummy_x= (1:1:65)';
    dummy_cleanchrom= ones(1,65)';
    [curveFit,goodness] = fit(dummy_x,dummy_cleanchrom,model.ft,model.fo);
  end  
  Coef=coeffvalues(curveFit);
  SSE=goodness.sse;
  adjrsquare=goodness.adjrsquare;
else
  % Failed to find a reasonable model.
  
  Coef = 9999;
  SSE = 9999;
  adjrsquare = 9999;
end

