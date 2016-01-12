function [model,ics] = choosemodel_AIC(cleanchrom,x,aicstring)

% Selects the number of Gaussians that best fits a single chromatogram.
%
% This code fits 5 models (1 Gaussian, 2 Gaussians, ... , 5 Gaussians) and picks the most
% appropriate model via AIC. Each model is fit using different initial parameter guesses until a
% good fit is achieved.
%
% Input:
%   cleanchrom: nx1 cleaned chromatogram, where n is the number of fractions.
%   maxIter: integer, the number of iterations to run until a good fit is achieved
% Output:
%   model: structure with three fields.
%       model.Ngauss - integer from 1 to 5, which is the number of Gaussians for the best model
%       model.fo - fit options for the best model, including ICs
%       model.ft - fit type of the best model
%
%
% Adapted from Nichollas Scott's Gaus_build_24_1.m.
% Made by Greg Stacey on Nov 25 2015.
%
% AIC reference, http://statweb.stanford.edu/~jtaylo/courses/stats203/notes/selection.pdf



%% 0. Initialize

if nargin<1
  error('Must input a chromatogram.')
end

if nargin<2
  x = 1:length(cleanchrom);
end
if isempty(x)
  x = 1:length(cleanchrom);
end

if nargin<3
  aicstring = 'AIC';
end

[n1,n2] = size(cleanchrom);
if n1~=1 && n2~=1 || (n1==1 && n2==1)
  error('Chromatogram must be a 1-dimensional vector.')
end

% transform chromatogram into column vector if necessary
if n2>n1
  cleanchrom = cleanchrom';
end

MaxIter = 2;



%% 1. Fit the five models

% Define fit type
ft{1} = fittype('gauss1');
ft{2} = fittype('gauss2');
ft{3} = fittype('gauss3');
ft{4} = fittype('gauss4');
ft{5} = fittype('gauss5');

% Define fit options
fo{1} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',zeros(1,3),'MaxIter',400);
fo{2} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',zeros(1,6),'MaxIter',400);
fo{3} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',zeros(1,9),'MaxIter',400);
fo{4} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',zeros(1,12),'MaxIter',400);
fo{5} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',zeros(1,15),'MaxIter',400);

I = cleanchrom>0;
msse = nan(5,1);
curveFit = cell(5,1);
R2 = nan(5,1);
for ii=1:5
  iter = 0;
  fitAgain = 1;
  while fitAgain
    ics = makeics(cleanchrom,ii); % make initial conditions
    set(fo{ii},'startpoint',ics(1:3*ii)); % Include ICs in fit options
    curveFit{ii} = fit(x,cleanchrom,ft{ii},fo{ii}); 
    yhat = feval(curveFit{ii},x);
    msse(ii) = sum((yhat(I) - cleanchrom(I)).^2);
    tmp = corrcoef(yhat,cleanchrom);
    R2(ii) = tmp(1,2).^2;
    
    iter = iter+1;
    fitAgain = R2(ii)<0.5 & iter<=MaxIter;
  end
end



%% 2. Pick the best model (AIC / AICc / BIC)
% AIC = (n)log(SSE/n)+2p
% BIC = (n)log(SSE/n)+(p)log(n)

%N = length(cleanchrom);
N = sum(cleanchrom>0);
AIC = zeros(5,1);
AICc = zeros(5,1);
BIC = zeros(5,1);
for ii=1:5
  k = length(coeffvalues(curveFit{ii}));
  AIC(ii) = N*log(msse(ii)/N) + 2*k + 5*k;
  AICc(ii) = N*log(msse(ii)/N) + 2*k + 2*k*(k+1)/(N-k-1) + 5*k;
  BIC(ii) = N*log(msse(ii)/N) + k*log(N) + 5*k;
end

% Throw out any models where the number of non-imputed data points, i.e. the real data, is less than
% the number of parameters.
Ngood = sum(cleanchrom>0);
maxModelSize = floor(Ngood/3);
AIC(maxModelSize:end) = inf;
AICc(maxModelSize:end) = inf;
BIC(maxModelSize:end) = inf;

if strcmp(aicstring,'AIC')
  [~,Ngauss] = min(AIC);
elseif strcmp(aicstring,'AICc')
  [~,Ngauss] = min(AICc);
elseif strcmp(aicstring,'BIC')
  [~,Ngauss] = min(BIC);
else
  error('Input aicstring must be either AIC, AICc, or BIC')
end

model.coeffs = coeffvalues(curveFit{Ngauss});
model.SSE = msse(Ngauss);
model.adjrsquare = R2(Ngauss);
model.curveFit = curveFit{Ngauss};
model.AIC = AIC;
model.AICc = AICc;
model.BIC = BIC;

model.Ngauss = Ngauss;
if model.Ngauss>0
  model.fo = fo{Ngauss};
  model.ft = ft{Ngauss};
else
  model.fo = [];
  model.ft = [];
end



function ics = makeics(cleanchrom,Ngauss)

% Make initial conditions, ics
[pks, pksI] = findpeaks(cleanchrom);
d = diff(pksI);
[~,I] = sort(d,'descend');
pks = pks([1; I+1]);
pksI = pksI([1; I+1]);

ss = 2;
ics = nan(1,15);
for ii=1:5
  if ii<=length(pks)
    ics((ii-1)*3 + 1) = pks(ii); % A, height
    ics((ii-1)*3 + 2) = pksI(ii); % mu, center
    ics((ii-1)*3 + 3) = ss; % sigma, width
  else
    ics((ii-1)*3 + 1) = max(cleanchrom); % A, height
    ics((ii-1)*3 + 2) = rand*length(cleanchrom); % mu, center
    ics((ii-1)*3 + 3) = ss; % sigma, width
  end
end

ics(Ngauss*3+1 : end) = 0;

