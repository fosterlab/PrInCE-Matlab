function [model,ics] = choosemodel_AIC(cleanchrom,x,aicstring,Ngaussmax)

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
  aicstring = 'AICc';
end

if nargin<4
  Ngaussmax = 5;
end

[n1,n2] = size(cleanchrom);
if n1~=1 && n2~=1 || (n1==1 && n2==1)
  error('Chromatogram must be a 1-dimensional vector.')
end

% transform chromatogram into column vector if necessary
if n2>n1
  cleanchrom = cleanchrom';
end

MaxIter = 3;



%% 1. Fit the five models

% Define fit type
ft{1} = fittype('gauss1');
ft{2} = fittype('gauss2');
ft{3} = fittype('gauss3');
ft{4} = fittype('gauss4');
ft{5} = fittype('gauss5');

% Define fit options
LB = [0.1 0 1]; % H, C, W
UB = [inf length(cleanchrom)+10 inf]; % H, C, W
fo{1} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,1),...
  'Upper',repmat(UB,1,1),'MaxIter',400,'Display','off');
fo{2} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,2),...
  'Upper',repmat(UB,1,1),'MaxIter',400,'Display','off');
fo{3} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,3),...
  'Upper',repmat(UB,1,1),'MaxIter',400,'Display','off');
fo{4} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,4),...
  'Upper',repmat(UB,1,1),'MaxIter',400,'Display','off');
fo{5} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,5),...
  'Upper',repmat(UB,1,1),'MaxIter',400,'Display','off');

I = cleanchrom>0;
msse = nan(5,1);
curveFit = cell(5,1);
R2 = nan(5,1);
for ii=1:Ngaussmax
  iter = 0;
  fitAgain = 1;
  while fitAgain
    iter = iter+1;

    ics = makeics(cleanchrom,ii); % make initial conditions
    set(fo{ii},'startpoint',ics(1:3*ii)); % Include ICs in fit options
    curveFit{ii} = fit(x,cleanchrom,ft{ii},fo{ii}); 
    yhat = feval(curveFit{ii},x);
    msse(ii) = sum((yhat(I) - cleanchrom(I)).^2);
    tmp = corrcoef(yhat,cleanchrom);
    R2(ii) = tmp(1,2).^2;
    
    %fitAgain = R2(ii)<0.5 & iter<=MaxIter;
    cf = coeffvalues(curveFit{ii});
    fitAgain = sum(cf(1:3:end))<0.5 & iter<=MaxIter;
  end
end



%% 2. Pick the best model (AIC / AICc / BIC)
% AIC = (n)log(SSE/n)+2p
% BIC = (n)log(SSE/n)+(p)log(n)

%N = length(cleanchrom);
N = sum(cleanchrom>0);
AIC = ones(5,1) * 10^9;
AICc = ones(5,1) * 10^9;
BIC = ones(5,1) * 10^9;
for ii=1:Ngaussmax
  k = length(coeffvalues(curveFit{ii}));
  AIC(ii) = N*log(msse(ii)/N) + 2*k + 2*k;
  AICc(ii) = N*log(msse(ii)/N) + 2*k + 2*k*(k+1)/(N-k-1) + 2*k;
  BIC(ii) = N*log(msse(ii)/N) + k*log(N) + 2*k;
end

% Throw out any models where the number of non-imputed data points, i.e. the real data, is less than
% the number of parameters.
Ngood = sum(cleanchrom>0.01);
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
  error('Model selection method must be either AIC, AICc, or BIC')
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
    ics((ii-1)*3 + 1) = pks(ii) + rand-.5; % A, height
    ics((ii-1)*3 + 2) = pksI(ii) + rand*3-1.5; % mu, center
    ics((ii-1)*3 + 3) = ss + rand-.5; % sigma, width
  else
    ics((ii-1)*3 + 1) = max(cleanchrom) + rand-.5; % A, height
    ics((ii-1)*3 + 2) = rand*length(cleanchrom) + rand*3-1.5; % mu, center
    ics((ii-1)*3 + 3) = ss + rand-.5; % sigma, width
  end
end

ics(Ngauss*3+1 : end) = 0;

