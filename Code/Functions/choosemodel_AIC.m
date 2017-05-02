function [model,ics] = choosemodel_AIC(coFracProfile,x,aicstring,Ngaussmax)
%CHOOSEMODEL_AIC Fit and choose Gaussian mixture model.
%   model=CHOOSEMODEL_AIC(coFracProfile,x,aicstring,Ngaussmax) fits 5 or 
%   fewer Gaussian mixture models to a single co-fractionation profile. 
%   These models are a mixture of 1, 2, 3, 4, or 5 Gaussians, respectively. 
%   The best model, chosen via AIC, AICc, or BIC, adequately describes the 
%   profile with a minimum of Gaussians. model is a structure describing the 
%   fitted models with fields:
%     coeffs      coefficients of chosen model
%     CIS         95% confidence intervals of coefficients
%     SSE         Sum of squared errors of chosen model
%     adjrsquare  Adjusted R^2 of chosen model
%     Ngauss      Number of Gaussians in chosen model
%     AIC         AIC values for all fitted models
%     BIC         BIC values for all fitted models
%     AICc        AICc values for all fitted models
%
%   coFracProfile is an Nx1 vector of a single co-fractionation profile.
%   
%   x is an Nx1 vector of fraction numbers corresponding to coFracProfile.
%   Default is x = 1:length(coFracProfile).
%
%   aicstring is the string 'AIC' (default), 'BIC', or 'AICc' and determines 
%   the model selection method.
% 
%   Ngaussmax is an integer that determines the maximum complexity Gaussian 
%   mixture model (default 5).
%
%   See also GAUSSBUILD.



% 0. Initialize

if nargin<1
  error('Must input a chromatogram.')
end

if nargin<2
  x = 1:length(coFracProfile);
  x = x';
end
if isempty(x)
  x = 1:length(coFracProfile);
  x = x';
end

if nargin<3
  aicstring = 'AICc';
end

if nargin<4
  Ngaussmax = 5;
end

[n1,n2] = size(coFracProfile);
if n1~=1 && n2~=1 || (n1==1 && n2==1)
  error('Chromatogram must be a 1-dimensional vector.')
end

% transform chromatogram into column vector if necessary
if n2>n1
  coFracProfile = coFracProfile';
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
LB = [-inf -inf -inf]; % H, C, W
UB = [inf length(coFracProfile)+10 inf]; % H, C, W
fo{1} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,1),...
  'Upper',repmat(UB,1,1),'MaxIter',400,'Display','off');
fo{2} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,2),...
  'Upper',repmat(UB,1,2),'MaxIter',400,'Display','off');
fo{3} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,3),...
  'Upper',repmat(UB,1,3),'MaxIter',400,'Display','off');
fo{4} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,4),...
  'Upper',repmat(UB,1,4),'MaxIter',400,'Display','off');
fo{5} = fitoptions('method','NonlinearLeastSquares','Robust','LAR','Lower',repmat(LB,1,5),...
  'Upper',repmat(UB,1,5),'MaxIter',400,'Display','off');


I = coFracProfile>0;
msse = nan(5,1);
curveFit = cell(5,1);
R2 = nan(5,1);
for ii=1:Ngaussmax
  iter = 0;
  fitAgain = 1;
  while fitAgain
    iter = iter+1;

    ics = makeics(coFracProfile,ii); % make initial conditions
    set(fo{ii},'startpoint',ics(1:3*ii)); % Include ICs in fit options
    curveFit{ii} = fit(x,coFracProfile,ft{ii},fo{ii}); 
    yhat = feval(curveFit{ii},x);
    msse(ii) = sum((yhat(I) - coFracProfile(I)).^2);
    tmp = corrcoef(yhat,coFracProfile);
    R2(ii) = tmp(1,2).^2;
    
    %fitAgain = R2(ii)<0.5 & iter<=MaxIter;
    cf = coeffvalues(curveFit{ii});
    fitAgain = sum(cf(1:3:end))<0.5 & iter<=MaxIter;
  end
end



%% 2. Pick the best model (AIC / AICc / BIC)
% AIC = (n)log(SSE/n)+2p
% BIC = (n)log(SSE/n)+(p)log(n)

%N = length(coFracProfile);
N = sum(coFracProfile>0);
AIC = ones(5,1) * 10^9;
AICc = ones(5,1) * 10^9;
BIC = ones(5,1) * 10^9;
for ii=1:Ngaussmax
  k = length(coeffvalues(curveFit{ii}));
  AIC(ii) = N*log(msse(ii)/N) + 2*k;
  AICc(ii) = N*log(msse(ii)/N) + 2*k + 2*k*(k+1)/(N-k-1);
  BIC(ii) = N*log(msse(ii)/N) + k*log(N);
end

% Throw out any models where the number of non-imputed data points, i.e. the real data, is less than
% the number of parameters.
% Ngood = sum(coFracProfile>0.01);
% maxModelSize = floor(Ngood/3);
% AIC(maxModelSize:end) = inf;
% AICc(maxModelSize:end) = inf;
% BIC(maxModelSize:end) = inf;

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
model.CIs = confint(curveFit{Ngauss});
model.SSE = msse(Ngauss);
model.adjrsquare = R2(Ngauss);
model.curveFit = curveFit{Ngauss};
model.AIC = AIC;
model.AICc = AICc;
model.BIC = BIC;

% Throw out Gaussians whose height is less than 15% of the max?
if 1
  % height < 15% of max | Center is out of bounds
  I = find(model.coeffs(1:3:end) < max(coFracProfile)*0.15 ...
      | model.coeffs(2:3:end)<x(1) | model.coeffs(2:end:end)>x(end)) * 3 - 2;
  I2 = zeros(length(I)*3,1);
  for ii = 1:length(I)
    I2((ii-1)*3 + 1 : ii*3) = I(ii) + [0 1 2];
  end
  model.coeffs(I2) = [];
  model.CIs(:,I2) = [];
  Ngauss = Ngauss - length(I);
  xfit = x;
  yfit = zeros(size(xfit));
  for kk = 1:Ngauss
    H = model.coeffs(1 + (kk-1)*3);
    C = model.coeffs(2 + (kk-1)*3);
    W = model.coeffs(3 + (kk-1)*3);
    yfit = yfit + H*exp(-((xfit-C)/W).^2);
  end
  tmp = corrcoef(yfit,coFracProfile);
  model.adjrsquare = tmp(1,2).^2;
end

model.Ngauss = Ngauss;
if model.Ngauss>0
  model.fo = fo{Ngauss};
  model.ft = ft{Ngauss};
else
  model.fo = [];
  model.ft = [];
end



function ics = makeics(coFracProfile,Ngauss)

% Make initial conditions, ics
[pks, pksI] = findpeaks([0;coFracProfile;0],[0 1:length(coFracProfile)+1]');
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
    ics((ii-1)*3 + 1) = max(coFracProfile) + rand-.5; % A, height
    ics((ii-1)*3 + 2) = rand*length(coFracProfile) + rand*3-1.5; % mu, center
    ics((ii-1)*3 + 3) = ss + rand-.5; % sigma, width
  end
end

ics(Ngauss*3+1 : end) = 0;

