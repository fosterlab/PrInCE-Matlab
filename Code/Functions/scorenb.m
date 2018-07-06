function [score,feats] = scorenb(Dist,possibleInts,TP_Matrix)
%SCORENB Predict protein interactions
%   [prob,feats]=SCORENB(Dist,P,TPFP) calculates the probability that
%   protein-protein pairs are interactions. SCORENB requires two inputs:
%   features describing the similarity of pairs co-fractionation profiles
%   (Dist), and interacting or non-interacting labels for a subset of
%   protein pairs (P, TPFP). Dist is an NxM  matrix of M features ("distance
%   measures") for N protein-protein pairs. P is an Nx1 logical vector 
%   corresponding to Dist. A value of P=1 denotes a protein pair where both 
%   proteins are in the reference database, and P=0 denotes a protein pair
%   where one or both proteins are not in the reference database. TPFP is
%   an Nx1 logical vector corresponding to Dist, where TPFP=1 labels a
%   gold-standard interaction, and TPFP=0 & P=0 labels a gold-standard
%   non-interaction. Therefore, TPFP==1 & P==1 denotes rows in Dist that
%   are known interactions, and TPFP==0 & P==1 denotes rows in Dist that
%   are known non-interactions.
%
%   prob is the Nx1 vector of interaction probabilities.
%
%   Features (distance measures) are selected using the criterion Fisher
%   Ratio > 2. feats is a 15xM matrix of Fisher Ratio values used for
%   feature selection, where each row is one iteration of the 15-fold
%   cross-validation and M is the number of features.
%
%   See also INTERACTIONS, INDFEAT.


% Use a Naive Bayes classifier to predict protein interactions.
% Can take NxN square matrices OR a 2D table.

% Get data
if isstruct(Dist)
  % Dist is a structure with NxN fields
  fn = fieldnames(Dist);
  X = nan(length(Dist.(fn{1})(:)),length(fn));
  for ii = 1:length(fn)
    X(:,ii) = Dist.(fn{ii})(:);
  end
  I = possibleInts(:)==1;
else
  % Dist is a table
  X = Dist;
  I = possibleInts(:)==1;
end

y = TP_Matrix(:);
y(y>0) = 1;
y(y~=1) = -1;
Nd = size(X,2);

% Impute missing values
for ii = 1:size(X,2)
  im = isnan(X(:,ii)); % rows with missing values
  X(im & y==1,ii) = randsample(X(~im & y==1,ii),sum(im & y==1),1);
  X(im & y==-1,ii) = randsample(X(~im & y==-1,ii),sum(im & y==-1),1);
end

% Soft whiten data
eps = 2e-16;
wmu = zeros(Nd,1);
wstd = zeros(Nd,1);
for ii = 1:Nd
  wmu(ii) = nanmean(X(:,ii));
  wstd(ii) = nanstd(X(:,ii));
  X(:,ii) = (X(:,ii) - wmu(ii)) / 2 / (wstd(ii) + eps);
end

% Divide labelled data into folds
Npos = sum(y(I)==1);
Nneg = sum(y(I)==-1);
Ilabel = find(I);
if Npos + Nneg <= 15
    ss = ['Few training labels (' num2str(Npos+Nneg) '). Gold standard reference may be too small.'];
    disp(ss)
    disp('Using Leave One Out cross-validation...')
    
    Ifold = randperm(length(Ilabel));
    Nfold = length(Ilabel);
else
    Nfold = 15;
    Ifold = ceil( randperm(length(Ilabel)) / length(Ilabel) * Nfold );
end


score = nan(size(X,1), Nfold);
feats = nan(Nfold, size(X,2));
for iter = 1:Nfold
  %non-balanced training data
  Itest = Ilabel(Ifold==iter);
  Itest = ismember(1:length(y),Itest);
  Itest = Itest | I'==0;
  Itrain = ~Itest;
  
  Xtr = X(Itrain,:);
  ytr = y(Itrain);
  Xnew = X(Itest,:);
  
  % Ensure no class-variable pair has zero variance
  testvar = nan(2,Nd);
  for jj = 1:Nd
    testvar(1,jj) = nanvar(Xtr(ytr>0,jj));
    testvar(2,jj) = nanvar(Xtr(ytr<0,jj));
  end
  
  % Feature selection
  feats(iter,:) = IndFeat(Xtr,ytr);
  if sum(feats(iter,:)<=2)==size(feats,2)
    feats(iter,:) = 3;
  end
  f2consider = find(feats(iter,:) > 2 & testvar(1,:)>0 & testvar(2,:)>0 );
  if isempty(f2consider)
    f2consider = find(testvar(1,:)>0 & testvar(2,:)>0 );
  end
  if isempty(f2consider)
    warning('Not enough data to classify (all features have zero variance).')
    continue
  end
  
  % Fit Naive Bayes model
  nab = fitcnb(Xtr(:,f2consider),ytr);
  [~,scoretmp] = predict(nab,Xnew(:,f2consider));
  
  score(Itest,iter) = scoretmp(:,2);
end
