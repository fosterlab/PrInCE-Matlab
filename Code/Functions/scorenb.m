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
  for ii = 1:length(fn);
    X(:,ii) = Dist.(fn{ii})(:);
  end
  I = possibleInts(:)==1;
else
  % Dist is a table
  X = Dist;
  I = possibleInts(:)==1 & sum(isnan(Dist),2)==0;
end

y = TP_Matrix(:);
y(y>0) = 1;
y(y~=1) = -1;
Nd = size(X,2);

% Impute missing values
for ii = 1:size(X,2)
 I_missing = isnan(X(:,ii));
 MU = nanmean(X(:,ii));
 SD = nanstd(X(:,ii));
 X(I_missing,ii) = nanmedian(X(:,ii)) + normrnd(MU, SD*.05);
end

% what training length do we need to expect to get 5 class == 1?
classratio = sum(y==1) / sum(y==-1);
trainingLength = max([8000 round(5/classratio)]);
trainingLength = min([trainingLength ceil(sum(I)/2)]);

Nmodel = 15;
score = nan(size(X,1),Nmodel);
feats = nan(Nmodel,size(X,2));
for iter = 1:Nmodel
  % Make training and testing data
  
  if 0
    % balance training data
    I1 = find(y==1 & I);
    I1 = I1(randsample(length(I1),round(trainingLength/2)));
    I0 = find(y==-1);
    I0 = I0(randsample(length(I0),round(trainingLength/2)));
    Itrain = ismember(1:length(y),[I1;I0]);
    Ipred = ~Itrain;
    Xtr = X(Itrain,:);
    ytr = y(Itrain);
    Xnew = X(Ipred,:);
    
    % Ensure no class-variable pair has zero variance
    testvar = nan(2,Nd);
    for jj = 1:Nd
      testvar(1,jj) = nanvar(Xtr(ytr>0,jj));
      testvar(2,jj) = nanvar(Xtr(ytr<0,jj));
    end
    
  else
    % non-balanced training data
    go = 1;
    train_iter = 0;
    iterMax = 1000;
    while go
      train_iter = train_iter+1;
      Iall = find(I);
      Iall = Iall(randsample(length(Iall),trainingLength));
      Itrain = ismember(1:length(y),Iall);
      
      Ipred = ~Itrain;
      Xtr = X(Itrain,:);
      ytr = y(Itrain);
      Xnew = X(Ipred,:);
      
      % Ensure there are at least 3 data points in each class
      go1 = sum(ytr>0)<2 | sum(ytr<0)<2;
      
      % Ensure no class-variable pair has zero variance
      testvar = nan(2,Nd);
      for jj = 1:Nd
        testvar(1,jj) = nanvar(Xtr(ytr>0,jj));
        testvar(2,jj) = nanvar(Xtr(ytr<0,jj));
      end
      go2 = sum(testvar(1,:)~=0) == 0 | sum(testvar(2,:)~=0) == 0;
      
      go = (go1 | go2) & train_iter<=iterMax;
    end
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
  
  % Fit Naive Bayes model
  nab = fitcnb(Xtr(:,f2consider),ytr);
  [~,scoretmp] = predict(nab,Xnew(:,f2consider));
  
  score(Ipred,iter) = scoretmp(:,2);
end

