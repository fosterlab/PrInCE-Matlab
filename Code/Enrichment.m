


%% 0. Initialize
tic
fprintf('\n    0. Initialize')

% Define folders, i.e. define where everything lives.
datadir = [user.maindir 'Data/']; % where data files live
figdir = [user.maindir 'Figures/Enrichment/']; % where figures live
%tmpdir = '/Users/Mercy/Academics/Foster/Jenny_PCPSILAC/PCPSILAC_Analysis/Data/Alignment/';
% Make folders if necessary
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist([datadir '/Enrichment/'], 'dir'); mkdir([datadir '/Enrichment/']); end
if ~exist(figdir, 'dir'); mkdir(figdir); end


InteractionIn = cell(length(user.desiredPrecision),1);
countPrec = 0;
for ii = 1:length(user.desiredPrecision)
  s1 = ['Final_Interactions_list_' num2str(user.desiredPrecision(ii)*100) '_precision.csv'];
  s = [user.maindir 'Data/ROC/CombinedResults/' s1];
  
  if ~exist(s,'file')
    fprintf('\n    Error: Enrichment: Following interaction file not found:')
    fprintf('\n        %s\n',s1)
  else
    countPrec = countPrec+1;
    InteractionIn{countPrec} = s;
  end
end
if countPrec == 0
  fprintf('\n    Error: Enrichment: No interaction files found!')
end
InteractionIn = InteractionIn(1:countPrec);


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)




%% 1. Read input data
tic
fprintf('    1. Read input data')

% Load interactions detected by ROC_PCPSILAC
interactionPairs = cell(countPrec,1);
for ii = 1:countPrec
  tmp = importdata(InteractionIn{ii});
  
  % remove header from textdata if necessary
  if size(tmp.textdata,1)==size(tmp.data,1)+1
    tmp.textdata = tmp.textdata(2:end,:);
  elseif size(tmp.textdata,1)==size(tmp.data,1)
  else
    disp('Error: Enrichment: Mismatch between imported numerical and text data sizes')
  end
  
  % only need the Protein-A, Protein-B, and number-of-replicates columns
  % Hard code these to columns 2 and 3
  interactionPairs{ii} = tmp.textdata(:,[2:3 6]);
end

% List of N unique protein names
uniqueProteins = [];
for ii = 1:countPrec
  tmp = interactionPairs{ii}(:,1:2);
  uniqueProteins = unique([uniqueProteins; tmp(:)]);
end
N = length(uniqueProteins);

% Load CORUM interactions
corumPairs_tmp = importdata(user.corumfile);
corumPairs_tmp = unique(corumPairs_tmp);

% Check that corumPairs is an mx1 cell
if ~iscell(corumPairs_tmp) || size(corumPairs_tmp,2)~=1
  error('Error: Enrichment: Incorrectly formatted CORUM file')
end

% split corumPairs into a mx2 cell
% this method is surprisingly faster than cellfun(@strsplit)
corumPairs = cell(length(corumPairs_tmp),2);
for ii = 1:length(corumPairs_tmp)
  fn = corumPairs_tmp{ii};
  for jj = 1:length(fn)
    if fn(jj) == ','
      corumPairs{ii,1} = fn(1:jj-1);
      corumPairs{ii,2} = fn(jj+1:end);
    end
  end
end
clear corumPairs_tmp

% Turn interactionPairs,corumPairs into indices of uniqueProteins
% (numerical pairs are easier to work with!)
% In interaction index pairs, also include how many replicates it was seen in.
interactionPairs2 = cell(countPrec,1);
for kk = 1:countPrec
  interactionPairs2{kk} = nan(size(interactionPairs{kk},1),3);
end
corumPairs2 = nan(size(corumPairs));
for ii = 1:length(uniqueProteins)
  ii
  % corumPairs
  [Ic1,Ic2] = find(ismember(corumPairs,uniqueProteins{ii}));
  for jj = 1:length(Ic1)
    corumPairs2(Ic1(jj),Ic2(jj)) = ii;
  end
  
  % interactionPairs
  for kk = 1:countPrec
    [I1,I2] = find(ismember(interactionPairs{kk},uniqueProteins{ii}));
    for jj = 1:length(I1)
      interactionPairs2{kk}(I1(jj),I2(jj)) = ii;
      % include the number of replicates
      tmp = interactionPairs{kk}{I1(jj),3};
      nrep = length(unique(tmp(isstrprop(tmp,'digit'))));
      interactionPairs2{kk}(I1(jj),3) = nrep;
    end
  end
end

% corumPairs has both A-B and B-A interactions.
% Remove the redundant entries.
corumPairs2 = sort(corumPairs2,2);
[corumPairs2,I1,I2] = unique(corumPairs2,'rows');
% Use the same indices to remove redundant entries from the string list.
corumPairs = corumPairs(I1,:);


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



% %% 2. Build complexes
% tic
% fprintf('    2. Build complexes')
% 
% % Pre-allocate
% ComplexList.Members = cell(10000,1);
% ComplexList.Members{1} = 'FakeString-1';
% ComplexList.Connections = cell(10000,1);
% 
% 
% % Algorithm 1: Slow, for-loop
% for ii = 1:1%count
%   cmplxcount = 0;
%   for jj = 1:length(proteinPairs{ii})
%     protA = proteinPairs{ii}{jj,1};
%     protB = proteinPairs{ii}{jj,2};
%     
%     % Search for protA and protB in the complex members
%     I = ~cellfun('isempty',ComplexList.Members);
%     Ia = find(~cellfun('isempty',strfind(ComplexList.Members(I),protA)));
%     Ib = find(~cellfun('isempty',strfind(ComplexList.Members(I),protB)));
%     
%     
%     % 1.
%     % If Ia and Ib are not in complexes, make the complex protA-protB
%     if isempty(Ia) && isempty(Ia)
%       cmplxcount = cmplxcount+1;
%       ComplexList.Members{cmplxcount} = strjoin({protA,protB}, ', ');
%       ComplexList.Connections{cmplxcount} = [0 1; 1 0];
%     end
%     
%     % 2.
%     % If only one is in complexes, add the other protein to that complex
%     if length(Ia)==1 && isempty(Ib)
%       ComplexList.Members{Ia} = strjoin({ComplexList.Members{Ia},protB}, ', ');
%       
%       % what index is protA in the complex?
%       % count the number of commas before it, add one
%       I = strfind(ComplexList.Members{Ia}, protA);
%       Icomma = strmatch(',', ComplexList.Members{Ia});
%       idx = sum(Icomma < I) +1
%       
%       n = size(ComplexList.Connections{Ia},1);
%       tmp = zeros(n+1);
%       tmp(1:n,1:n) = ComplexList.Connections{Ia};
%       tmp(idx,n+1) = 1;
%       tmp(n+1,idx) = 1;
%       ComplexList.Connections{Ia} = tmp;
%     end
%     
%     % 3.
%     % If both are in one or more complexes, merge the complexes
%     if length(Ia)>1 || length(Ib)>1
%       
%     end
%     
%     % 4.
%     % If both are in the same complex, do nothing
% 
%     
%   end
% end
% ComplexList.Members = ComplexList.Members(1:cmplxcount);
% ComplexList.Connections = ComplexList.Connections(1:cmplxcount);
% 
% 
% tt = toc;
% fprintf('  ...  %.2f seconds\n',tt



%% 2. Build complexes
% Algorithm: 
% for each row of the interaction matrix
%   find interactions that are not in a complex yet
%   for each of these interactions
%     note that that interactions is in a complex
%     again, find indices in that row/colum in IM that are i) non-zero and ii) binary vector = 0
%     keep doing this until no new indices are found
%     assign all indices to a complex
%     make ComplexList.Members = indices, ComplexList.Connections = intMatrix(indices)
%
% Basically, start at first protein and explore all interactions. When no new interactions are
% found, call that a complex. Move on to the next protein not in a complex.
%
% One trick used here: interaction matrix values are number of replicates + 100*inCorum. This will
% give a values like 102 (in 2 replicates and Corum) and 3 (in 3 replicates, not corum).
% This is designed to use a single interaction matrix, since this can be very large for >3000
% proteins.
%
% Inputs: N1x2 interaction list, N2x2 CORUM list, list of N unique protein names, Nx1 binary
% Outputs: Complex(i).Members, Complex(i).Connections
tic
fprintf('    2. Build complexes')


clear ComplexList
for ii = 1:countPrec
  % Make NxN interaction matrix
  intMatrix = zeros(N,N);
  % add nrep
  for jj = 1:size(interactionPairs2{ii},1)
    x = interactionPairs2{ii}(jj,1:2);
    nrep = interactionPairs2{ii}(jj,3);
    intMatrix(x(1),x(2)) = intMatrix(x(1),x(2)) + nrep;
    intMatrix(x(2),x(1)) = intMatrix(x(2),x(1)) + nrep;
  end
  
  % add 100*inCorum
  I = find(~isnan(corumPairs2(:,1)) & ~isnan(corumPairs2(:,2)));
  for jj = 1:length(I)
    x = corumPairs2(I(jj),1:2);
    intMatrix(x(1),x(2)) = intMatrix(x(1),x(2)) + 100;
    intMatrix(x(2),x(1)) = intMatrix(x(2),x(1)) + 100;
  end
  
  % Reduce intMatrix to just upper-triangular
  intMatrix = triu(intMatrix);
  
  % Binary vector
  % Used to keep track of which proteins have already been assigned to a complex
  bv = zeros(N,1);

  % Start finding complexes
  countcmplx = 0;
  for jj = 1:N
    % check if this protein is already in a complex
    if bv(jj)~=0
      continue;
    end
    
    I = find(intMatrix(jj,:)>1 | intMatrix(:,jj)'>1);
    if ~isempty(I)
      % initialize this complex
      countcmplx = countcmplx+1;
      bv(jj) = 1;
      membs = [jj I];
      
      % BOOM! Another option is to IGNORE open_branches. Just keep exploring ALL branches, adding
      % interactions to the list, and pruning with unique. When the list stops growing, you've
      % stopped finding new interactions, and the complex is complete. Ignore the binary vector!
      
      % start exploring connections
      open_branches = membs(bv(membs)==0);
      while ~isempty(open_branches)
        for kk = 1:length(open_branches)
          I2 = find(intMatrix(open_branches(kk),:)>1 | intMatrix(:,open_branches(kk))'>1);
          membs = [membs I2];
        end
        membs = unique(membs);
        open_branches = membs(bv(membs)==0);
        bv(open_branches) = 1;
        length(open_branches)
      end
      
      ComplexList(ii).Members{countcmplx} = membs;
      ComplexList(ii).Connections{countcmplx} = intMatrix(membs,membs);
    end
  end
  
  % Check that all proteins are in a complex
  membcount = 0;
  for jj = 1:countcmplx
    membcount = membcount+length(ComplexList(ii).Members{jj});
  end
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



