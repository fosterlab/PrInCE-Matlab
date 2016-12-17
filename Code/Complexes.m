%COMPLEXES Predict protein complexes
%   Using the ClusterONE algorithm, COMPLEXES builds protein complexes from
%   the protein-protein interactions predicted by INTERACTIONS. Two free
%   parameters (p, dens) are optimized through a grid search in order to
%   predict complexes that most closely match the reference complexes, as
%   quantified by their matching ratio. Predicted complexes that have at
%   least two proteins in common with one or more known complexes are
%   matched to these known complexes and reported in output tables.
%
%   Master script PRINCE and modules GAUSSBUILD and INTERACTIONS must be run 
%   before COMPLEXES.
%
%   COMPLEXES produces two output folders: Output/Data/Complexes, which
%   contains all output csv tables, and Output/Figures/Complexes, which
%   contains all figures.
%
%   Workflow:
%   1. Read Interactions output
%   2. Pre-process
%   3. Optimize parameters p and dens
%   4. Predict complexes
%   5. Match each predicted complex to a reference complex, if possible
%   6. Write output tables
%   7. Make figures
%
%   See also PRINCE, GAUSSBUILD, INTERACTIONS, MYCLUSTERONE, MATCHINGRATIO.

%   References:
%      Tamas Nepusz, Haiyuan Yu, and Alberto Paccanaro. Detecting 
%         overlapping protein complexes in protein-protein interaction 
%         networks. Nature methods, 9(5):471-472, 2012.



% To fix:
%
% - Because geometric accuracy doesn?t penalize novel complex members, it can lead to very large complexes. So matching ratio is better?
% - Don?t worry about optimizing best_p
% - Do worry about optimizing the intMatrix. The best one so far is
%       intMatrix2 = intMatrix-0.75;
%       intMatrix2(intMatrix2<0) = 0;
%       intMatrix2 = intMatrix2.^10;
%   i.e. relatively high precision (75%), cutoff-subtracted, and a high exponent.
%
% To Do
% 1. Fix optimization.
% 2. Make sure that no-isoform versions are matched to CORUM.


%% 0. Initialize
tic
fprintf('\n    0. Initialize')

minrep = user.minrep; % minimum number of replicates an interaction has to be in
desiredPrecision = user.desiredPrecision;
Nchannels = length(user.silacratios);

% Define folders, i.e. define where everything lives.
datadir = [user.maindir '/Output/Data/Complexes/']; % where data files live
figdir = [user.maindir '/Output/Figures/Complexes/']; % where data files live
% Make folders if necessary
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist(figdir, 'dir'); mkdir(figdir); end

% pairwise interactions, found in Final_Interactions_list
s1 = ['Final_Interactions_list_' num2str(desiredPrecision*100) '_precision.csv'];
s = [user.maindir 'Output/Data/Interactions/' s1];

countPrec = 0;
if ~exist(s,'file')
  % desiredPrecision was likely reset in ROC_PCPSILAC
  % look for other Final_Interactions_lists
  dd = dir([user.maindir 'Output/Data/Interactions/Final_Interactions_list_*_precision.csv']);
  dd_name_length = nan(length(dd),1);
  for jj = 1:length(dd)
    dd_name_length(jj) = length(dd(jj).name);
  end
  if sum(dd_name_length==40)>0
    I = find(dd_name_length==40,1,'last');
    warning('\n    Complexes: Could not find %s',['Final_Interactions_list_' num2str(desiredPrecision(ii)*100) '_precision.csv'])
    warning('\n    Using %s instead. Continuing...\n\n',dd(I).name)
    s = [user.maindir 'Output/Data/Interactions/' dd(I).name];
    countPrec = countPrec+1;
    InteractionIn = s;
    
    % change desiredPrecision to reflect the file you're using
    desiredPrecision(ii) = str2double(s(end-15:end-14)) / 100;
  else
    error('\n    Following interaction file not found:\n        %s\n    No replacement Final_Interactions_list found. Aborting...',s1)
  end
else
  countPrec = countPrec+1;
  InteractionIn = s;
end
if countPrec == 0
  fprintf('\n    Error: Complexes: No interaction files found!')
end

% scoreMatrix+Protein, found in score_repx.mat
scoreIn = cell(length(user.Nreplicate),1);
countScore = 0;
for ii = 1:length(user.silacratios)
  fn = ['score_chan' num2str(ii) '.mat'];
  s = [user.maindir 'Output/tmp/' fn];
  
  if ~exist(s,'file')
    fprintf('\n    Error: Complexes: Following scoreMatrix file not found:')
    fprintf('\n        %s\n',s)
  else
    countScore = countScore+1;
    scoreIn{ii} = s;
  end
end
if countScore < length(user.silacratios)
  fprintf('\n    Error: Complexes: Missing a scoreMatrix file!')
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)




%% 1. Read input data
% Final_Interactions_list
% scoreMatrix+Protein
% corum
tic
fprintf('    1. Read input data')

% Final_Interactions_list
[tmp,tmp_head] = readFinalInteractionList(InteractionIn);

% Need 4 columns: Protein-A, Protein-B, channel, precision-cutoff
Itext(1) = find(ismember(tmp_head.text,'Protein A'));
Itext(2) = find(ismember(tmp_head.text,'Protein B'));
Itext(3) = find(ismember(tmp_head.text,'Channel'));
Idata = find(ismember(tmp_head.data,'Precision level') | ...
    ismember(tmp_head.data,'Precision (avg.)'));
interactionPairs = [tmp.text(:,Itext) num2cell(tmp.data(:,Idata))];
for ii = 1:size(interactionPairs,1)
  I1 = find(ismember(interactionPairs{ii,1},'-'));
  I2 = find(ismember(interactionPairs{ii,2},'-'));
  % Remove isoform tags
  if ~isempty(I1)
    interactionPairs{ii,1} = interactionPairs{ii,1}(1:I1-1);
  end
  if ~isempty(I2)
    interactionPairs{ii,2} = interactionPairs{ii,2}(1:I2-1);
  end
  interactionPairs(ii,1:2) = sort(interactionPairs(ii,1:2));
end

clear tmp tmp_head

% Protein structure
protIndex = cell(length(scoreIn),1);
protAll = cell(length(scoreIn),1);
protAll_noIsoforms = cell(size(protAll));
for ii = 1:length(scoreIn)
  % score, proteinAll, possList, classList, indexList, feats_new
  load(scoreIn{ii})
  protIndex{ii} = indexList;
  protAll{ii} = proteinAll;
  
  protAll_noIsoforms{ii} = cell(size(protAll{ii}));
  for jj = 1:length(protAll{ii})
    I1 = find(ismember(protAll{ii}{jj},'-'));
    if isempty(I1)
      I1 = length(protAll{ii}{jj})+1;
    end
    protAll_noIsoforms{ii}{jj} = protAll{ii}{jj}(1:I1-1);
  end
  clear score proteinAll possList classList feats_new
end

% Corum complexes
corumComplex = importdata(user.corumcomplexfile, ',');
if ~isvector(corumComplex)
  error('Error: Complexes: Incorrectly formatted CORUM complex file')
end

% Corum pairwise interactions
fid=fopen(user.corumpairwisefile, 'rt');    %Input corum data base.
Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
fclose(fid);
No=length(Corum_Import{1})/2;
Corum_Protein_names = (reshape(Corum_Import{1,1},2,No)');
Unique_Corum = unique(Corum_Protein_names);


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 2. Pre-process pairwise interactions

tic
fprintf('    2. Pre-process pairwise interactions')

% Make uniqueProteins, a single list of all proteins
uniqueProteins = unique(cat(2,protAll_noIsoforms{:})); % count the number of proteins in interactions
Nprot_pred = length(uniqueProteins);
uniqueProteins = unique([uniqueProteins Unique_Corum']); % then include proteins from reference
Nproteins = length(uniqueProteins);

% Convert protein IDs to indices
% (numerical values are faster to work with!)
% interactionPairs --> interactionPairs2, indices of uniqueProteins
interactionPairs2 = nan(size(interactionPairs,1),6);
channels = nan(size(interactionPairs,1),length(user.silacratios));
for jj = 1:size(interactionPairs,1)
  I1 = find(ismember(uniqueProteins,interactionPairs{jj,1}));
  I2 = find(ismember(uniqueProteins,interactionPairs{jj,2}));
  
  % what channels was this interaction seen in?
  chans = interactionPairs{jj,4};
  for kk = 1:length(user.silacratios)
    if ismember(user.silacratios{kk}, chans)
      channels(jj,kk) = kk;
    end
  end
  
  tmp2 = zeros(1,3);
  for kk = 1:3
    tmp2(kk) = double(ismember(num2str(kk),interactionPairs{jj,3}));
  end
  
  precdrop = interactionPairs{jj,4};
  
  interactionPairs2(jj,:) = [I1 I2 precdrop tmp2];
end

corumComplex2 = cell(size(corumComplex));
for ii = 1:length(corumComplex)
  cmplx = corumComplex{ii};
  
  Idelim = [0 strfind(cmplx, ',') length(cmplx)+1];
  Nprot = length(Idelim) - 1;
  
  I = nan(1, Nprot);
  for jj = 1:Nprot
    prot1 = cmplx(Idelim(jj)+1 : Idelim(jj+1)-1);
    I(jj) = find(ismember(uniqueProteins, prot1));
  end
  
  if sum(isempty(I))>0 || sum(isnan(I))>0
    error('Error: Complexes: CORUM complex protein not found in pairwise proteins.')
  end
  
  % Check that complex includes at least one predicted-interaction protein.
  if sum(I<=Nprot_pred)>0
    corumComplex2{ii} = I;
  end
end
corumComplex2(cellfun('isempty',corumComplex2)) = [];

% Determine how to split up files
% column 1 = replicate index
% column 2 = channel index
% 0 indicates use all replicates / channels.
if user.separateByReplicate == 1 && Nreplicates ==1
  warning('user.separateByReplicate set to 1, but only one replicate detected.')
  warning('Setting user.separateByReplicate to 0 and proceeding...')
  user.separateByChannel = 0;
end
if user.separateByChannel == 1 && Nchannels ==1
  warning('user.separateByChannel set to 1, but only one channel detected.')
  warning('Setting user.separateByChannel to 0 and proceeding...')
  user.separateByChannel = 0;
end
clear csplit
if user.separateByReplicate == 0 && user.separateByChannel == 0
  % All interactions
  csplit = [0 0];
elseif user.separateByReplicate == 1 && user.separateByChannel == 0
  % All interactions + per-replicate
  csplit = nan(Nreplicates, 2);
  cc = 0;
  for ii = 1:Nreplicates
    cc = cc+1;
    csplit(cc,:) = [ii 0];
  end
  csplit = [csplit; 0 0];
elseif user.separateByReplicate == 0 && user.separateByChannel == 1
  % All interactions + per-channel
  csplit = nan(Nchannels, 2);
  cc = 0;
  for jj = 1:Nchannels
    cc = cc+1;
    csplit(cc,:) = [0 jj];
  end
  csplit = [csplit; 0 0];
elseif user.separateByReplicate == 1 && user.separateByChannel == 1
  % All interactions + per-replicate + per-channel
  csplit = nan(Nchannels * Nreplicates, 2);
  cc = 0;
  for ii = 1:Nreplicates
    for jj = 1:Nchannels
      cc = cc+1;
      csplit(cc,:) = [ii,jj];
    end
  end
  csplit = [csplit; 0 0];
else
  error('Badly formatted user.separateByReplicates, user.separateByChannel')
end
csplit = unique(csplit,'rows');

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


%% 3. Optimize complex-building parameters
% For each replicate, optimize:
%   i) score cutoff
%   ii) p, the uncertainty parameter in clusterONE
%   iii) density threshold, how "loose" is too loose for clusterONE?
tic
fprintf('    3. Optimize parameters for complex-building')

pRange = [0 100 500 2000 10000 50000 10^6];
densRange = [0 0.1 0.2 0.3 0.4 0.5 0.75 1];

ga = nan(length(pRange),length(densRange));
mr = nan(length(pRange),length(densRange));
dens = nan(length(pRange),length(densRange));

best_p = nan(1,size(csplit,1));
best_dens = nan(1,size(csplit,1));
user.optimizeHyperParameters = 0;
if user.optimizeHyperParameters ==1
  for ii = 1:size(csplit,1)
    
    if csplit(ii,1) == 0
      I1 = ones(size(interactionPairs2,1),1);
    else
      I1 = sum(replicates == csplit(ii,1),2)>0;
    end
    if csplit(ii,2) == 0
      I2 = ones(size(interactionPairs2,1),1);
    else
      I2 = sum(channels == csplit(ii,2),2)>0;
    end
    I = I1==1 & I2==1;
    ip2_subset = interactionPairs2(I,:);
    
    
    % This makes intMatrix from interactionPairs2.
    intMatrix = zeros(Nprot,Nprot);
    for jj = 1:size(ip2_subset,1)
      x = ip2_subset(jj,1:2);
      nrep = ip2_subset(jj,3);
      intMatrix(x(1),x(2)) = intMatrix(x(1),x(2)) + nrep;  % add nrep
      intMatrix(x(2),x(1)) = intMatrix(x(2),x(1)) + nrep;
    end
    %intMatrix(intMatrix>0) = 1;
    
    for jj = 1:length(pRange)
      s = ['\n    p=' num2str(pRange(jj))];
      fprintf(s)
      
      % ii) make complexes
      [Members, Density] = myclusterone(intMatrix, pRange(jj), 0);
      if isempty(Members); continue; end
      
      for mm = 1:length(densRange)
        %s = ['\n        densThresh=' num2str(densRange(mm))];
        %fprintf(s)
        
        Members2 = Members;
        
        % Remove small or low-density complexes
        dens_comp = nan(size(Members2));
        for kk = 1:size(Members2)
          I = Members2{kk};
          m = intMatrix(I,I);
          n = length(I);
          dens_comp(kk) = sum(m(:)) / (n * (n-1)/2);
          if n<3 || dens_comp(kk) < densRange(mm)
            Members2{kk} = [];
          end
        end
        Members2(cellfun('isempty',Members2)) = [];
        if isempty(Members2); continue; end
        
        % iv) assess how good the complexes are
        %ga(jj,mm) = geomacc(Members2,corumComplex2);
        mr(jj,mm) = matchingratio(Members2,corumComplex2);
        %dens(jj,mm) = mean(Density);
        
      end
    end
    
    % va) Make matrix to optimize on
    %ga = (ga - min(ga(:))) / (max(ga(:)) - min(ga(:)));
    %mr = (mr - min(mr(:))) / (max(mr(:)) - min(mr(:)));
    %dens = (dens - min(dens(:))) / (max(dens(:)) - min(dens(:)));
    optMatrix = mr;
    
    % vb) choose parameters
    [~,I] = max(optMatrix(:));
    [i1, i2] = ind2sub(size(optMatrix),I);
    best_p(ii) = pRange(i1);
    best_dens(ii) =  densRange(i2);
    
  end
else
  fprintf(' ... skipping optimization ')
  best_p = 500 * ones(size(csplit,1),1);
  best_dens = 0.1 * ones(size(csplit,1),1);
end


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 4. Build final complex list

tic
fprintf('    4. Build final complex list')

Ncomplex = zeros(countPrec,1);
for ii = 1:size(csplit,1)
  
  if csplit(ii,1) == 0
    I1 = ones(size(interactionPairs2,1),1);
  else
    I1 = sum(replicates == csplit(ii,1),2)>0;
  end
  if csplit(ii,2) == 0
    I2 = ones(size(interactionPairs2,1),1);
  else
    I2 = sum(channels == csplit(ii,2),2)>0;
  end
  I = I1==1 & I2==1;
  ip2_subset = interactionPairs2(I,:);
  
  % This makes intMatrix from interactionPairs2.
  intMatrix = zeros(Nproteins,Nproteins);
  for jj = 1:size(ip2_subset,1)
    x = ip2_subset(jj,1:2);
    nrep = ip2_subset(jj,3);
    intMatrix(x(1),x(2)) = nrep;  % add nrep
    intMatrix(x(2),x(1)) = nrep;
  end
  intMatrix = intMatrix - user.desiredPrecision;
  intMatrix(intMatrix<0) = 0;
  %intMatrix = intMatrix-0.75;
  %intMatrix(intMatrix<0) = 0;
  %intMatrix = intMatrix.^10;
  
  % Make complex list
  [CL(ii).Members,CL(ii).Density] = myclusterone(intMatrix, best_p(ii), 0);
  
  % Remove small or low-density complexes
  tmp = size(CL(ii).Members);
  if tmp(2)==0
    NN = 0;
  else
    NN = tmp(1);
  end
  for kk = 1:NN
    I = CL(ii).Members{kk};
    m = intMatrix(I,I);
    n = length(I);
    dens_comp = sum(m(:)) / (n * (n-1)/2);
    if n<3 || dens_comp<best_dens(ii)
      CL(ii).Members{kk} = [];
    end
  end
  I = cellfun('isempty',CL(ii).Members);
  CL(ii).Members(I) = [];
  CL(ii).Density(I) = [];
  NN = sum(~I);
  
  % make mini-interaction matrices
  for jj = 1:NN
    CL(ii).Connections{jj} = zeros(length(CL(ii).Members{jj}),length(CL(ii).Members{jj}));
    for kk = 1:length(CL(ii).Members{jj})
      for mm = 1:length(CL(ii).Members{jj})
        CL(ii).Connections{jj}(kk,mm) = intMatrix(CL(ii).Members{jj}(kk),CL(ii).Members{jj}(mm));
      end
    end
  end
  CL(ii).GeomAcc = geomacc(CL(ii).Members,corumComplex2);
  CL(ii).MatchRatio = matchingratio(CL(ii).Members,corumComplex2);
  Ncomplex(ii) = length(CL(ii).Members);
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 5. Match each predicted complex to a CORUM complex
% for each predicted complex
%   calculate the overlap with each CORUM complex
%   report the CORUM complex with the highest overlap

tic
fprintf('    5. Match each complex to a CORUM complex')

corumMatches = cell(countPrec,1);
for ii = 1:size(csplit,1)
  
  % columns: 1. predicted complex number (Members),
  %          2. corum complex number (corumComplex2),
  %          3. match rank (e.g. best, 2nd best),
  %          4. N overlap,
  %          5. overlap ratio (N_overlap/size_of_corum)
  corumMatches{ii} = zeros(Ncomplex(ii)*5, 2);
  cc = 0;
    
  overlap = zeros(Ncomplex(ii),length(corumComplex2));
  overlap_ratio = zeros(Ncomplex(ii),length(corumComplex2));
  for jj = 1:Ncomplex(ii)
    for kk = 1:length(corumComplex2)
      overlap(jj,kk) = length(intersect(CL(ii).Members{jj},corumComplex2{kk}));
      overlap_ratio(jj,kk) = length(intersect(uniqueProteins_noIsoform(CL(ii).Members{jj}),...
        uniqueProteins_noIsoform(corumComplex2{kk}))) / length(corumComplex2{kk});
    end
    overlap_ratio(overlap<2) = 0;
    overlap(overlap<2) = 0;
    
    % sort first by N overlap
    Nover = sort(unique(overlap(jj,:)), 'descend');
    Nover = Nover(Nover>0);
    for mm = 1:min(2,length(Nover))
      % in case of ties, sort by overlap_fraction
      lapratio = overlap_ratio(jj,overlap(jj,:) == Nover(mm));
      [~,I] = sort(lapratio,'descend');
      matches = find(overlap(jj,:) == Nover(mm));
      matches = matches(I);
      for nn = 1:length(matches)
        cc = cc+1;
        corumMatches{ii}(cc,1) = jj;
        corumMatches{ii}(cc,2) = matches(nn);
        corumMatches{ii}(cc,3) = mm;
        corumMatches{ii}(cc,4) = length(intersect(CL(ii).Members{jj},corumComplex2{matches(nn)}));
        corumMatches{ii}(cc,5) = overlap_ratio(jj,matches(nn));
      end
    end
  end
  corumMatches{ii} = corumMatches{ii}(1:cc,:);
  
end


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 6. Write output

tic
fprintf('    6. Write output')

writeOutput_complexes

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 7. Make figures

tic
fprintf('    7. Make figures')

makeFigures_complexes

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


