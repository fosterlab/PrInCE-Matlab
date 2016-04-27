%%%%%%%%%%%%%%% Instructions for running Gauss_Build.m:
%
%
%
%%%%%%%%%%%%%%% Logic:
%   0. Initialize
%   for replicate
%       1. Read input
%       2. Pre-process Gaussians, make Chromatograms and Dist.
%       3. Make Protein structure
%       4. Make TP, FP, possibleInts
%       5. Score each interaction (Naive Bayes)
%   6. Calculate score threshold with all replicates
%   for desiredPrecision
%       7. Find and concatenate interactions at desired precision
%       8. Create a list of all the unique interactions
%       9. Calculate final TP, FP, FN and TN
%       10. Find Treatment interactions
%       11. Write output files
%       12. Make figures
%
%
%%%%%%%%%%%%%%% To do:
% - Why is this script pre-processing chromatograms? Shouldn't that have been done earlier?
% - What's going on in 8b? We keep interactions with Delta_Center>2... so why check for it?
% - FN_binary_interaction_list can get up to 3 GB. Can we make it smaller?
%
%
%%%%%%%%%%%%%%% Fix the logic:
% - 'Look up Major Protein group Information' can be seriously improved with a few findstr
% - FN and TN are approximated. Instead of looking at all the combinations of Gaussians, it just
%   treats them all the same. This is for speed reasons: there were too many combinations to go
%   through individually. A fix is to treat things exactly how we treat TP, that is go through and
%   assess each Gaussian-Gaussian pair.
% - FN is calculated using only Int_matrix==1, but Int_matrix has values up to 4. Is it better to
%   just use Int_matrix>0 ?
%
%
%%%%%%%%%%%%%%% Bugs in my (this) code:
% - My Corum_Dataset is 11516x1. Nick's Corum_Dataset is 1x16393. Why are they different sizes?


diary([user.maindir 'logfile.txt'])
disp('ROC_PCPSILAC.m')



%% 0. Initialize
tic
fprintf('\n    0. Initialize')



% Load user settings
maindir = user.maindir;
Experimental_channels = user.silacratios;
desiredPrecision = user.desiredPrecision;
number_of_channels = length(user.silacratios);
%InputFile{1} = user.majorproteingroupsfile;
%InputFile{2} = user.corumfile;
%InputFile{4} = user.mastergaussian;



% Define folders, i.e. define where everything lives.
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir = [maindir 'Data/']; % where data files live
datadir1 = [datadir 'ROC/']; % where data files live
datadir2 = [datadir 'ROC/tmp/']; % where data files live
datadir3 = [datadir 'ROC/CombinedResults/']; % where data files live
figdir1 = [maindir 'Figures/ROC/']; % where figures live
%tmpdir = '/Users/Mercy/Academics/Foster/Jenny_PCPSILAC/PCPSILAC_Analysis/Data/Alignment/';
% Make folders if necessary
if ~exist(codedir, 'dir'); mkdir(codedir); end
if ~exist(funcdir, 'dir'); mkdir(funcdir); end
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist(datadir1, 'dir'); mkdir(datadir1); end
if ~exist(datadir2, 'dir'); mkdir(datadir2); end
if ~exist(datadir3, 'dir'); mkdir(datadir3); end
if ~exist(figdir1, 'dir'); mkdir(figdir1); end

if user.nickflag==1
  %define Raw SILAC ratios data, This is the output from the alignment script
  tmpdir = '/Users/Mercy/Academics/Foster/NickCodeData/4B_ROC_homologue_DB/Combined Cyto new DB/';
  MvsL_filename_Raw_rep1=[tmpdir 'MvsL_alignment/Realignment/Adjusted_MvsL_Raw_for_ROC_analysis_rep1.csv'];
  MvsL_filename_Raw_rep2=[tmpdir 'MvsL_alignment/Realignment/Adjusted_MvsL_Raw_for_ROC_analysis_rep2.csv'];
  MvsL_filename_Raw_rep3=[tmpdir 'MvsL_alignment/Realignment/Adjusted_MvsL_Raw_for_ROC_analysis_rep3.csv'];
  HvsL_filename_Raw_rep1=[tmpdir 'HvsL_alignment/Realignment/Adjusted_HvsL_Raw_for_ROC_analysis_rep1.csv'];
  HvsL_filename_Raw_rep2=[tmpdir 'HvsL_alignment/Realignment/Adjusted_HvsL_Raw_for_ROC_analysis_rep2.csv'];
  HvsL_filename_Raw_rep3=[tmpdir 'HvsL_alignment/Realignment/Adjusted_HvsL_Raw_for_ROC_analysis_rep3.csv'];
  MvsL_filename_gaus_rep1=[tmpdir 'MvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep1.csv'];
  MvsL_filename_gaus_rep2=[tmpdir 'MvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep2.csv'];
  MvsL_filename_gaus_rep3=[tmpdir 'MvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep3.csv'];
  HvsL_filename_gaus_rep1=[tmpdir 'HvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep1.csv'];
  HvsL_filename_gaus_rep2=[tmpdir 'HvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep2.csv'];
  HvsL_filename_gaus_rep3=[tmpdir 'HvsL_alignment/Realignment/Adjusted_Combined_OutputGaus_rep3.csv'];
  ChromatogramIn={MvsL_filename_Raw_rep1,MvsL_filename_Raw_rep2,MvsL_filename_Raw_rep3,...
    HvsL_filename_Raw_rep1,HvsL_filename_Raw_rep2,HvsL_filename_Raw_rep3};
  GaussIn={MvsL_filename_gaus_rep1,MvsL_filename_gaus_rep2,MvsL_filename_gaus_rep3,...
    HvsL_filename_gaus_rep1,HvsL_filename_gaus_rep2,HvsL_filename_gaus_rep3};
  
else
  dd = dir([datadir 'Alignment/Adjusted*Raw_for_ROC_analysis*rep*csv']);
  
  % Define Raw SILAC ratios data. Dynamically find filenames
  if user.skipalignment==1 || isempty(dd)
    % If Alignment was skipped, use raw data + Gauss_Build output
    
    dd = dir([datadir 'GaussBuild/*_Raw_data_maxquant_rep*.csv']);
    ChromatogramIn = cell(size(dd));
    for di = 1:length(dd)
      ChromatogramIn{di} = [datadir 'GaussBuild/' dd(di).name];
    end
    
    GaussIn = cell(size(dd));
    dd = dir([datadir 'GaussBuild/*Combined_OutputGaus*rep*csv']);
    for di = 1:length(dd)
      GaussIn{di} = [datadir 'GaussBuild/' dd(di).name];
    end
  else
    % If Alignment was not skipped, use Alignment output
    
    dd = dir([datadir 'Alignment/Adjusted*Raw_for_ROC_analysis*rep*csv']);
    ChromatogramIn = cell(size(dd));
    for di = 1:length(dd)
      ChromatogramIn{di} = [datadir 'Alignment/' dd(di).name];
    end
    
    GaussIn = cell(size(dd));
    dd = dir([datadir 'Alignment/Adjusted_Combined_OutputGaus*rep*csv']);
    for di = 1:length(dd)
      GaussIn{di} = [datadir 'Alignment/' dd(di).name];
    end
  end
end

% Define which replicates are treatment
number_of_replicates = length(ChromatogramIn) / number_of_channels;
I = find(strcmp(user.silacratios,user.treatmentcondition));
treatment_replicates = [];
if ~isempty(I)
  treatment_replicates = (1:number_of_replicates)+(I-1)*number_of_replicates;
end
untreatment_replicates = [];
I = 1:number_of_replicates*number_of_channels;
if ~isempty(I)
  untreatment_replicates = I(~ismember(I,treatment_replicates));
end

% pre-allocate final results
Recall = nan(number_of_replicates,number_of_channels);
Precision = nan((size(Recall)));
TPR = nan((size(Recall)));
FPR = nan((size(Recall)));
No_Int = nan((size(Recall)));

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)

%%
for replicate_counter = 1:number_of_replicates*number_of_channels
  
  s = ['\n        Replicate ' num2str(replicate_counter)];
  fprintf(s)
  
  %% 1. Read input data
  tic
  s = '\n        1. Read input data';
  fprintf(s)
  
  Protein_IDs = readproteingroupsfile(user);
  
  % Corum binary interactions
  fid=fopen(user.corumpairwisefile, 'rt');    %Input corum data base.
  Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
  fclose(fid);
  No=length(Corum_Import{1})/2;
  Corum_Protein_names = (reshape(Corum_Import{1,1},2,No)');
  Unique_Corum = unique(Corum_Protein_names);
  
  %     %Raw data for analysis of euclidean distance
  %     %import data files summary data
  %     Maxquant = fopen (ChromatogramIn{replicate_counter});
  %     Maxquant_raw = textscan(Maxquant, '%s', 'Delimiter',',');
  %     fclose(Maxquant);
  %
  %     %Reshape data for use in analysis
  %     Raw_Dimension=size(Maxquant_raw{:});
  %     Raw_data_reshaped=reshape(Maxquant_raw{:},68,(Raw_Dimension(1)/68))';
  %
  %     %divide up reshaped data in correct form
  %     Chromatograms_raw = Raw_data_reshaped(2:end,3:end);
  %     Chromatograms_raw = cellfun(@str2num,Chromatograms_raw);
  
  % Raw SILAC ratios
  tmp = importdata(ChromatogramIn{replicate_counter});
  if isfield(tmp.data,'Sheet1')
    Chromatograms_raw = tmp.data.Sheet1;
  else
    Chromatograms_raw = tmp.data;
  end
  % if the first column is Replicate, remove it
  tmp = length(unique(Chromatograms_raw(:,1)));
  if tmp < user.Nreplicate*2+1 % i.e. the column has unique elements in the ballpark of user.Nreplicate
    Chromatograms_raw = Chromatograms_raw(:,2:end);
  end
  
  %   %Gaus data data for analysis of center analysis
  %   Gaus_import_name = fopen (GaussIn{replicate_counter});
  %   Gaus_import = textscan(Gaus_import_name, '%s', 'Delimiter',',');
  %   fclose(Gaus_import_name);
  %
  %   %Reshape data for use in analysis
  %   Gaus_import_Dimension=size(Gaus_import{:});
  %   Gaus_import_reshaped=reshape(Gaus_import{:},7,(Gaus_import_Dimension(1)/7))';
  
  % Gaussian parameters
  tmp = importdata(GaussIn{replicate_counter}, ',');
  if isfield(tmp.data,'Sheet1')
    tmp1 = tmp.data.Sheet1;
    tmp2 = tmp.textdata.Sheet1;
  else
    tmp1 = tmp.data;
    tmp2 = tmp.textdata;
  end
  % Ensure tmp2 is a single column of protein names
  if size(tmp2,2)~=1
    disp('Comparison: Error: violated assumption, textdata is not single column of protein names')
  end
  tmp2 = tmp2(:,1);
  % Ensure protein names are formatted correctly
  for ii = 2:size(tmp2,1)
    tmp2{ii} = strrep(tmp2{ii},',','');
  end
  % Correct bug found for Nick's tissue data.
  % Isoform number (should be textdata) becomes the first column of data.
  if size(tmp1,2)==7
    I = find(~isnan(tmp1(:,7)))';
    if ~isempty(I)
      for jj = I
        tmp2{jj+1} = [tmp2{jj+1} num2str(tmp1(jj,1))];
        tmp1(jj,1:6) = tmp1(jj,2:end);
        tmp1(jj,end) = nan;
      end
    end
  end
  tmp1 = tmp1(:,1:6);
  Gaus_import = [tmp2 num2cell( cat(1,zeros(1,size(tmp1,2)),tmp1) )];
  % Gaus_import(:,1): Protein name
  % Gaus_import(:,2): Height
  % Gaus_import(:,3): Center
  % Gaus_import(:,4): Width
  % Gaus_import(:,5): SSE
  % Gaus_import(:,6): adjrsquare
  % Gaus_import(:,7): Complex size
  
  
  % do a bit of housekeeping
  clear Maxquant_raw Raw_data_reshaped Corum_Import tmp
  
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 2. Pre-process gaussians, chromatograms. Get Dist
  % - In chromatograms, replaces nans with 0.05
  % - Remove gaussians and chromatograms with C<5
  % - Make normalized gaussians, i.e. set height to 1
  tic
  fprintf('        2. Clean gaussians and chromatograms, calculate Dist.')
  
  % Pre-process data to remove Gaussian above fraction five
  % Create array to use for comparsions
  %H_raw = cellfun(@str2num,Gaus_import(2:end,2));
  %C_raw = cellfun(@str2num,Gaus_import(2:end,3));
  %W_raw = cellfun(@str2num,Gaus_import(2:end,4));
  H_raw = cell2mat(Gaus_import(2:end,2));
  C_raw = cell2mat(Gaus_import(2:end,3));
  W_raw = cell2mat(Gaus_import(2:end,4));
  
  % The data is nan-padded. Find where the real data starts and stops.
  nanmax = size(Chromatograms_raw,1);
  tmp = find(sum(isnan(Chromatograms_raw))==size(Chromatograms_raw,1));
  if isempty(tmp)
    tmp = -1;
  end
  frac1 = max([1 tmp(find(tmp<user.Nfraction/2,1,'last'))]); % start of real data
  frac2 = min([size(Chromatograms_raw,2) tmp(find(tmp>user.Nfraction/2,1,'first'))]); % end of real data
  
  % Replace NaN values with 0.05
  Chromatograms = Chromatograms_raw;
  I = find(isnan(Chromatograms));
  Chromatograms(I)= 0.05 * rand(size(I));
  
  % remove gaussian with centers below frac1
  %Ibad = C_raw<frac1;
  Ibad = zeros(size(C_raw));
  H = H_raw(~Ibad);
  C = C_raw(~Ibad);
  W = W_raw(~Ibad);
  Chromatograms = Chromatograms(~Ibad,:);
  Chromatograms_raw = Chromatograms_raw(~Ibad,:);
  
  % Make normalized gaussians
  %   x = 0.25:0.25:100;
  %   gaussian_fit=zeros(length(H),length(x));
  %   for ii = 1:length(H)
  %     gaussian_fit(ii,:) = 1*exp(-(x- C(ii)).^2 /2/(W(ii).^2));
  %   end
  
  % How many Gaussians?
  I = find(~Ibad);
  unique_names = Gaus_import(I+1,1);
  X = repmat(unique_names,1,length(unique_names));
  Ngauss = sum(strcmp(X,X'),1)';
  clear X unique_names
  
  % Co-Apex score = norm(C) / sqrt(length(C))
%   I = find(abs(diff(Ngauss))>0 | Ngauss(2:end)==1);
%   I = [1; I+1];
%   I0 = 0;
%   CoApex = zeros(size(Ngauss));
%   for ii = 1:length(I)
%     I2 = I0+1 : I0 + Ngauss(I(ii));
%     CoApex(I2) = norm(C(I2)) / sqrt(length(I2));
%     I0 = max(I2);
%   end
  
  % Co-Apex score = norm(C) / sqrt(length(C))
  CoApex = zeros(size(Ngauss));
  for ii = 2:size(Gaus_import,1)
    protName = Gaus_import{ii,1};
    I = find(ismember(Gaus_import(:,1),protName)) - 1;
    CoApex(ii-1) = norm(C(I)) / sqrt(length(I));
  end
  
  % Area under the chromatogram
  auc = sum(Chromatograms,2);
  
  % Reduce variables from gaussian-level to protein-level
  %   Ireduce = find(abs(diff(Ngauss))>0 | Ngauss(2:end)==1);
  %   Ireduce = [1; Ireduce+1];
  %   Chromatograms = Chromatograms(Ireduce,:);
  %   C = C(Ireduce,:);
  %   Ngauss = Ngauss(Ireduce);
  %   CoApex = CoApex(Ireduce);
  %   auc = auc(Ireduce,:);
  %   Gaus_import_reshaped(2:end,:) = Gaus_import_reshaped(Ireduce+1,:);
  
  % Calculate distance matrices
  clear Dist
  Dist.Euc = squareform(pdist(Chromatograms,'euclidean'));
  Dist.Center = squareform(pdist(C,'euclidean'));
  Dist.R2 = 1 - corr(Chromatograms').^2; % one minus R squared
  Dist.Ngauss = squareform(pdist(Ngauss));
  Dist.CoApex = squareform(pdist(CoApex));
  Dist.AUC = squareform(pdist(auc));
  %[R,p] = corrcoef(Chromatograms_raw','rows','pairwise');
  %Dist.R2raw = 1 - R.^2;
  %Dist.Rpraw = p;
  
  % GRAVEYARD OF RELEGATED DIST MATRICES
  %   Dist.RawOverlap = nan(size(Chromatograms_raw,1),size(Chromatograms_raw,1));
  %   for ii = 1:size(Chromatograms_raw,1)
  %     for jj = 1:size(Chromatograms_raw,1)
  %       Dist.RawOverlap(ii,jj) = sum(~isnan(Chromatograms_raw(ii,:)) & ~isnan(Chromatograms_raw(jj,:)));
  %     end
  %   end
  %   Dist.R2raw = corr(Chromatograms_raw','rows','pairwise').^2; % one minus R squared
  %   Dist.R2raw2 = Dist.R2raw.*Noverlap;
  %   Dist.Height = squareform(pdist(H,'euclidean'));
  %   Dist.Width = squareform(pdist(W,'euclidean'));
  %   Dist.Gaussian_fits = squareform(pdist(gaussian_fit,'euclidean'));
  %   Dist.dtw = ones(size(Dist.R2))*100;
  %    for ii = 1:size(Chromatograms,1)
  %      if mod(ii,100)==0;disp(num2str(ii));end
  %      for jj = 1:size(Chromatograms,1)
  %        if Dist.R2(ii,jj)<0.4
  %          Dist.dtw(ii,jj) = dtw(Chromatograms(ii,:)',Chromatograms(jj,:)');
  %          %Dist.dtw(ii,jj) = dtw_old(Chromatograms(ii,:),Chromatograms(jj,:));
  %        end
  %      end
  %    end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 3. Make Protein structure
  % This summarizes the proteins in our sample
  tic
  fprintf('        3. Make Protein structure.')
  
  clear Protein
  
  % Combine Majority Protein ID with Proteins Identified in replicate
  Protein.Isoform=Gaus_import(find(~Ibad)+1,1);
  %Protein.Height=cellfun(@str2num,Gaus_import(2:end,2));
  %Protein.Width=cellfun(@str2num,Gaus_import(2:end,3));
  %Protein.Center=cellfun(@str2num,Gaus_import(2:end,4));
  Protein.Height = H;
  Protein.Width = W;
  Protein.Center = C;
  
  % remove protein name with center less then 5
  %Protein.Isoform(Ibad)=[];
  %Protein.Height(Ibad)=[];
  %Protein.Width(Ibad)=[];
  %Protein.Center(Ibad)=[];
  
  Dimensions_Gaus_import = length(Protein.Isoform);
  Dimension_of_Protein_IDs = size(Protein_IDs);
  
  % Create Protein.NoIsoform by removing the isoform tag from Protein.Isoform
  for Gaus_import_counter1 = 1:Dimensions_Gaus_import(1)
    Find_Gaus_isoforms=strfind(Protein.Isoform(Gaus_import_counter1), '-');
    if ~isempty(Find_Gaus_isoforms{1}) == 1
      length_of_Find_Gaus_Protein=length(Protein.Isoform{Gaus_import_counter1});
      number_of_characters_to_remove2=Find_Gaus_isoforms{1}-1;
      Protein.NoIsoform(Gaus_import_counter1,1)={Protein.Isoform{Gaus_import_counter1}(1:number_of_characters_to_remove2)};
    elseif isempty(Find_Gaus_isoforms{1}) == 1
      Protein.NoIsoform(Gaus_import_counter1,1)=Protein.Isoform(Gaus_import_counter1);
    end
  end
  
  % Create Protein_IDS_no_isoform by removing the isoform tag from Protein_IDs
  Protein_IDS_no_isoform = cell(size(Protein_IDs));
  for Protein_ID_counter1 = 2:Dimension_of_Protein_IDs(1)
    for Protein_ID_counter2 = 1:Dimension_of_Protein_IDs(2)-sum(cellfun('isempty',Protein_IDs(Protein_ID_counter1,:)))
      Find_isoforms=strfind(Protein_IDs(Protein_ID_counter1,Protein_ID_counter2), '-');
      if ~isempty(Find_isoforms{1}) == 1
        length_of_Protein_ID=length(Protein_IDs{Protein_ID_counter1,Protein_ID_counter2});
        number_of_characters_to_remove=Find_isoforms{1}-1;
        Protein_IDS_no_isoform(Protein_ID_counter1,Protein_ID_counter2)=cellstr(Protein_IDs{Protein_ID_counter1,Protein_ID_counter2}(1:number_of_characters_to_remove));
      elseif isempty(Find_isoforms{1}) == 1
        Protein_IDS_no_isoform(Protein_ID_counter1,Protein_ID_counter2)=Protein_IDs(Protein_ID_counter1,Protein_ID_counter2);
      end
    end
  end
  
  % Remove non-unique entries from each row of Protein_IDS_no_isoform
  Protein_IDS_no_isoform_no_dup_processed = cell(size(Protein_IDs));
  for Protein_ID_counter3 = 2:Dimension_of_Protein_IDs(1)
    Array_to_check1=Protein_IDS_no_isoform(Protein_ID_counter3,:);
    Array_to_check_empty=cellfun(@isempty,Array_to_check1);
    Array_to_check1(Array_to_check_empty)=[]; % Replace empty cells with NaN
    uniqueCells=unique(Array_to_check1);
    for Protein_ID_counter4 = 1:length(uniqueCells)
      Protein_IDS_no_isoform_no_dup_processed{Protein_ID_counter3,Protein_ID_counter4}=uniqueCells(Protein_ID_counter4);
    end
  end
  
  %Look up Major Protein group Information
  % For each protein that was fitted by a Gaussian
  %   Find the row containing that protein name in Protein_IDs
  %   Store all protein names from that row of Protein_IDs
  Protein.MajorIDs=cell(Dimensions_Gaus_import,Dimension_of_Protein_IDs(2));
  Protein.MajorID_NoIsoforms=cell(Dimensions_Gaus_import,Dimension_of_Protein_IDs(2));
  f = 0;
  for Gaus_import_counter1 = 1:Dimensions_Gaus_import(1)
    
    %disp(Lookup_protein_name)
    Lookup_protein_name = Protein.Isoform(Gaus_import_counter1);
    
    tmp = find(strcmp(Lookup_protein_name,Protein_IDs));
    [Check_counter1, Check_counter2] = ind2sub(size(Protein_IDs),tmp);
    
    %copy names to MajorID from Protein_IDs raw data
    Array_to_check1 = Protein_IDs(Check_counter1,:);
    Array_to_check_empty1 = cellfun(@isempty,Array_to_check1);
    for Check_counter3 = 1:nnz(~Array_to_check_empty1)
      Protein.MajorIDs{Gaus_import_counter1,Check_counter3} = {Protein_IDs{Check_counter1,Check_counter3}};
    end
    
    %copy names to MajorID_NoIsoforms from Protein_IDS_no_isoform_no_dup_processed
    Array_to_check2=Protein_IDS_no_isoform_no_dup_processed(Check_counter1,:);
    Array_to_check_empty2=cellfun(@isempty,Array_to_check2);
    for Check_counter4 = 1:nnz(~Array_to_check_empty2)
      a = Protein_IDS_no_isoform_no_dup_processed{Check_counter1,Check_counter4};
      Protein.MajorID_NoIsoforms{Gaus_import_counter1,Check_counter4}=a{1};
    end
  end
  
  dimension_Protein_MajorID_NoIsoforms=size(Protein.MajorID_NoIsoforms);
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 4. Make TP, FP matrices
  tic
  fprintf('        4. Make TP, FP, FN matrices')
  
  % Make Pos_Corum_proteins, which is all individual proteins that are in corum and our dataset
  try
    Unique_MajorProteinIDs = unique(vertcat(Protein.MajorID_NoIsoforms{:}),'rows');
  catch
    Unique_MajorProteinIDs = unique(vertcat(Protein.MajorID_NoIsoforms(:)),'rows');
  end
  Pos_Corum_proteins = Unique_Corum(ismember(Unique_Corum, Unique_MajorProteinIDs));
  
  % Make Corum_Dataset, which is all interactions that are in corum and our dataset
  cc = 0;
  Corum_Dataset = cell(1000,1);
  for jj = 1:length(Corum_Protein_names); %Write out values as i
    Prot1 = Corum_Protein_names(jj,1);
    Prot2 = Corum_Protein_names(jj,2);
    
    % find these protein names in our sample
    tmp1 = sum(strcmp(Prot1,Pos_Corum_proteins));
    tmp2 = sum(strcmp(Prot2,Pos_Corum_proteins));
    
    if tmp1>0 && tmp2>0
      cc = cc+1;
      Corum_Dataset{cc} = [Prot1{1},'-',Prot2{1}];
    end
  end
  
  % Expand Corum_Dataset to deal with protein groups
  % This all has to do with Major_protein_groups.xlsx.
  % In effect, some proteins have multiple names.
  % e.g. Interaction A-B is not in corum, but A-C is, and B and C are in the same protein group, so
  %   treat A-B as a true positive interaction.
  
  TP_Matrix= zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
  %Corum_Dataset_expanded = cell(size(Corum_Dataset));
  ii = 0;
  for cc=1:size(Corum_Dataset,1)
    I_hyphen = strfind(Corum_Dataset{cc},'-');
    Prot1 = Corum_Dataset{cc}(1:I_hyphen-1);
    Prot2 = Corum_Dataset{cc}(I_hyphen+1:end);
    
    % find where Prot1 and Prot2 are in Protein.MajorID_NoIsoforms
    ind1 = find(strcmp(Prot1,Protein.MajorID_NoIsoforms));
    [i1,j1] = ind2sub(size(Protein.MajorID_NoIsoforms),ind1);
    ind2 = find(strcmp(Prot2,Protein.MajorID_NoIsoforms));
    [i2,j2] = ind2sub(size(Protein.MajorID_NoIsoforms),ind2);
    
    TP_Matrix(i1,i2) = 1;
  end
  %I = unique(cell2mat(Corum_Dataset_expanded),'rows');
  %Corum_Dataset_expanded = unique(Corum_Dataset_expanded{ii});
  
  
  % Make other matrices
  % - Possible interaction matrix
  % - Self interaction matrix
  
  % check which individual proteins are in corum
  % Important! Also check whether any member of that protein's group is in corum
  Ngroup = size(Protein.MajorID_NoIsoforms,2);
  inCorum = zeros(Ngroup,Dimensions_Gaus_import);
  for ni = 1:Ngroup
    tmp = Protein.MajorID_NoIsoforms(:,ni);
    % ismember can't handle empty cells, so as a workaround fill with a string that won't be in Corum
    isemptyCellArray = cellfun('isempty',tmp);
    tmp(isemptyCellArray) = {'--fake'};
    inCorum(ni,:) = ismember(tmp,Unique_Corum);
  end
  inCorum = sum(inCorum);
  
  % Calculate possible interaction matrix, asks are both proteins in corum
  Int_matrix = inCorum'*inCorum;
  
  % Create self interaction matrix
  Self_Int_matrix = zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
  for i=1:Dimensions_Gaus_import
    for ii=1:Dimensions_Gaus_import
      if isequal(Protein.NoIsoform{i},Protein.NoIsoform{ii})
        Self_Int_matrix(i,ii)=1;
      end
    end
  end
  inverse_self = ~Self_Int_matrix; %creation of inverse of self interactions
  clear Self_Int_matrix
  
  % Possible interactions, i.e. both in corum and not a self interaction
  possibleInts = inverse_self & Int_matrix;
  
  % FP matrix
  Inverse_TP_matrix = ~TP_Matrix;
  FP_Matrix = (Int_matrix & Inverse_TP_matrix);
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  %% 5. Calculate score, the likelihood that a protein pair is interacting
  tic
  fprintf('        5. Score each interaction via Naive Bayes classifier')
  % score is proportional the likelihood that a protein pair is interacting
  
  scoreMatrix = scorenb(Dist,possibleInts,TP_Matrix);
  scoreMatrix = nanmedian(scoreMatrix,2);
  
  sf = [datadir2 'score_rep' num2str(replicate_counter) '.mat'];
  save(sf,'scoreMatrix','TP_Matrix','possibleInts','Protein','inverse_self','Chromatograms','Dist')
  clear scoreMatrix TP_Matrix possibleInts Protein inverse_self Chromatograms Dist
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
end
%%%%% Replicate counter ends




% do a bit of housekeeping
clear Int_matrix


%% 6. Find the score thresholds

tic
fprintf('    6. Find the score cutoff across all replicates')

% 6a. Concatenate scores across replicates.
% The following code converts Protein1 and Protein2 to unique (hopefully!) numeric values.
% This greatly reduces the memory used, as otherwise allScores is ~2.5GB.
% NB: It's possible that two protein names will be mapped to the same numeric value.
allScores = zeros(10^7,5);
xcutoff_rep = nan(number_of_replicates*number_of_channels, length(desiredPrecision));
calcprec_rep = nan(number_of_replicates*number_of_channels, length(desiredPrecision));
calcrec_rep = nan(number_of_replicates*number_of_channels, length(desiredPrecision));
ff = rand(1,15);
kk = 0;
for rr = 1:(number_of_replicates*number_of_channels)
  sf = [datadir2 'score_rep' num2str(rr) '.mat'];
  load(sf)
  a = find(triu(possibleInts));
  [prot1i, prot2i] = ind2sub(size(possibleInts),a);
  I = kk+1 : kk+length(a);
  allScores(I,1) = TP_Matrix(a);           % Class label
  allScores(I,4) = rr;         % Replicate
  allScores(I,5) = median(scoreMatrix(a),2);  % Median score (equivalent to ensemble voting)
  
  % For space reasons, let's try converting Protein1, Protein2 into numeric
  kk0 = kk;
  for ii = 1:length(prot1i)
    kk = kk+1;
    tmp = double(Protein.Isoform{prot1i(ii)});
    allScores(kk,2) = sum(tmp .* ff(1:length(tmp)));
    tmp = double(Protein.Isoform{prot2i(ii)});
    allScores(kk,3) = sum(tmp .* ff(1:length(tmp)));
  end
  
  % Check the name-to-numeric conversion
  % Did we lose any unique names in the conversion?
  N1a = length(unique(Protein.Isoform(prot1i)));
  N2a = length(unique(Protein.Isoform(prot2i)));
  N1b = length(unique(allScores(kk0+1:kk,2)));
  N2b = length(unique(allScores(kk0+1:kk,3)));
  if N1a~=N1b || N2a~=N2b
    disp('Two protein names stored as the same number!')
  end
  
  % Make non-redundant replicate-specific score/class
  scorestmp = allScores(kk0+1:kk,:);
  [~,I1,I2] = unique(scorestmp(:,2:3),'rows');
  allScores2 = nan(length(I1),100);
  countV = zeros(length(I1),1);
  for ii = 1:length(I2)
    countV(I2(ii)) = countV(I2(ii))+1; % how many times have we seen this interaction?
    allScores2(I2(ii),1) = scorestmp(ii,1);
    allScores2(I2(ii),countV(I2(ii))+1) = scorestmp(ii,5);
  end
  class_rep = allScores2(:,1);
  score_rep = max(allScores2(:,2:end),[],2);      % i) max score
  
  % Find replicate-specific threshold
  for di = 1:length(desiredPrecision)
    [xcutoff_rep(rr,di), tmp] = calcPPIthreshold(score_rep,class_rep,desiredPrecision(di));
    calcprec_rep(rr,di) = tmp.calcprec;
    calcrec_rep(rr,di) = tmp.calcrec;
  end
  
  clear scoreMatrix inverse_self Protein possibleInts TP_Matrix scorestmp Dist
end
allScores = allScores(1:kk,:);

% Remove redundant interactions
% Options: i) max score, ii) mean score, iii) median score
[~,I1,I2] = unique(allScores(:,2:3),'rows');
allScores2 = nan(length(I1),100);
countV = zeros(length(I1),1);
for ii = 1:length(I2)
  countV(I2(ii)) = countV(I2(ii))+1; % how many times have we seen this interaction?
  allScores2(I2(ii),1) = allScores(ii,1);
  allScores2(I2(ii),countV(I2(ii))+1) = allScores(ii,5);
end
class = allScores2(:,1);
score = max(allScores2(:,2:end),[],2);      % i) max score
%score = nanmean(allScores2(:,2:end),2);    % ii) mean score
%score = nanmedian(allScores2(:,2:end),2);  % iii) median score
clear allScores2 allScores TP_Matrix I1 I2


% 6b. Calculate global threshold
xcutoff = nan(size(desiredPrecision));
calcprec = zeros(size(xcutoff));
calcrec = zeros(size(xcutoff));
scoreRange = [];                       % Save
precRange = [];                        % these
recRange = [];                         % for
tprRange = [];                         % plotting
fprRange = [];                         % later.
Ninteract = [];                        % ...
for di = 1:length(desiredPrecision)
  [xcutoff(di), tmp] = calcPPIthreshold(score,class,desiredPrecision(di));
  calcprec(di) = tmp.calcprec;
  calcrec(di) = tmp.calcrec;
  
  % save plotting variables
  I = ~isnan(tmp.scoreRange);
  scoreRange = [scoreRange; tmp.scoreRange(I)'];
  precRange = [precRange; tmp.precRange(I)'];
  recRange = [recRange; tmp.recRange(I)'];
  tprRange = [tprRange; tmp.tprRange(I)'];
  fprRange = [fprRange; tmp.fprRange(I)'];
  Ninteract = [Ninteract; tmp.Ninteractions(I)'];
end
[scoreRange,I] = sort(scoreRange);
precRange = precRange(I);
recRange = recRange(I);
tprRange = tprRange(I);
fprRange = fprRange(I);
Ninteract = Ninteract(I);


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)




%% 7. Find and concatenate interactions at desired precision

firstFlag = 1;
for pri = 1:length(desiredPrecision)
  tic
  s = ['        7. Find and concatenate interactions at precision = ' num2str(desiredPrecision(pri)) '\n'];
  fprintf(s)
  
  % Set these here, since we want to concatenate across replicates
  fn_count = 0; % false negative interaction counter, used in 7b
  interaction_count = 0; % interaction counter, used in 7c
  
  binary_interaction_list = cell(2*10^6,13);
  Neg_binary_interaction_list = zeros(2*10^6,3);
  
  for replicate_counter = 1:number_of_replicates*number_of_channels
    
    % 7b. Calculate recall and precision, negative interactions
    %Determine the false "binary" interaction list
    tic
    fprintf(['            Replicate ' num2str(replicate_counter)])
    
    % Re-load scoreMatrix, TP_Matrix, possibleInts, and Protein
    sf = [datadir2 'score_rep' num2str(replicate_counter) '.mat'];
    load(sf)
    
    % determine final interaction matrices
    pos = scoreMatrix>xcutoff(pri); % positive hits
    pos = reshape(pos,size(inverse_self,1),size(inverse_self,1));
    Final_interactions = pos & inverse_self;
    
    % Calculate TP, FP, TN, FN
    TP = sum(Final_interactions(:) & TP_Matrix(:) & possibleInts(:));
    FP = sum(Final_interactions(:) & ~TP_Matrix(:) & possibleInts(:));
    TN = sum(~Final_interactions(:) & ~TP_Matrix(:) & possibleInts(:));
    FN = sum(~Final_interactions(:) & TP_Matrix(:) & possibleInts(:));
    
    %Generate negative matrix for final interaction analysis
    %Negative2_matrix = ~Final_interactions & possibleInts;
    
    % final results
    %Recall(replicate_counter,pri) = TP/(TP+FN);
    %Precision(replicate_counter,pri) = TP/(TP+FP);
    %TPR(replicate_counter,pri) = TP/(TP+FN);
    %FPR(replicate_counter,pri) = FP/(FP+TN);
    %No_Int(replicate_counter,pri) = length(find(triu(Final_interactions==1)));
    
    % Find "binary" NON interactions, both proteins in corum
    scoreMatrix = nanmedian(scoreMatrix,2);
    scoreMatrix = reshape(scoreMatrix,size(Final_interactions,1),size(Final_interactions,1));
    Protein_elution = strcat(num2str(Protein.Center),'*',Protein.Isoform);
    U1 = triu(~Final_interactions & possibleInts); 
    %U1 = triu(~Final_interactions & TP_Matrix); % Take half of matrix therefor unique interactions
    [ia,ib] = find(U1 == 1); % protein1, protein2 indices in false negatives
    for ri=1:length(ia)
      ri1 = ia(ri);
      ri2 = ib(ri);
      Int1=Protein_elution(ri1);
      Int2=Protein_elution(ri2);
      if Protein.Center(ri1)==Protein.Center(ri2) && ri1==ri2 % make sure it's not a self interaction
      else
        fn_count = fn_count+1;
        tmp = double(Protein.Isoform{ri1});
        tmpname1 = sum(tmp .* ff(1:length(tmp)));
        tmp = double(Protein.Isoform{ri2});
        tmpname2 = sum(tmp .* ff(1:length(tmp)));
        Neg_binary_interaction_list(fn_count,1) = tmpname1; % Protein_interactions1
        Neg_binary_interaction_list(fn_count,2) = tmpname2; % Protein_interactions2
        Neg_binary_interaction_list(fn_count,3) = TP_Matrix(ri1,ri2); % Interaction_in_Corum
      end
      % sanity check
      if ~strcmp(Int1, Int2) && Protein.Center(ri1)==Protein.Center(ri2) && ri1==ri2
        disp('uh oh!')
      end
    end
    
    
    % Find "binary" interactions
    Final_interactions = triu(Final_interactions);
    [I1,I2] = find(Final_interactions>0);
    for ii = 1:length(I1)
      vii = I1(ii);
      viii = I2(ii);
      if vii == viii
        continue;
      end
      Int1=Protein_elution{vii};
      Int2=Protein_elution{viii};
      interaction_count=1+interaction_count;
      binary_interaction_list{interaction_count,1} = strcat(Int1,'-',Int2); % Interactions
      %binary_interaction_list{interaction_count,2} = norm(Chromatograms(vii,:)-Chromatograms(viii,:)); % Delta_EucDist
      %binary_interaction_list{interaction_count,3} = abs(Protein.Center(vii)-Protein.Center(viii)); % Delta_Center
      %binary_interaction_list{interaction_count,4} = abs(Protein.Height(vii)-Protein.Height(viii)); % Delta_Height
      %binary_interaction_list{interaction_count,5} = abs(Protein.Width(vii)-Protein.Width(viii)); % Delta_Width
      binary_interaction_list{interaction_count,2} = Dist.Euc(vii,viii); % Euc
      binary_interaction_list{interaction_count,3} = Dist.Center(vii,viii); % Co-apex
      binary_interaction_list{interaction_count,4} = Dist.R2(vii,viii); % R^2
      binary_interaction_list{interaction_count,5} = 0; % empty for now
      binary_interaction_list{interaction_count,6} = possibleInts(vii,viii); % Both_proteins_Corum
      binary_interaction_list{interaction_count,7} = TP_Matrix(vii,viii); % Interaction_in_Corum
      binary_interaction_list{interaction_count,8} = Protein.Isoform{vii}; % Protein_interactions1
      binary_interaction_list{interaction_count,9} = Protein.Isoform{viii}; % Protein_interactions2
      binary_interaction_list{interaction_count,10} = Protein.Center(vii); % Protein_interactions_center1
      binary_interaction_list{interaction_count,11} = Protein.Center(viii); % Protein_interactions_center2
      binary_interaction_list{interaction_count,12} = replicate_counter; % Protein_interactions_center2
      binary_interaction_list{interaction_count,13} = scoreMatrix(vii,viii); % Interaction score
    end
    
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
    
  end
  
  binary_interaction_list = binary_interaction_list(1:interaction_count,:);
  Neg_binary_interaction_list = Neg_binary_interaction_list(1:fn_count,:);
  
  
  
  tic
  fprintf('        8. Create a list of all the unique interactions')
  
  interaction_pairs = cell(size(binary_interaction_list,1),1);
  for ri=1:size(binary_interaction_list,1)
    interaction_pairs{ri} = [binary_interaction_list{ri,8} '_' binary_interaction_list{ri,9}];
  end
  unique_interaction = unique(interaction_pairs);
  length_unique_inter = length(unique_interaction);
  
  %Create array to store values
  interaction_final.unique_interactions=cell(length_unique_inter,0);
  interaction_final.replicate_numbers=cell(length_unique_inter,0);
  interaction_final.proteinA=cell(length_unique_inter,0);
  interaction_final.proteinB=cell(length_unique_inter,0);
  interaction_final.CenterA=cell(length_unique_inter,0);
  interaction_final.CenterB=cell(length_unique_inter,0);
  %interaction_final.DeltaHeight=cell(0,0);
  interaction_final.DeltaCenter=cell(length_unique_inter,0);
  interaction_final.R2=cell(length_unique_inter,0);
  interaction_final.DeltaEucDist=cell(length_unique_inter,0);
  interaction_final.proteinInCorum=cell(length_unique_inter,0);
  interaction_final.interactionInCorum=cell(length_unique_inter,0);
  interaction_final.score=nan(length_unique_inter,60);
  Unique_interaction_counter=0;
  for interaction_counter2A = 1:length_unique_inter
    %Find location of unique protein interaction
    location_interaction_pairs = find(strcmp(unique_interaction(interaction_counter2A), interaction_pairs));
    
    %find which replicate the interactions were seen in
    diffC = abs(cellfun(@minus,binary_interaction_list(location_interaction_pairs,10),...
      binary_interaction_list(location_interaction_pairs,11)));
    position_within_two = find(diffC<2);
    
    %Check if multiple centers are withing two fractions of each other
    if ~nnz(position_within_two)==0
      %Add one to counters and remove location already checked
      %row_counter=row_counter+1;
      Unique_interaction_counter=Unique_interaction_counter+1;
      
      interaction_final.unique_interactions(Unique_interaction_counter,1) = unique_interaction(interaction_counter2A);
      interaction_final.proteinInCorum(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(position_within_two(1)),6);
      interaction_final.interactionInCorum(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(position_within_two(1)),7);
      
      for row = 1:length(position_within_two)
        %check size of replicate_numbers
        interaction_final.replicate_numbers{Unique_interaction_counter,row} = ...
          binary_interaction_list{location_interaction_pairs(position_within_two(row)),12};
        
        % copy the values for the first unique interaction to the list
        interaction_final.proteinA(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),8);
        interaction_final.proteinB(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),9);
        interaction_final.CenterA(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),10);
        interaction_final.CenterB(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),11);
        %interaction_final.DeltaHeight(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),4);
        interaction_final.DeltaCenter(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),3);
        interaction_final.R2(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),4);
        interaction_final.DeltaEucDist(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),2);
        interaction_final.score(Unique_interaction_counter,row) = binary_interaction_list{location_interaction_pairs(position_within_two(row)),13};
        
      end
      %location_interaction_pairs(position_within_two)=[];
      
      %If center not within two fractions
    elseif nnz(position_within_two)==0
      %Add one to counters and remove location already checked
      %location_interaction_pairs(1)=[];
      Unique_interaction_counter=Unique_interaction_counter+1;
      
      interaction_final.unique_interactions(Unique_interaction_counter,1)= unique_interaction(interaction_counter2A);
      %check size of replicate_numbers
      interaction_final.replicate_numbers{Unique_interaction_counter,1}=...
        binary_interaction_list{location_interaction_pairs(1),12};
      
      %copy the values for the first unique interaction to the list
      interaction_final.proteinA(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),8);
      interaction_final.proteinB(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),9);
      interaction_final.CenterA(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),10);
      interaction_final.CenterB(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),11);
      %interaction_final.DeltaHeight(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),4);
      interaction_final.DeltaCenter(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),3);
      interaction_final.R2(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),4);
      interaction_final.DeltaEucDist(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),2);
      interaction_final.proteinInCorum(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),6);
      interaction_final.interactionInCorum(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),7);
      interaction_final.score(Unique_interaction_counter,1) = binary_interaction_list{location_interaction_pairs(1),13};
    end
  end
  
  % Calculate the average score for each pairwise interaction
  interaction_final.score = nanmean(interaction_final.score,2);
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  % Calculate final TP, FP, FN and TN
  tic
  fprintf('        9. Calculate final TP, FP, FN and TN')
  % True interactions Observation across replicate
  % Count how many times a true positive interactions was observed
  Total_unique_interactions=length(interaction_final.unique_interactions);
  number_observation_pre_interaction=zeros(Total_unique_interactions,1);
  number_unique_interaction=zeros(Total_unique_interactions,1);
  %Record the number of times a interaction was detected in irrespective of replicate number
  interactions_prep_replicate=zeros(Total_unique_interactions,replicate_counter);
  for counter=1:Total_unique_interactions
    
    %convert cell array to double
    a = interaction_final.replicate_numbers(counter,:);
    ix = ~cellfun('isempty',a);
    testA = [a{ix}];
    %Save values
    number_observation_pre_interaction(counter) = nnz(testA);
    
    %Determine if multiple gaussians within the same replicate were indentified
    unique_testA = unique(testA);
    number_of_unique = length(unique_testA);
    number_unique_interaction(counter) = number_of_unique;
    
    %Which replicate is the interaction in
    for replicate_counter = 1:(number_of_replicates*number_of_channels)
      interactions_prep_replicate(counter,replicate_counter) = nnz(find(testA(:)==replicate_counter));
    end
  end
  
  %Copy number of unique_interaction across replicate to interaction array
  %interaction_final.times_observed_across_reps = number_unique_interaction;
  
  %total number of interactions of unique replicates
  %total_observation = sum(number_observation_pre_interaction);
  %total number of interactions of unique replicates
  %total_unique_observation = sum(number_unique_interaction);
  
  
  % create unique identifier for false interactions
  % ##### THIS CODE IS AN APPROXIMATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  % ##### If protein A and B are not interacting in corum, this sets ALL THE GAUSSIAN COMBINATIONS
  %       between A and B to not interacting. See Dec 3 notes for more info.
  %       The proper fix is to individually treat all Gaussian combinations exactly like we do for
  %       the positive interactions.
  %interaction_neg_pairs = Neg_binary_interaction_list(:,1) + Neg_binary_interaction_list(:,2);
  %length_unique_neg_inter = length(unique(interaction_neg_pairs));
  
  % Determine the number of unique false negatives
  %length_unique_FN_incorum = length(unique(interaction_neg_pairs( Neg_binary_interaction_list(:,3)==1 )));

  % Caliculate TP, FP, FN and TN
  %Caliculate values from true interaction data set
  True_proteinInCorum = cell2mat(interaction_final.proteinInCorum);
  True_interactionInCorum = cell2mat(interaction_final.interactionInCorum);
%  
%   %Caliculate values from true interaction data set
%   False_proteinInCorum = length_unique_neg_inter;
%   False_interactionInCorum = length_unique_FN_incorum;
  
  %caliculate TP, FP, FN and TN #PARSE_HERE
  final_TN = sum(Neg_binary_interaction_list(:,3)==0);
  final_FP = sum(True_proteinInCorum & ~True_interactionInCorum);
  final_FN = sum(Neg_binary_interaction_list(:,3)==1);
  final_TP = sum(True_proteinInCorum & True_interactionInCorum);
  
  % final results
  final_Recall = final_TP/(final_TP+final_FN);
  final_Precision = final_TP/(final_TP+final_FP);
  final_TPR = final_TP/(final_TP+final_FN);
  final_FPR =final_FP/(final_FP+final_TN);
  
  clearvars Neg_binary_interaction_list
  
  
  %Format data to be written out
  %Crete list of centerA, note ensure strjoin function is avalible
  interaction_final.centerA_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    centerA_2bformated1=interaction_final.CenterA(format_loop,:);
    centerA_2bformated1(cellfun('isempty',centerA_2bformated1)) = [];
    length_centerA=length(centerA_2bformated1);
    
    %Test if the array is longer then one entry
    if length_centerA<2
      interaction_final.centerA_formated(format_loop)=centerA_2bformated1(1);
    elseif length_centerA>=2
      for jj= 1:length(centerA_2bformated1)
        centerA_2bformated1{jj} = num2str(centerA_2bformated1{jj});
      end
      interaction_final.centerA_formated(format_loop)= {strjoin(centerA_2bformated1,' ; ')};
    end
  end
  
  %Crete list of centerB, note ensure strjoin function is avalible
  interaction_final.centerB_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    centerB_2bformated1=interaction_final.CenterB(format_loop,:);
    centerB_2bformated1(cellfun('isempty',centerB_2bformated1)) = [];
    length_centerB=length(centerB_2bformated1);
    
    %Test if the array is longer then one entry
    if length_centerB<2
      interaction_final.centerB_formated(format_loop)=centerB_2bformated1(1);
    elseif length_centerB>=2
      for jj= 1:length(centerB_2bformated1)
        centerB_2bformated1{jj} = num2str(centerB_2bformated1{jj});
      end
      interaction_final.centerB_formated(format_loop)= {strjoin(centerB_2bformated1,' ; ')};
    end
  end
  
  %Crete list of DeltaEuc, note ensure strjoin function is avalible
  interaction_final.DeltaEuc_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    DeltaEuc_2bformated1=interaction_final.DeltaEucDist(format_loop,:);
    DeltaEuc_2bformated1(cellfun('isempty',DeltaEuc_2bformated1)) = [];
    length_DeltaEuc=length(DeltaEuc_2bformated1);
    %Test if the array is longer then one entry
    if length_DeltaEuc<2
      interaction_final.DeltaEuc_formated(format_loop)=DeltaEuc_2bformated1(1);
    elseif length_DeltaEuc>=2
      for jj= 1:length(DeltaEuc_2bformated1)
        DeltaEuc_2bformated1{jj} = num2str(DeltaEuc_2bformated1{jj});
      end
      interaction_final.DeltaEuc_formated(format_loop)= {strjoin(DeltaEuc_2bformated1,' ; ')};
    end
  end
  
  %Crete list of DeltaCenter, note ensure strjoin function is avalible
  interaction_final.DeltaCenter_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    DeltaCenter_2bformated1=interaction_final.DeltaCenter(format_loop,:);
    DeltaCenter_2bformated1(cellfun('isempty',DeltaCenter_2bformated1)) = [];
    length_DeltaCenter=length(DeltaCenter_2bformated1);
    
    %Test if the array is longer then one entry
    if length_DeltaCenter<2
      interaction_final.DeltaCenter_formated(format_loop)=DeltaCenter_2bformated1(1);
    elseif length_DeltaCenter>=2
      for jj= 1:length(DeltaCenter_2bformated1)
        DeltaCenter_2bformated1{jj} = num2str(DeltaCenter_2bformated1{jj});
      end
      interaction_final.DeltaCenter_formated(format_loop)= {strjoin(DeltaCenter_2bformated1,' ; ')};
    end
  end
  
  %Crete list of R^2, note ensure strjoin function is avalible
  interaction_final.R2_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    R2_2bformated1=interaction_final.R2(format_loop,:);
    R2_2bformated1(cellfun('isempty',R2_2bformated1)) = [];
    length_R2=length(R2_2bformated1);
    %Test if the array is longer then one entry
    if length_R2<2
      interaction_final.R2_formated(format_loop)=R2_2bformated1(1);
    elseif length_R2>=2
      for jj= 1:length(R2_2bformated1)
        R2_2bformated1{jj} = num2str(R2_2bformated1{jj});
      end
      interaction_final.R2_formated(format_loop)= {strjoin(R2_2bformated1,' ; ')};
    end
  end
  
  %Crete list of replicates, note ensure strjoin function is avalible
  interaction_final.Replicates_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    Replicates_2bformated1=interaction_final.replicate_numbers(format_loop,:);
    Replicates_2bformated1(cellfun('isempty',Replicates_2bformated1)) = [];
    length_Replicates=length(Replicates_2bformated1);
    %Test if the array is longer then one entry
    if length_Replicates<2
      interaction_final.Replicates_formated{format_loop}=mat2str(Replicates_2bformated1{1});
    elseif length_Replicates>=2
      for jj= 1:length(Replicates_2bformated1)
        Replicates_2bformated1{jj} = num2str(Replicates_2bformated1{jj});
      end
      interaction_final.Replicates_formated{format_loop}= strjoin(Replicates_2bformated1,' ; ');
    end
  end
  
  %Crete list of scores, note ensure strjoin function is avalible
  interaction_final.scores_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    tmp = interaction_final.score(format_loop);
    tmp(isnan(tmp)) = [];
    length_Replicates=length(tmp);
    %Test if the array is longer then one entry
    if length_Replicates<2
      interaction_final.scores_formated{format_loop}=mat2str(tmp(1));
    elseif length_Replicates>=2
      Replicates_2bformated1 = cell(length_Replicates,1);
      for jj= 1:length(tmp)
        Replicates_2bformated1{jj} = num2str(tmp(jj));
      end
      interaction_final.scores_formated{format_loop} = strjoin(Replicates_2bformated1,' ; ');
    end
  end
  
  % Create scoreRank, scorePrctile
  [~, ~, rankDescend] = unique(interaction_final.score);
  interaction_final.scoreRank = (1-rankDescend) + max(rankDescend);
  interaction_final.scorePrctile = interaction_final.scoreRank/length(interaction_final.scoreRank);
  
  
  % Create precision-dropout
  interaction_final.precisionDropout = nan(Total_unique_interactions,1);
  tmpProtinCor = cell2mat(interaction_final.proteinInCorum);
  tmpIntinCor = cell2mat(interaction_final.interactionInCorum);
  for ii = 1:Total_unique_interactions
    cutoff = interaction_final.score(ii);
    I = interaction_final.score >= cutoff;
    protInCor = tmpProtinCor(I);
    intInCor = tmpIntinCor(I);
    TP = sum(protInCor==1 & intInCor==1);
    FP = sum(protInCor==1 & intInCor==0);
    interaction_final.precisionDropout(ii) = TP / (TP+FP);
  end
  
  %#PARSE_HERE
  %create figure of precision replicates over replicates
  %create array to fill with number of interaction with both proteins in
  %corum and the number of interactions for these proteins in corum
  Precision_array=zeros((number_of_replicates*number_of_channels),3);
  Precision_redundant=zeros((number_of_replicates*number_of_channels),3);
  
  %Determine the precision for interactions observed across replicates
  % Precision_redundant(:,1) = protein in corum
  % Precision_redundant(:,2) = interaction in corum
  for precision_counter = 1:(number_of_replicates*number_of_channels)
    %Find position of interaction observed within the required number of
    %replicates
    position_to_sum=find(interaction_final.times_observed_across_reps>=precision_counter);
    
    %Count the number values for non redundant list
    Precision_array(precision_counter,1) = sum(cell2mat(interaction_final.proteinInCorum(position_to_sum)));
    Precision_array(precision_counter,2) = sum(cell2mat(interaction_final.interactionInCorum(position_to_sum)));
    Precision_array(precision_counter,3) = length(position_to_sum);
    %Count the number values for redundant list
    Precision_redundant(precision_counter,1) = (sum(cell2mat(interaction_final.proteinInCorum(position_to_sum)))*precision_counter);
    Precision_redundant(precision_counter,2) = (sum(cell2mat(interaction_final.interactionInCorum(position_to_sum)))*precision_counter);
    Precision_redundant(precision_counter,3) = (length(position_to_sum)*precision_counter);
  end
  
  %Generation precision as a precentage for non redundant interactions
  Precent_precision_R=zeros((number_of_replicates*number_of_channels),1);
  for precision_counter=1:(number_of_replicates*number_of_channels)
    Precent_precision_R(precision_counter)=Precision_array(precision_counter,2)/Precision_array(precision_counter,1);
  end
  
  %Generate precision as a precentage for redundant interactions
  Precent_precision_NR=sum(Precision_redundant(:,2))/sum(Precision_redundant(:,1));
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  if ~isempty(treatment_replicates)
    tic
    fprintf('        10. Interactions observed only under treatment')
    Observed_treatment=zeros(Total_unique_interactions,1);
    Observed_untreatment=zeros(Total_unique_interactions,1);
    treatment_specific_interaction=zeros(Total_unique_interactions,1);
    untreatment_specific_interaction=zeros(Total_unique_interactions,1);
    
    %Determine If there are condition specific interactions detected
    for treatment_counter=1:Total_unique_interactions
      
      Observed_treatment(treatment_counter,1)=nnz(interactions_prep_replicate(treatment_counter,treatment_replicates));
      Observed_untreatment(treatment_counter,1)=nnz(interactions_prep_replicate(treatment_counter,untreatment_replicates));
      
      %Determine if there are treatment specific interactions detected
      if Observed_treatment(treatment_counter,1)>=2 & Observed_untreatment(treatment_counter,1)==0
        treatment_specific_interaction(treatment_counter,1)=1;
      end
      
      %Determine if there are untreatment specific interactions detected
      if Observed_untreatment(treatment_counter,1)>=2 & Observed_treatment(treatment_counter,1)==0
        untreatment_specific_interaction(treatment_counter,1)=1;
      end
      
    end
    
    %Copy treatment specific interaction to list to write out
    treatment_specific.unique_interactions=cell(0,0);
    treatment_specific.proteinA=cell(0,0);
    treatment_specific.proteinB=cell(0,0);
    treatment_specific.centerA_formated=cell(0,0);
    treatment_specific.centerB_formated=cell(0,0);
    treatment_specific.Replicates_formated=cell(0,0);
    %treatment_specific.DeltaHeight_formated=cell(0,0);
    treatment_specific.DeltaCenter_formated=cell(0,0);
    treatment_specific.R2_formated=cell(0,0);
    treatment_specific.DeltaEuc_formated=cell(0,0);
    treatment_specific.proteinInCorum=cell(0,0);
    treatment_specific.interactionInCorum=cell(0,0);
    
    %Create counter
    treatment_row_counter=1;
    for treatment_counter=1:Total_unique_interactions
      if  treatment_specific_interaction(treatment_counter,1)==1;
        treatment_specific.unique_interactions(treatment_row_counter,1)=interaction_final.unique_interactions(treatment_counter,1);
        treatment_specific.proteinA(treatment_row_counter,1)=interaction_final.proteinA(treatment_counter,1);
        treatment_specific.proteinB(treatment_row_counter,1)=interaction_final.proteinB(treatment_counter,1);
        treatment_specific.centerA_formated(treatment_row_counter,1)=interaction_final.centerA_formated(treatment_counter,1);
        treatment_specific.centerB_formated(treatment_row_counter,1)=interaction_final.centerB_formated(treatment_counter,1);
        treatment_specific.Replicates_formated(treatment_row_counter,1)=interaction_final.Replicates_formated(treatment_counter,1);
        %treatment_specific.DeltaHeight_formated(treatment_row_counter,1)=interaction_final.DeltaHeight_formated(treatment_counter,1);
        treatment_specific.DeltaCenter_formated(treatment_row_counter,1)=interaction_final.DeltaCenter_formated(treatment_counter,1);
        treatment_specific.R2_formated(treatment_row_counter,1)=interaction_final.Deltawidth_formated(treatment_counter,1);
        treatment_specific.DeltaEuc_formated(treatment_row_counter,1)=interaction_final.DeltaEuc_formated(treatment_counter,1);
        treatment_specific.proteinInCorum(treatment_row_counter,1)=interaction_final.proteinInCorum(treatment_counter,1);
        treatment_specific.interactionInCorum(treatment_row_counter,1)=interaction_final.interactionInCorum(treatment_counter,1);
        treatment_row_counter=1+treatment_row_counter;
      end
    end
    
    %Determine the number of interactions
    treament_specific_interactionsnum=length(treatment_specific.unique_interactions);
    
    
    %Copy treatment specific interaction to list to write out
    untreatment_specific.unique_interactions=cell(0,0);
    untreatment_specific.proteinA=cell(0,0);
    untreatment_specific.proteinB=cell(0,0);
    untreatment_specific.centerA_formated=cell(0,0);
    untreatment_specific.centerB_formated=cell(0,0);
    untreatment_specific.Replicates_formated=cell(0,0);
    %untreatment_specific.DeltaHeight_formated=cell(0,0);
    untreatment_specific.DeltaCenter_formated=cell(0,0);
    untreatment_specific.R2_formated=cell(0,0);
    untreatment_specific.DeltaEuc_formated=cell(0,0);
    untreatment_specific.proteinInCorum=cell(0,0);
    untreatment_specific.interactionInCorum=cell(0,0);
    
    %Create counter
    untreatment_row_counter=1;
    for treatment_counter=1:Total_unique_interactions
      if  untreatment_specific_interaction(treatment_counter,1)==1;
        untreatment_specific.unique_interactions(untreatment_row_counter,1)=interaction_final.unique_interactions(treatment_counter,1);
        untreatment_specific.proteinA(untreatment_row_counter,1)=interaction_final.proteinA(treatment_counter,1);
        untreatment_specific.proteinB(untreatment_row_counter,1)=interaction_final.proteinB(treatment_counter,1);
        untreatment_specific.centerA_formated(untreatment_row_counter,1)=interaction_final.centerA_formated(treatment_counter,1);
        untreatment_specific.centerB_formated(untreatment_row_counter,1)=interaction_final.centerB_formated(treatment_counter,1);
        untreatment_specific.Replicates_formated(untreatment_row_counter,1)=interaction_final.Replicates_formated(treatment_counter,1);
        %untreatment_specific.DeltaHeight_formated(untreatment_row_counter,1)=interaction_final.DeltaHeight_formated(treatment_counter,1);
        untreatment_specific.DeltaCenter_formated(untreatment_row_counter,1)=interaction_final.DeltaCenter_formated(treatment_counter,1);
        untreatment_specific.R2_formated(untreatment_row_counter,1)=interaction_final.R2_formated(treatment_counter,1);
        untreatment_specific.DeltaEuc_formated(untreatment_row_counter,1)=interaction_final.DeltaEuc_formated(treatment_counter,1);
        untreatment_specific.proteinInCorum(untreatment_row_counter,1)=interaction_final.proteinInCorum(treatment_counter,1);
        untreatment_specific.interactionInCorum(untreatment_row_counter,1)=interaction_final.interactionInCorum(treatment_counter,1);
        untreatment_row_counter=1+untreatment_row_counter;
      end
    end
    
    %Determine the number of interactions
    untreament_specific_interactionsnum=length(untreatment_specific.unique_interactions);
    
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
  end
  
  
  
  tic
  fprintf('        11. Write output files')
  writeOutput_ROC
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  tic
  fprintf('        12. Make figures')
  makeFigures_ROC
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  firstFlag = 0;
end



diary('off')


