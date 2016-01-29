%%%%%%%%%%%%%%% Instructions for running Gauss_Build.m:
%
%
%
%%%%%%%%%%%%%%% Logic:
%   Initialize
%   for replicate
%       read input
%       pre-process Gaussians, chroms. Make Dist.
%       make Protein structure
%       make TP, FP matrix from Corum
%       for Dist.i
%           figure ROC curve
%       for precision = 1:100
%           calculate precision as a function of Euc, Center
%       for precision = [50 60 70]
%           get the Euc,Dist for this precision
%           calculate recall, FN interactions
%           find binary interactions
%           combine binary interactions
%   for precision
%       combine interactions across replicates
%
%
%%%%%%%%%%%%%%% Custom functions called:
% -
%
%
%%%%%%%%%%%%%%% To do:
% - Why is this script pre-processing chromatograms? Shouldn't that have been done earlier?
% - What's going on in 8b? We keep interactions with Delta_Center>2... so why check for it?
% - FN_binary_interaction_list can get up to 3 GB. Can we make it smaller?
% - Both 5. and 6. calculate precision as a function of Euc distance. Can we combine them?
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


disp('ROC_PCPSILAC.m')



%% 0. Initialize
tic
fprintf('\n    0. Initialize')

% user defined variables
number_of_replicates=3;
number_of_channels=2;
replicate_counter = 4; % Do a single replicate. This is Nick's loop index.
%desiredPrecision = [.50 .60 .70]; % 50%, 60%, 70% precision
desiredPrecision = .70; % 70% precision
treatment_replicates=[4,5,6];
untreatment_replicates=[1,2,3];

% Define folders, i.e. define where everything lives.
maindir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/'; % where everything lives
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir = [maindir 'Data/']; % where data files live
datadir1 = [datadir 'ROC/']; % where data files live
datadir2 = [datadir 'ROC/tmp/']; % where data files live
datadir3 = [datadir 'ROC/CombinedResults/']; % where data files live
figdir1 = [maindir 'Figures/ROC/']; % where figures live
tmpdir = '/Users/Mercy/Academics/Foster/NickCodeData/4B_ROC_homologue_DB/Combined Cyto new DB/';
% Make folders if necessary
if ~exist(codedir, 'dir'); mkdir(codedir); end
if ~exist(funcdir, 'dir'); mkdir(funcdir); end
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist(datadir1, 'dir'); mkdir(datadir1); end
if ~exist(datadir2, 'dir'); mkdir(datadir2); end
if ~exist(datadir3, 'dir'); mkdir(datadir3); end
if ~exist(figdir1, 'dir'); mkdir(figdir1); end

% List all input files. These contain data that will be read by this script.
InputFile{1} = [datadir '/Major_protein_groups.xlsx'];
InputFile{2} = [datadir '/Corum_correctly_formated_Uniprot_IDs.csv'];
%InputFile{3} = [datadir 'Input/Corum_2012_human.xlsx'];
InputFile{4} = [datadir '/Master_guassian_list.csv'];

%define Raw SILAC ratios data, This is the output from the alignment script
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
List_of_Raw_filename={MvsL_filename_Raw_rep1,MvsL_filename_Raw_rep2,MvsL_filename_Raw_rep3,...
  HvsL_filename_Raw_rep1,HvsL_filename_Raw_rep2,HvsL_filename_Raw_rep3};
List_of_Gaus_filename={MvsL_filename_gaus_rep1,MvsL_filename_gaus_rep2,MvsL_filename_gaus_rep3,...
  HvsL_filename_gaus_rep1,HvsL_filename_gaus_rep2,HvsL_filename_gaus_rep3};

% pre-allocate final results
Recall = nan(number_of_replicates,number_of_channels);
Precision = nan((size(Recall)));
TPR = nan((size(Recall)));
FPR = nan((size(Recall)));
No_Int = nan((size(Recall)));

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


%%%%% Replicate counter starts
for replicate_counter = 1:number_of_replicates*number_of_channels
  
  s = ['\n        Replicate ' num2str(replicate_counter)];
  fprintf(s)
  
  %% 1. Read input data
  tic
  s = ['\n        1. Read input data, replicate ' num2str(replicate_counter)];
  fprintf(s)
  
  [~, Protein_IDs] = xlsread(InputFile{1});
  
  % Corum binary interactions
  fid=fopen(InputFile{2}, 'rt');    %Input corum data base.
  Corum_Import= textscan(fid,'%s\t', 'Delimiter',',');
  fclose(fid);
  No=length(Corum_Import{1})/2;
  Corum_Protein_names = (reshape(Corum_Import{1,1},2,No)');
  Unique_Corum = unique(Corum_Protein_names);
  
  Master_Gaussian_list=importdata(InputFile{4},',');
  
  %Raw data for analysis of euclidean distance
  %import data files summary data
  Maxquant = fopen (List_of_Raw_filename{replicate_counter});
  Maxquant_raw = textscan(Maxquant, '%s', 'Delimiter',',');
  fclose(Maxquant);
  
  %Reshape data for use in analysis
  Raw_Dimension=size(Maxquant_raw{:});
  Raw_data_reshaped=reshape(Maxquant_raw{:},68,(Raw_Dimension(1)/68))';
  Raw_data_reshaped_size=size(Raw_data_reshaped);
  
  %divide up reshaped data in correct form
  Chromatograms_raw = Raw_data_reshaped(2:end,3:end);
  Chromatograms_raw = cellfun(@str2num,Chromatograms_raw);
  
  %Gaus data data for analysis of center analysis
  Gaus_import_name = fopen (List_of_Gaus_filename{replicate_counter});
  Gaus_import = textscan(Gaus_import_name, '%s', 'Delimiter',',');
  fclose(Gaus_import_name);
  
  %Reshape data for use in analysis
  Gaus_import_Dimension=size(Gaus_import{:});
  Gaus_import_reshaped=reshape(Gaus_import{:},7,(Gaus_import_Dimension(1)/7))';
  
  % do a bit of housekeeping
  clear Maxquant_raw Raw_data_reshaped Corum_Import

  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 2. Pre-process gaussians, chromatograms. Get Dist
  % - In chromatograms, replaces nans with 0.05
  % - Remove gaussians and chromatograms with C<5
  % - Make normalized gaussians, i.e. set height to 1
  tic
  fprintf('        2. Pre-process gaussians, chromatograms. Get Dist.')
  
  % Pre-process data to remove Gaussian above fraction five
  % Create array to use for comparsions
  H_raw = cellfun(@str2num,Gaus_import_reshaped(2:end,2));
  C_raw = cellfun(@str2num,Gaus_import_reshaped(2:end,3));
  W_raw = cellfun(@str2num,Gaus_import_reshaped(2:end,4));
  
  % Replace NaN values with 0.05
  Chromatograms_raw(isnan(Chromatograms_raw))= 0.05;
  
  % remove gaussian with centers below five
  Ibad = C_raw<5;
  H = H_raw(~Ibad);
  C = C_raw(~Ibad);
  W = W_raw(~Ibad);
  Chromatograms = Chromatograms_raw(~Ibad,:);
  
  % Make normalized gaussians
  x = 0.25:0.25:100;
  gaussian_fit=zeros(length(H),length(x));
  for ii = 1:length(H)
    gaussian_fit(ii,:) = 1*exp(-(x- C(ii)).^2 /2/(W(ii).^2));
  end
  
  %Calculate distances
  Dist.Euc = squareform(pdist(Chromatograms,'euclidean'));
  Dist.Center = squareform(pdist(C,'euclidean'));
  Dist.Height = squareform(pdist(H,'euclidean'));
  Dist.Width = squareform(pdist(W,'euclidean'));
  Dist.Gaussian_fits = squareform(pdist(gaussian_fit,'euclidean'));
  Dist.R2 = 1 - corr(Chromatograms').^2; % one minus R squared
  Dist.dtw = ones(size(Dist.R2))*100;
%   for ii = 1:size(Chromatograms,1)
%     if mod(ii,100)==0;disp(num2str(ii));end
%     for jj = 1:size(Chromatograms,1)
%       if Dist.R2(ii,jj)<0.4
%         Dist.dtw(ii,jj) = dtw(Chromatograms(ii,:)',Chromatograms(jj,:)');
%         %Dist.dtw(ii,jj) = dtw_old(Chromatograms(ii,:),Chromatograms(jj,:));
%       end
%     end
%   end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 3. Make Protein structure
  % This summarizes the proteins in our sample
  tic
  fprintf('        3. Make Protein structure.')
  
  % Combine Majority Protein ID with Proteins Identified in replicate
  Protein.Isoform=Gaus_import_reshaped(2:end,1);
  Protein.Height=cellfun(@str2num,Gaus_import_reshaped(2:end,2));
  Protein.Width=cellfun(@str2num,Gaus_import_reshaped(2:end,3));
  Protein.Center=cellfun(@str2num,Gaus_import_reshaped(2:end,4));
  
  % remove protein name with center less then 5
  Protein.Isoform(Ibad)=[];
  Protein.Height(Ibad)=[];
  Protein.Width(Ibad)=[];
  Protein.Center(Ibad)=[];
  
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
  Unique_MajorProteinIDs = unique(vertcat(Protein.MajorID_NoIsoforms{:}),'rows');
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
  Corum_Dataset_expanded = cell(size(Corum_Dataset));
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
  
  
  %% 5. Make ROC curves for different Dist fieldnames
  % Aside from making FPR_out and Precision_out, this can be skipped
  tic
  fprintf('        5. Make figure for ROC curves.')
  
  if 0
  fn = fieldnames(Dist);
  for i=[1 6]
    
    count=0;
    nr=100;
    Recall_out.(fn{i})=zeros(nr,1);
    Precision_out.(fn{i})=zeros(nr,1);
    TPR_out.(fn{i})=zeros(nr,1);
    FPR_out.(fn{i})=zeros(nr,1);
    
    % dummy vectors
    dist = Dist.(fn{i})(possibleInts(:));
    tp = TP_Matrix(possibleInts(:));
    
    xDist = linspace(0,max(Dist.(fn{i})(:)),nr);    
    for Distance = xDist %Loop to calculate performance of variables  
      count=1+count;
      
      TP = sum(dist<Distance & tp);
      FP = sum(dist<Distance & ~tp);
      TN = sum(dist>Distance & ~tp);
      FN = sum(dist>Distance & tp);

      Recall = TP/(TP+FN);
      Precision = TP/(TP+FP);
      TPR = TP/(TP+FN);
      FPR = FP/(FP+TN);
      
      Recall_out.(fn{i})(count,1)= Recall;
      Precision_out.(fn{i})(count,1)=Precision;
      TPR_out.(fn{i})(count,1)=TPR;
      FPR_out.(fn{i})(count,1)=FPR;
    end
  end
  
  
  % plot all this
  if 1
    cols = {'b' 'g' 'c' [.5 .5 .5] 'r' 'm' [1 .5 0]};
    
    figure,hold on
    %for i=1:length(fn)
    for i=[1 6]
      plot(FPR_out.(fn{i}),TPR_out.(fn{i}),'color',cols{i})
    end
    plot([0 1],[0 1],'--k')
    axis([0 1 0 1])
    %legend('Euc','Center','Height','Width','Gauss fits','R2','location','northwest')
    legend('Euc','R2','location','northeast')
    xlabel('FPR')
    ylabel('TPR')
    
    figure,hold on
    %for i=1:length(fn)
    for i=[1 6]
      plot(Recall_out.(fn{i}),Precision_out.(fn{i}),'color',cols{i})
    end
    %axis([0 1 0 1])
    %legend('Euc','Center','Height','Width','Gauss fits','R2','location','northeast')
    legend('Euc','R2','location','northeast')
    xlabel('Recall')
    ylabel('Precision')
  end
  end
  
  % do a bit of housekeeping
  clear dist tp
  clear Int2_matrix Int3_matrix Int5_matrix Int6_matrix FP_Matrix_Int1 FP_Matrix_Int2 FP_Matrix_Int2_Triu
  clear MaxTP_matrix1 MaxTP_matrix2 TP_Matrix_Int2_Triu TP_Matrix_Int1 TP_Matrix_Int2 TP_Matrix_Int2_Triu
  clear Possible_Int1 Possible_Int2 Possible_Int2_Triu Total_Int Total_Int2_Triu
  clear Int6_matrix_Triu Int_matrix_minus_TP Inverse_TP_matrix Inverse_TP_matrix1 MaxTP_matrix2_Triu
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 6. Calculate score, the likelihood that a protein pair is interacting
  tic
  fprintf('        6. Calculate score, the likelihood that a protein pair is interacting')
  
  % Predict interactions through Euc
  %scoreMatrix = 1 - Dist.Euc;
  %score0 = scoreMatrix(possibleInts(:));
  
  % Predict interactions through R^2
  %scoreMatrix = 1 - Dist.R2;
  
  % Predict interactions with a classifier
  %scoreMatrix = scoresvm(Dist,possibleInts,TP_Matrix);
  %scoreMatrix = nanmean(scoreMatrix,2);
  %scoreMatrix = reshape(scoreMatrix,size(Dist.R2,1),size(Dist.R2,1));
  scoreMatrix = scorenb(Dist,possibleInts,TP_Matrix);
  scoreMatrix = nanmean(scoreMatrix,2);
  scoreMatrix = reshape(scoreMatrix,size(Dist.R2,1),size(Dist.R2,1));
  
  score = scoreMatrix(possibleInts(:));
        
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 7. Calculate precision as a function of score
  % Algorithm: Calculate precision as a function of score very coarsely. Identify at what score
  % value precision crosses the desired level. Zoom in on that score value. Iterate.
  tic
  fprintf('        7. Calculate precision as a function of score')
  
  class = TP_Matrix(possibleInts(:));
  nn = 25;
  ds = linspace(min(score),max(score),nn); % start off with coarse search
  Tol = 0.001; % get within 0.1% precision
  maxIter = 20; % zoom in 20 times at most
  
  xcutoff = nan(size(desiredPrecision));
  calcprec = zeros(size(xcutoff));
  calcrec = zeros(size(xcutoff));
  scoreRange = nan(nn*maxIter*length(desiredPrecision),1);  % Save these for
  precRange = nan(size(scoreRange));                        % plotting
  recRange = nan(size(scoreRange));                         % later.
  for di = 1:length(desiredPrecision)
    calcTol = 10^10;
    iter = 0;
    deltaPrec = 1;
    prec0 = zeros(nn,1);
    % Stop zooming in when one of three things happens:
    % i) you get close enough to the desired precision (within calcTol)
    % ii) you've been through maxIter iterations
    % iii) zooming in stops being useful (precision changes by less than deltaPrec b/w iterations)
    while calcTol>Tol && iter<maxIter && deltaPrec>1e-3
      iter=iter+1;
      rec = nan(nn,1);
      prec = nan(nn,1);
      for dd =1:length(ds)
        TP = sum(score>ds(dd) & class==1);
        FP = sum(score>ds(dd) & class==0);
        FN = sum(score<ds(dd) & class==1);
        prec(dd) = TP/(TP+FP);
        rec(dd) = TP/(TP+FN);
      end
      deltaPrec = nanmean(abs(prec - prec0));
      
      % Save vectors for plotting
      I = (iter-1)*nn+1 : iter*nn;
      scoreRange(I) = ds;
      precRange(I) = prec;
      recRange(I) = rec;
      
      % Calculate how close to desiredPrecision(di) you got
      [calcTol,I] = min(abs(prec - desiredPrecision(di)));
            
      % Zoom in on region of interest
      i1 = find(prec>desiredPrecision(di));
      if isempty(i1);
        mx = max(score);
      else
        mx = ds(i1(1));
      end
      i2 = find(prec<desiredPrecision(di));
      if isempty(i2);
        mn = min(score);
      else
        mn = ds(i2(end));
      end
      ds = linspace(mn,mx,nn);
      prec0 = prec;
    end
    xcutoff(di) = ds(I);
    calcprec(di) = prec(I);
    calcrec(di) = rec(I);
  end
  [scoreRange,I] = sort(scoreRange);
  precRange = precRange(I);
  recRange = recRange(I);

  
%   % i. Make interactions_based_only_Ecldian_distance at 80% precision
%   xDist = linspace(0,max(Dist.Euc(:)),nr);
%   prec_thresh = 0.80;
%   Precision = Precision_out.Euc; % Precision_out.Euc
%   
%   %Test if the Euclidian distance ever gives a precision above 0.8
%   if max(Precision) < prec_thresh
%     idx=1;
%     goodDist = xDist(idx); %determines the Eucldian distance which gives a precision closest to 0.8
%   else    
%     % interpolate to find goodDist
%     p1 = Precision - prec_thresh;
%     segments = find(diff(sign(p1)) == -2);
%     segments = segments(end);
%     x(1) = xDist(segments);
%     x(2) = xDist(segments+1);
%     y(1) = p1(segments);
%     y(2) = p1(segments+1);
%     m = (y(2)-y(1))/(x(2)-x(1));
%     n = y(2) - m * x(2);
%     goodDist = -n/m;
%   end
%   
%   %Determine interactions within precision of >0.80 based on Eu distance
%   Interactions_detected_due_to_High_Eucldian_distance = Dist.Euc < goodDist;
%   Interactions_Based_only_Eucldian_distance = (Interactions_detected_due_to_High_Eucldian_distance & inverse_self & (Dist.Center < 2));
%   
%   %create matrix for Eu dis interactions with precision >0.8
%   Final_Interactions_Based_only_Eucldian_distance = Interactions_Based_only_Eucldian_distance;
%   
%   
%   % ii. Make optimisation_matrix_precision
%   Height_Distance = 20.0; % As height has very little effect on interaction determination, this has been set very large
%   Width_Distangce = 20.0; % As width has very little effect on interaction determination, this has been set very large
%   MatrixHeight = Dist.Height < Height_Distance;
%   MatrixWidth = Dist.Width < Width_Distangce;
%   %Determine interactions for a FPR of 0.5%
%   [~, idx] = min(abs(FPR_out.Euc - 0.05)); %determines the postion which gives a FPR of 0.5%
%   Eucldian_distance_with_0_05FPR = ceil(xDist(idx)); %determines the Eucldian distance which gives a FPR of 0.5%
%   %Define number of data points
%   Optimisation_Dimension = Eucldian_distance_with_0_05FPR/0.025;
%   %optimisation_matrix_number_interaction = zeros(Optimisation_Dimension,Optimisation_Dimension);
%   %optimisation_matrix_precision = zeros(Optimisation_Dimension,Optimisation_Dimension);
%   Euclidian_Distance_optimisation = 0.025:0.025:Eucldian_distance_with_0_05FPR;
%   Gaussian_fit_Distance_optimisation = 0.025:0.025:Eucldian_distance_with_0_05FPR;
%   if 0
%     for iii=1:Optimisation_Dimension
%       for iv=1:Optimisation_Dimension
%         NS1 = Euclidian_Distance_optimisation(iii); % e-6
%         NS2 = Gaussian_fit_Distance_optimisation(iv); % e-6
%         MatrixGaussian = Dist.Gaussian_fits < NS1; % e-2
%         Matrix_Distance = Dist.Euc < NS2; % e-2
%         Interactions_optimisation_counter = (MatrixGaussian & MatrixHeight & MatrixWidth & Matrix_Distance & inverse_self); % e-2
%         total_interaction_number = length(find(triu(Interactions_optimisation_counter==1))); % e-2
%         optimisation_matrix_number_interaction(iii,iv) = total_interaction_number;
%         %Calculate true positives
%         TP_Matrix_Int_optimisation = (Interactions_optimisation_counter & TP_Matrix);
%         TP_optimisation = length(find(TP_Matrix_Int_optimisation==1));
%         %Calculate false positives
%         Int2_matrix_optimisation = (Interactions_optimisation_counter & Int_matrix);
%         Int3_matrix_optimisation = (Int2_matrix_optimisation & inverse_self);
%         FP_optimisation = length(find(Int3_matrix_optimisation==1))-TP_optimisation;
%         %Precision
%         Precision = TP_optimisation/(TP_optimisation+FP_optimisation);
%         optimisation_matrix_precision(iii,iv) = Precision;
%       end
%     end
%   else
%     % Takes a while! Just load it for now
%     a = ['/Users/Mercy/Academics/Foster/NickCodeData/4B_ROC_homologue_DB/Combined Cyto new DB/Replicate' ...
%       num2str(replicate_counter) '/ROC analysis data/optimisation_matrix_interaction_numbers_replicate' ...
%       num2str(replicate_counter) '.csv'];
%     optimisation_matrix_number_interaction = importdata(a);
%     
%     b = ['/Users/Mercy/Academics/Foster/NickCodeData/4B_ROC_homologue_DB/Combined Cyto new DB/Replicate' ...
%       num2str(replicate_counter) '/ROC analysis data/optimisation_matrix_precision_numbers_replicate' ...
%       num2str(replicate_counter) '.csv'];
%     optimisation_matrix_precision = importdata(b);
%   end
  
%   % iii. Make Precision_replicate_level
%   minimum_precision=(ceil(min(min(optimisation_matrix_precision))*100))/100;
%   test_precision = fliplr(minimum_precision:0.01:1);
%   number_test_precisions = length(test_precision);
%   
%   Recall_replicate_level = zeros(1,number_test_precisions);
%   Precision_replicate_level = zeros(1,number_test_precisions);
%   TPR_replicate_level = zeros(1,number_test_precisions);
%   FPR_replicate_level = zeros(1,number_test_precisions);
%   No_Int_replicate_level = zeros(1,number_test_precisions);
%   for test_precision_counter = 1:number_test_precisions
%     %  test_precision_counter
%     Centre_Euc_require_precision = optimisation_matrix_precision > test_precision(test_precision_counter);
%     precision_matrix = nan(size(optimisation_matrix_number_interaction));
%     I = optimisation_matrix_precision > test_precision(test_precision_counter);
%     precision_matrix(I) = optimisation_matrix_number_interaction(I);
%     
%     %maximum_counter=Optimisation_Dimension; %define based on the optimisation_matrix number
%     maximum_counter = size(precision_matrix,1); % define this way, otherwise can go out of bounds
%     Interactions_based_one_summed_precisions = zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
%     parfor vi = 1:maximum_counter   % scan down the Euc row and find which values results the required precision and maximum number of interactions
%       Maximum_interactions_at_given_Center = max(precision_matrix(vi,:));
%       if ~isnan(Maximum_interactions_at_given_Center)
%         inds = find(precision_matrix(vi,:) == Maximum_interactions_at_given_Center);
%         [col] = inds(1); %this takes the first number of inds thus ensure the most conservative Center will be used
%         optimised_center = Gaussian_fit_Distance_optimisation(vi);
%         optimised_Euc = Euclidian_Distance_optimisation(col);
%         
%         %use optimisation code to maximum out
%         MatrixGaussian = Dist.Gaussian_fits < optimised_center;
%         Matrix_Distance = Dist.Euc < optimised_Euc;
%         
%         % determine final interaction matrices
%         Interactions__precisions = (MatrixGaussian & MatrixHeight & MatrixWidth & Matrix_Distance & inverse_self);
%         Interactions_based_one_summed_precisions = (Interactions_based_one_summed_precisions | Interactions__precisions);
%       end
%     end
%     
%     % determine final interaction matrices
%     Final_interactions = (Final_Interactions_Based_only_Eucldian_distance | Interactions_based_one_summed_precisions); %test just interactions on Euc dist
%     
%     % Calculate True possitives
%     TP_Matrix_Int = (Final_interactions & TP_Matrix);
%     TP= length(find(TP_Matrix_Int==1));
%     
%     % Calculate false negatives
%     MaxTP = length(find(TP_Matrix==1));
%     FN = MaxTP-TP;
%     
%     % Calculate false positives
%     Int2_matrix = (Final_interactions & Int_matrix);
%     Int3_matrix = (Int2_matrix & inverse_self);
%     FP = length(find(Int3_matrix==1))-TP;
%     
%     % Find True negatives
%     TN = (length(find(Int_matrix==1)))-FP;
%     
%     Recall_replicate_level(test_precision_counter) = TP/(TP+FN);
%     Precision_replicate_level(test_precision_counter) = TP/(TP+FP);
%     TPR_replicate_level(test_precision_counter) = TP/(TP+FN);
%     FPR_replicate_level(test_precision_counter) =FP/(FP+TN);
%     No_Int_replicate_level(test_precision_counter) = length(find(triu(Final_interactions==1)));
%   end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 8. Determine interactions at a pre-set precision
  % What variables need to be saved in the precision loop?
  % - Sorted_of_final_interactions_Protein_70pc
  % - Sorted_of_final_interactions_Center_70pc
  % - List_of_final_interactions_Protein_70pc
  
  tic
  fprintf('        8. Determine interactions at a pre-set precision.')
  
  for pri = 1:length(desiredPrecision)
    
    % pre-allocate
    clear Sorted_of_final_interactions_Protein Sorted_of_final_interactions_Center List_of_final_interactions_Protein
    
    fn_count = 0; % false negative interaction counter, used in 7b
    interaction_count = 0; % interaction counter, used in 7c
    
%     %Determine interactions within precision of >0.70 based on Eu distance
%     prec_thresh = desiredPrecision(pri);
%     
%     % Define precision required
%     %note This module has been to designed to use the settings which give the highest number of high precision interactions with the widest Euc
%     Result_matrix = nan(size(optimisation_matrix_number_interaction));
%     I = optimisation_matrix_precision > prec_thresh;
%     Result_matrix(I) = optimisation_matrix_number_interaction(I);
%     
%     s = ['\n            ' num2str(prec_thresh*100) ' %% precision'];
%     fprintf(s)
%     
%     % 7a. Get the Euc,Dist for this precision; calculate recall
%     tic
%     fprintf('\n            7a. Find binary interactions')
%     
%     %determine the Euc and Center which give the user defined precision
%     %maximum_counter = Optimisation_Dimension; %define based on the optimisation_matrix number
%     maximum_counter = size(Result_matrix,1);
%     Interactions_based_one_summed_precisions = zeros(Dimensions_Gaus_import,Dimensions_Gaus_import);
%     parfor vi=1:maximum_counter   % scan down the Euc row and find which values results the required precision and maximum number of interactions
%       Maximum_interactions_at_given_Center = max(Result_matrix(vi,:));
%       if ~isnan(Maximum_interactions_at_given_Center)
%         inds = find(Result_matrix(vi,:)==Maximum_interactions_at_given_Center);
%         [col] = inds(1); %this takes the first number of inds thus ensure the most conservative Center will be used
%         optimised_center = Gaussian_fit_Distance_optimisation(vi);
%         optimised_Euc = Euclidian_Distance_optimisation(col);
%         
%         %use optimisation code to maximum out
%         MatrixGaussian = Dist.Gaussian_fits < optimised_center;
%         Matrix_Distance = Dist.Euc < optimised_Euc;
%         
%         % determine final interaction matrices
%         Interactions__precisions = (MatrixGaussian & MatrixHeight & MatrixWidth & Matrix_Distance & inverse_self);
%         Interactions_based_one_summed_precisions = (Interactions_based_one_summed_precisions | Interactions__precisions);
%       end
%     end
%     clear MatrixGaussian Matrix_Distance MatrixHeight MatrixWidth
%     
%     tt = toc;
%     fprintf('  ...  %.2f seconds\n',tt)
    
    
    % 7b. Calculate recall and precision, negative interactions
    %Determine the false "binary" interaction list
    tic
    fprintf('\n            b. Calculate recall and precision, negative interactions')
    
    % determine final interaction matrices
    %Final_interactions = (Final_Interactions_Based_only_Eucldian_distance | Interactions_based_one_summed_precisions); %test just interactions on Euc dist
    Final_interactions = scoreMatrix>xcutoff(pri) & inverse_self;

    % Calculate True possitives
    TP_Matrix_Int = (Final_interactions & TP_Matrix);
    TP = length(find(TP_Matrix_Int==1));
    
    % Calculate false negatives
    MaxTP = length(find(TP_Matrix==1));
    FN = MaxTP-TP;
    
    % Calculate false positives
    Int2_matrix = (Final_interactions & Int_matrix);
    Int3_matrix = (Int2_matrix & inverse_self);
    FP = length(find(Int3_matrix==1))-TP;
    clear Int2_matrix Int3_matrix
    
    % Find True negatives
    TN = (length(find(Int_matrix==1)))-FP;
    
    %Generate negative matrix for final interaction analysis
    Negative1_matrix = (Int_matrix  & ~Final_interactions);
    Negative2_matrix = (Negative1_matrix & inverse_self);
    
    % final results
    Recall(replicate_counter,pri) = TP/(TP+FN);
    Precision(replicate_counter,pri) = TP/(TP+FP);
    TPR(replicate_counter,pri) = TP/(TP+FN);
    FPR(replicate_counter,pri) = FP/(FP+TN);
    No_Int(replicate_counter,pri) = length(find(triu(Final_interactions==1)));
    
    % Predicted negative interactions
    Protein_elution = strcat(num2str(C),'*',Protein.Isoform);
    U1 = triu(Negative2_matrix); %Take half of matrix therefor unique interactions
    [ia,ib] = find(U1 == 1); % protein1, protein2 indices in false negatives
    FNBinaryInteractions = cell(length(ia),3);
    for ri=1:length(ia)
      ri1 = ia(ri);
      ri2 = ib(ri);
      Int1=Protein_elution(ri1);
      Int2=Protein_elution(ri2);
      %if strcmp(Int1, Int2)
      if C(ri1)==C(ri2) && ri1==ri2 % make sure it's not a self interaction
      else
        fn_count = fn_count+1;
        %FNBinaryInteractions{fn_count,1} = strcat(Int1,'-',Int2); % Interactions
        %FNBinaryInteractions{fn_count,2} = Dist.Euc(ri1,ri2); % Delta_EucDist
        %FNBinaryInteractions{fn_count,3} = abs(C(ri1)-C(ri2)); % Delta_Center
        %FNBinaryInteractions{fn_count,4} = abs(H(ri1)-H(ri2)); % Delta_Height
        %FNBinaryInteractions{fn_count,5} = abs(W(ri1)-W(ri2)); % Delta_Width
        %FNBinaryInteractions{fn_count,6} = Int_matrix(ri1,ri2); % Both_proteins_Corum
        %FNBinaryInteractions{fn_count,7} = TP_Matrix(ri1,ri2); % Interaction_in_Corum
        FNBinaryInteractions{fn_count,1} = Protein.Isoform{ri1}; % Protein_interactions1
        FNBinaryInteractions{fn_count,2} = Protein.Isoform{ri2}; % Protein_interactions2
        FNBinaryInteractions{fn_count,3} = TP_Matrix(ri1,ri2); % Interaction_in_Corum
        %FNBinaryInteractions{fn_count,10} = C(ri1); % Protein_interactions_center1
        %FNBinaryInteractions{fn_count,11} = C(ri2); % Protein_interactions_center2
      end
      % sanity check
      if ~strcmp(Int1, Int2) && C(ri1)==C(ri2) && ri1==ri2
        disp('uh oh!')
      end
    end
    
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
    
    
    
    % 7c. Find "Binary" interactions
    % G: A-B, A-C, B-C, B-Q, ...
    % G: What happens with this output?
    % G: It's read in after the replicate loop and combined
    tic
    fprintf('            c. Find binary interactions')
    
    Protein_elution= strcat(num2str(C),'*',Protein.Isoform);
    [n,m]=size(Dist.Euc);
    U = triu(Final_interactions);    %Take half of matrix therefor unique interactions
    BinaryInteractions = cell(sum(U(:)==1),13);
    for vii=1:n
      for viii= 1:m
        if U(vii,viii)==1
          Int1=Protein_elution{vii};
          Int2=Protein_elution{viii};
          if strcmp(Int1, Int2)
          else
            interaction_count=1+interaction_count;
            BinaryInteractions{interaction_count,1} = strcat(Int1,'-',Int2); % Interactions
            BinaryInteractions{interaction_count,2} = Dist.Euc(vii,viii); % Delta_EucDist
            BinaryInteractions{interaction_count,3} = abs(C(vii)-C(viii)); % Delta_Center
            BinaryInteractions{interaction_count,4} = abs(H(vii)-H(viii)); % Delta_Height
            BinaryInteractions{interaction_count,5} = abs(W(vii)-W(viii)); % Delta_Width
            BinaryInteractions{interaction_count,6} = Int_matrix(vii,viii); % Both_proteins_Corum
            BinaryInteractions{interaction_count,7} = TP_Matrix(vii,viii); % Interaction_in_Corum
            BinaryInteractions{interaction_count,8} = Protein.Isoform{vii}; % Protein_interactions1
            BinaryInteractions{interaction_count,9} = Protein.Isoform{viii}; % Protein_interactions2
            BinaryInteractions{interaction_count,10} = C(vii); % Protein_interactions_center1
            BinaryInteractions{interaction_count,11} = C(viii); % Protein_interactions_center2
            BinaryInteractions{interaction_count,12} = replicate_counter; % Protein_interactions_center2
            BinaryInteractions{interaction_count,13} = scoreMatrix(vii,viii); % Interaction score
          end
        else
        end
      end
    end
    
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
    
    
    
    tic
    fprintf('            d. Combine binary interactions')
    % List of interactions
    % A - B,C,Q; B - A; C - A,Q,D; ...
    % G: What happens with this output?
    % G: Just written to file.
    
    List_of_final_interactions = Final_interactions;    %As we want the complete list do not use triu
    Protein_elution_1= strcat(num2str(C),'*',Protein.Isoform);
    [n,m]=size(Dist.Euc);
    row_counter=0;
    empty_cell=cell(1,1);
    for vii=1:n
      List_of_final_interactions_Protein.Guassian(vii,1)= Protein_elution_1(vii);
      List_of_final_interactions_Protein.Protein(vii,1)= Protein.Isoform(vii);
      List_of_final_interactions_Protein.Protein_NoIsoform(vii,1)= Protein.NoIsoform(vii);
      for MajorID_protein_counter = 1:dimension_Protein_MajorID_NoIsoforms(2)
        List_of_final_interactions_Protein.Protein_MajorID_NoIsoforms(vii,MajorID_protein_counter)= Protein.MajorID_NoIsoforms(vii,MajorID_protein_counter);
      end
      List_of_final_interactions_Protein.Center(vii,1)= C(vii);
      final_count=0;
      for viii= 1:m
        if List_of_final_interactions(vii,viii)== 1 %Determine if proteins interact
          final_count=1+final_count;
          List_of_final_interactions_Protein.Interaction_Gaussian(vii,final_count)=Protein_elution_1(viii);
          List_of_final_interactions_Protein.Interaction_Protein(vii,final_count)=Protein.Isoform(viii);
          List_of_final_interactions_Protein.Interaction_Protein_NoIsoform(vii,final_count)=Protein.NoIsoform(viii);
          List_of_final_interactions_Protein.Interaction_Center(vii,final_count)=C(viii);
          if  TP_Matrix(vii,viii)== 1 %Determine if interaction is known in Corum
            List_of_final_interactions_Protein.Interaction_known_in_corum(vii,final_count)= 1;
          else
            List_of_final_interactions_Protein.Interaction_known_in_corum(vii,final_count)= 0;
          end
        end
      end
      if final_count == 0
        List_of_final_interactions_Protein.Interaction_Gaussian(vii,1)= {NaN};
        List_of_final_interactions_Protein.Interaction_Protein(vii,1)={NaN};
        List_of_final_interactions_Protein.Interaction_Center(vii,1)= 0;
        List_of_final_interactions_Protein.Interaction_known_in_corum(vii,1)=0;
      end
    end
    
    %Remove values in List of final with NaN
    List_of_final_interactions_Protein.Interaction_Gaussian(cellfun(@(x) any(isnan(x)),List_of_final_interactions_Protein.Interaction_Gaussian)) = {''};
    List_of_final_interactions_Protein.Interaction_Protein(cellfun(@(x) any(isnan(x)),List_of_final_interactions_Protein.Interaction_Protein)) = {''};
    
    %Determine the dimension of List_of_final_interactions_Protein.Interaction_Protein array
    Maximum_observed_interactions=size(List_of_final_interactions_Protein.Interaction_Protein); % determine the dimension of the interaction array in structure
    
    %Count number of corum interactions and total interactions
    List_of_final_interactions_Protein.List_Interaction_Gaussians=cell(Maximum_observed_interactions(1),1);
    List_of_final_interactions_Protein.Total_Interaction=cell(Maximum_observed_interactions(1),1);
    for Interaction_number_counter1=1:Maximum_observed_interactions(1)
      Interaction_counter_for_protein=0;
      for Interaction_number_counter2=1:Maximum_observed_interactions(2)
        Guassian_to_combine(Interaction_number_counter2)=List_of_final_interactions_Protein.Interaction_Gaussian(Interaction_number_counter1,Interaction_number_counter2);
        if any(Guassian_to_combine{Interaction_number_counter2})== 1 % test if value is not zero
          String_name= strcat(Guassian_to_combine{Interaction_number_counter2},' ; ');
          List_of_final_interactions_Protein.List_Interaction_Gaussians{Interaction_number_counter1}=...
            strcat(String_name,List_of_final_interactions_Protein.List_Interaction_Gaussians{Interaction_number_counter1});
          Interaction_counter_for_protein=Interaction_counter_for_protein+1;
        end
      end
      List_of_final_interactions_Protein.Total_Interaction{Interaction_number_counter1,1}=Interaction_counter_for_protein;
    end
    
    %Sum number of known interactions in corum for protein
    for Interaction_number_counter1=1:Maximum_observed_interactions(1)
      for Interaction_number_counter2=1:Maximum_observed_interactions(2)
        Detected_Center_Interaction_known_in_corum(Interaction_number_counter2)= List_of_final_interactions_Protein.Interaction_known_in_corum(Interaction_number_counter1,Interaction_number_counter2);
      end
      List_of_final_interactions_Protein.Number_Interaction_known_in_corum(Interaction_number_counter1,1)= sum(Detected_Center_Interaction_known_in_corum);
    end
    
    %Protein list of interactions
    List_of_final_interactions_Protein.List_Interaction_Protein=cell(Maximum_observed_interactions(1),1);
    for protein_table1=1:Maximum_observed_interactions(1)
      for protein_table2=1:Maximum_observed_interactions(2)
        Proteins_to_combine(protein_table2)=List_of_final_interactions_Protein.Interaction_Protein(protein_table1,protein_table2);
        if any(Proteins_to_combine{protein_table2})== 1 % test if value is not zero
          String_name= strcat(Proteins_to_combine{protein_table2},' ; ');
          List_of_final_interactions_Protein.List_Interaction_Protein{protein_table1,1}=strcat(String_name,List_of_final_interactions_Protein.List_Interaction_Protein{protein_table1});
        end
      end
    end
    
    %Copy List_of_final_interactions_Protein_70pc to use for other experiuments
    % List_of_final_interactions_Protein_pc=List_of_final_interactions_Protein;
    %
    %   %Sort structure by Center/Gaussian write out table
    %   [~,location2]=sort([List_of_final_interactions_Protein.Center]); %Sort strcuture based on protein name
    %   for sort_counter1=1:Maximum_observed_interactions(1)
    %     Sorted_of_final_interactions_Center.Guassian{sort_counter1,1} = List_of_final_interactions_Protein.Guassian{location2(sort_counter1)};
    %     Sorted_of_final_interactions_Center.Protein{sort_counter1,1}= List_of_final_interactions_Protein.Protein{location2(sort_counter1)};
    %     Sorted_of_final_interactions_Center.Protein_NoIsoform{sort_counter1,1}=List_of_final_interactions_Protein.Protein_NoIsoform{location2(sort_counter1)};
    %     for MajorID_protein_counter = 1:dimension_Protein_MajorID_NoIsoforms(2)
    %       Sorted_of_final_interactions_Center.Protein_MajorID_NoIsoforms{sort_counter1,MajorID_protein_counter}=...
    %         List_of_final_interactions_Protein.Protein_MajorID_NoIsoforms{location2(sort_counter1),MajorID_protein_counter};
    %     end
    %     Sorted_of_final_interactions_Center.Center{sort_counter1,1}= List_of_final_interactions_Protein.Center(location2(sort_counter1));
    %     for sort_counter2=1:Maximum_observed_interactions(2)
    %       Sorted_of_final_interactions_Center.Interaction_Gaussian{sort_counter1,sort_counter2}= List_of_final_interactions_Protein.Interaction_Gaussian{location2(sort_counter1),sort_counter2};
    %       Sorted_of_final_interactions_Center.Interaction_Protein{sort_counter1,sort_counter2}=List_of_final_interactions_Protein.Interaction_Protein{location2(sort_counter1),sort_counter2};
    %       Sorted_of_final_interactions_Center.Interaction_Center{sort_counter1,sort_counter2}=List_of_final_interactions_Protein.Interaction_Center(location2(sort_counter1),sort_counter2);
    %       Sorted_of_final_interactions_Center.Interaction_known_in_corum{sort_counter1,sort_counter2}=List_of_final_interactions_Protein.Interaction_known_in_corum(location2(sort_counter1),sort_counter2);
    %     end
    %     Sorted_of_final_interactions_Center.Number_Interaction_known_in_corum(sort_counter1,1)= List_of_final_interactions_Protein.Number_Interaction_known_in_corum(location2(sort_counter1));
    %     Sorted_of_final_interactions_Center.Total_Interaction(sort_counter1,1)=List_of_final_interactions_Protein.Total_Interaction(location2(sort_counter1));
    %     Sorted_of_final_interactions_Center.List_Interaction_Gaussians(sort_counter1,1)=List_of_final_interactions_Protein.List_Interaction_Gaussians(location2(sort_counter1));
    %     Sorted_of_final_interactions_Center.List_Interaction_Protein(sort_counter1,1)=List_of_final_interactions_Protein.List_Interaction_Protein(location2(sort_counter1));
    %   end
    %
    %
    %   %Sort strcuture by Protein name
    %   [~,location]=sort([List_of_final_interactions_Protein.Protein]); %Sort strcuture based on protein name
    %   for sort_counter1=1:Maximum_observed_interactions(1)
    %     Sorted_of_final_interactions_Protein.Guassian{sort_counter1,1} = List_of_final_interactions_Protein.Guassian{location(sort_counter1)};
    %     Sorted_of_final_interactions_Protein.Protein{sort_counter1,1}= List_of_final_interactions_Protein.Protein{location(sort_counter1)};
    %     Sorted_of_final_interactions_Protein.Protein_NoIsoform{sort_counter1,1}=List_of_final_interactions_Protein.Protein_NoIsoform{location(sort_counter1)};
    %     for MajorID_protein_counter = 1:dimension_Protein_MajorID_NoIsoforms(2)
    %       Sorted_of_final_interactions_Protein.Protein_MajorID_NoIsoforms{sort_counter1,MajorID_protein_counter}=...
    %         List_of_final_interactions_Protein.Protein_MajorID_NoIsoforms{location(sort_counter1),MajorID_protein_counter};
    %     end
    %     Sorted_of_final_interactions_Protein.Center{sort_counter1,1}= List_of_final_interactions_Protein.Center(location(sort_counter1));
    %     for sort_counter2=1:Maximum_observed_interactions(2)
    %       Sorted_of_final_interactions_Protein.Interaction_Gaussian{sort_counter1,sort_counter2}= List_of_final_interactions_Protein.Interaction_Gaussian{location(sort_counter1),sort_counter2};
    %       Sorted_of_final_interactions_Protein.Interaction_Protein{sort_counter1,sort_counter2}=List_of_final_interactions_Protein.Interaction_Protein{location(sort_counter1),sort_counter2};
    %       Sorted_of_final_interactions_Protein.Interaction_Center{sort_counter1,sort_counter2}=List_of_final_interactions_Protein.Interaction_Center(location(sort_counter1),sort_counter2);
    %       Sorted_of_final_interactions_Protein.Interaction_known_in_corum{sort_counter1,sort_counter2}=List_of_final_interactions_Protein.Interaction_known_in_corum(location(sort_counter1),sort_counter2);
    %     end
    %     Sorted_of_final_interactions_Protein.Number_Interaction_known_in_corum(sort_counter1,1)= List_of_final_interactions_Protein.Number_Interaction_known_in_corum(location(sort_counter1));
    %     Sorted_of_final_interactions_Protein.Total_Interaction(sort_counter1,1)=List_of_final_interactions_Protein.Total_Interaction(location(sort_counter1));
    %     Sorted_of_final_interactions_Protein.List_Interaction_Gaussians(sort_counter1,1)=List_of_final_interactions_Protein.List_Interaction_Gaussians(location(sort_counter1));
    %     Sorted_of_final_interactions_Protein.List_Interaction_Protein(sort_counter1,1)=List_of_final_interactions_Protein.List_Interaction_Protein(location(sort_counter1));
    %   end
    
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
    
    
    
    
    tic
    s = ['            e. Write binary interactions, rep' num2str(replicate_counter) ' pc' num2str(round(desiredPrecision(pri)*100))];
    fprintf(s)
    
    p = round(desiredPrecision(pri)*100);
    sf = [datadir2 'Interaction_list_pc' num2str(p) '_rep' num2str(replicate_counter) '.mat'];
    save(sf,'FNBinaryInteractions','BinaryInteractions','List_of_final_interactions_Protein')
    clear FNBinaryInteractions BinaryInteractions List_of_final_interactions_Protein
    
    tt = toc;
    fprintf('  ...  %.2f seconds\n',tt)
  end
  
  xcut{replicate_counter} = xcutoff
  
end
%%%%% Replicate counter ends

% do a bit of housekeeping
clear Dist Int_matrix


%% 8. Combine interactions across replicates
% The only files Nick reads in here are:
% - Interactions_list_70pc_precision_replicate1.csv
%    fprintf(fid11B,'%6.3f,%s,%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,\n',replicate_counter,...
%      Protein_interactions{ix,1},Protein_interactions{ix,2},Protein_interactions_center(ix,1),Protein_interactions_center(ix,2),...
%      Delta_Height(ix),Delta_Center(ix), Delta_Width(ix), Delta_EucDist(ix),...
%      Both_proteins_Corum(ix), Interaction_in_Corum(ix));
%   --> This is now in the BinaryInteractions structure
% - Neg_proteins_interactions_70pc_replicate_1.csv
%    fprintf(fid_Neg,'%6f,%s,%s,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,\n',...
%      replicate_counter(1) ,FN_Protein_interactions1, FN_Protein_interactions2,...
%      FN_Protein_interactions_center1, FN_Protein_interactions_center2, FN_Delta_Height,...
%      FN_Delta_Center, FN_Delta_Width, FN_Delta_EucDist, FN_Both_proteins_Corum, FN_Interaction_in_Corum);
%   --> This is now in the FNBinaryInteractions structure

tic
fprintf('    8. Make final interaction list')

binary_interaction_list = [];
FN_binary_interaction_list = [];

for pri = 1:length(desiredPrecision)
  prec_thresh = desiredPrecision(pri);
  s = ['\n       ' num2str(prec_thresh*100) ' %% precision'];
  fprintf(s)
  
  fprintf('\n       8a. Combine interactions across replicates')
  
  for replicate_counter = 1:(number_of_replicates*number_of_channels)
    % read in binary interactions file
    % contains FNBinaryInteractions, BinaryInteractions, List_of_final_interactions_Protein
    sf = [datadir2 'Interaction_list_pc' num2str(p) '_rep' num2str(replicate_counter) '.mat'];
    load(sf)
    
    % concatenate interactions files
    %True interactions
    binary_interaction_list = vertcat(binary_interaction_list,BinaryInteractions);
    %False interactions
    FN_binary_interaction_list = vertcat(FN_binary_interaction_list,FNBinaryInteractions);
    
    clear FNBinaryInteractions BinaryInteractions List_of_final_interactions_Protein
  end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  tic
  fprintf('       8b. Create a list of all the unique interactions')
  
  interaction_pairs = cell(size(binary_interaction_list,1),1);
  for ri=1:size(binary_interaction_list,1)
    interaction_pairs{ri} = [binary_interaction_list{ri,8} '_' binary_interaction_list{ri,9}];
  end
  unique_interaction = unique(interaction_pairs);
  length_unique_inter = length(unique_interaction);
  
  %Create array to store values
  interaction_final.unique_interactions=cell(0,0);
  interaction_final.replicate_numbers=cell(0,0);
  interaction_final.proteinA=cell(0,0);
  interaction_final.proteinB=cell(0,0);
  interaction_final.CenterA=cell(0,0);
  interaction_final.CenterB=cell(0,0);
  interaction_final.DeltaHeight=cell(0,0);
  interaction_final.DeltaCenter=cell(0,0);
  interaction_final.DeltaWidth=cell(0,0);
  interaction_final.DeltaEucDist=cell(0,0);
  interaction_final.proteinInCorum=cell(0,0);
  interaction_final.interactionInCorum=cell(0,0);
  Unique_interaction_counter=0;
  for interaction_counter2A = 1:length_unique_inter
    %Find location of unique protein interaction
    %location_interaction_pairs = ind2sub(length_unique_inter, strmatch(unique_interaction(interaction_counter2A), interaction_pairs(2:end), 'exact'));
    location_interaction_pairs = strmatch(unique_interaction(interaction_counter2A), interaction_pairs(1:end), 'exact');
    
    %find which replicate the interactions were seen in
    %while ~isempty(location_interaction_pairs(:))
    
    diffC = abs(cellfun(@minus,binary_interaction_list(location_interaction_pairs,10),...
      binary_interaction_list(location_interaction_pairs,11)));
    position_within_two = find(diffC<100);
    
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
        interaction_final.DeltaHeight(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),4);
        interaction_final.DeltaCenter(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),3);
        interaction_final.DeltaWidth(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),5);
        interaction_final.DeltaEucDist(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),2);
        interaction_final.score(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(position_within_two(row)),13);
        
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
      interaction_final.DeltaHeight(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),4);
      interaction_final.DeltaCenter(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),3);
      interaction_final.DeltaWidth(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),5);
      interaction_final.DeltaEucDist(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),2);
      interaction_final.proteinInCorum(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),6);
      interaction_final.interactionInCorum(Unique_interaction_counter,1) = binary_interaction_list(location_interaction_pairs(1),7);
      interaction_final.score(Unique_interaction_counter,row) = binary_interaction_list(location_interaction_pairs(1),13);
    end
  end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  mySound
  % Calculate TP, FP, FN and TN
  tic
  fprintf('       8c. Calculate TP, FP, FN and TN')
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
  interaction_final.times_observed_across_reps = number_unique_interaction;
  
  
  %Determine how many times an interaction was observed
%   for replicate_counter=1:max(number_observation_pre_interaction)
%     number_observation_pre_replicate(replicate_counter)=nnz(find(number_observation_pre_interaction(:)==replicate_counter));
%   end
  
  %Determine how many times an interaction was observed
%   for replicate_counter=1:max(number_unique_interaction)
%     number_unique_observation_pre_replicate(replicate_counter)=nnz(find(number_unique_interaction(:)==replicate_counter));
%   end
  
  %total number of interactions of unique replicates
  total_observation=sum(number_observation_pre_interaction);
  %total number of interactions of unique replicates
  total_unique_observation=sum(number_unique_interaction);
  
  
  % create unique identifier for false interactions
  % ##### THIS CODE IS AN APPROXIMATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  % ##### If protein A and B are not interacting in corum, this sets ALL THE GAUSSIAN COMBINATIONS
  %       between A and B to not interacting. See Dec 3 notes for more info.
  %       The proper fix is to individually treat all Gaussian combinations exactly like we do for
  %       the positive interactions.
  interaction_FN_pairs = cell(size(FN_binary_interaction_list,1),1);
  for ri=1:size(FN_binary_interaction_list,1)
    interaction_FN_pairs{ri} = [FN_binary_interaction_list{ri,1} '_' FN_binary_interaction_list{ri,2}];
  end
  %interaction_FN_pairs=strcat(FN_binary_interaction_list(:,2),'_',FN_binary_interaction_list(:,3));
  unique_FN_interaction = unique(interaction_FN_pairs(2:end));
  length_unique_FN_inter = length(unique_FN_interaction);
  
  %Determine the number of unique false negatives
  interaction_as_number = cat(1,FN_binary_interaction_list{:,3});
  interaction_FN_incorum_pairs = interaction_FN_pairs(interaction_as_number==1);
  length_unique_FN_incorum = length(unique(interaction_FN_incorum_pairs));
    
  %clear varibles to speed up processing
  clearvars interaction_FN_incorum_pairs;
  clearvars interaction_FN_pairs;
  clearvars FN_binary_interaction_list;
  clear interaction_as_number
  
  % Caliculate TP, FP, FN and TN
  %Caliculate values from true interaction data set
  True_proteinInCorum=sum([cell2mat(interaction_final.proteinInCorum)]);
  True_interactionInCorum=sum([cell2mat(interaction_final.interactionInCorum)]);
  
  %Caliculate values from true interaction data set
  False_proteinInCorum=length_unique_FN_inter;
  False_interactionInCorum=length_unique_FN_incorum;
  
  %caliculate TP, FP, FN and TN
  final_TN=False_proteinInCorum-False_interactionInCorum;
  final_FP=True_proteinInCorum-True_interactionInCorum;
  final_FN=False_interactionInCorum;
  final_TP=True_interactionInCorum;
  
  % final results
  final_Recall = final_TP/(final_TP+final_TN);
  final_Precision = final_TP/(final_TP+final_FP);
  final_TPR = final_TP/(final_TP+final_FN);
  final_FPR =final_FP/(final_FP+final_TN);
  
  
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
  
  %Crete list of DeltaHeight, note ensure strjoin function is avalible
  interaction_final.DeltaHeight_formated=cell(Total_unique_interactions,1);
  
  for format_loop=1:Total_unique_interactions
    DeltaHeight_2bformated1=interaction_final.DeltaHeight(format_loop,:);
    DeltaHeight_2bformated1(cellfun('isempty',DeltaHeight_2bformated1)) = [];
    length_DeltaHeight=length(DeltaHeight_2bformated1);
    
    %Test if the array is longer then one entry
    if length_DeltaHeight<2
      interaction_final.DeltaHeight_formated(format_loop)=DeltaHeight_2bformated1(1);
    elseif length_DeltaHeight>=2
      for jj= 1:length(DeltaHeight_2bformated1)
        DeltaHeight_2bformated1{jj} = num2str(DeltaHeight_2bformated1{jj});
      end
      interaction_final.DeltaHeight_formated(format_loop)= {strjoin(DeltaHeight_2bformated1,' ; ')};
    end
  end
  
  %Crete list of Deltawidth, note ensure strjoin function is avalible
  interaction_final.Deltawidth_formated=cell(Total_unique_interactions,1);
  for format_loop=1:Total_unique_interactions
    Deltawidth_2bformated1=interaction_final.DeltaWidth(format_loop,:);
    Deltawidth_2bformated1(cellfun('isempty',Deltawidth_2bformated1)) = [];
    length_Deltawidth=length(Deltawidth_2bformated1);
    %Test if the array is longer then one entry
    if length_Deltawidth<2
      interaction_final.Deltawidth_formated(format_loop)=Deltawidth_2bformated1(1);
    elseif length_Deltawidth>=2
      for jj= 1:length(Deltawidth_2bformated1)
        Deltawidth_2bformated1{jj} = num2str(Deltawidth_2bformated1{jj});
      end
      interaction_final.Deltawidth_formated(format_loop)= {strjoin(Deltawidth_2bformated1,' ; ')};
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
    position_to_sum=find(interaction_final.times_observed_across_reps==precision_counter);
    
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
  
  
  
  tic
  fprintf('       8d. Interactions observed only under treatment')
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
  treatment_specific.DeltaHeight_formated=cell(0,0);
  treatment_specific.DeltaCenter_formated=cell(0,0);
  treatment_specific.Deltawidth_formated=cell(0,0);
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
      treatment_specific.DeltaHeight_formated(treatment_row_counter,1)=interaction_final.DeltaHeight_formated(treatment_counter,1);
      treatment_specific.DeltaCenter_formated(treatment_row_counter,1)=interaction_final.DeltaCenter_formated(treatment_counter,1);
      treatment_specific.Deltawidth_formated(treatment_row_counter,1)=interaction_final.Deltawidth_formated(treatment_counter,1);
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
  untreatment_specific.DeltaHeight_formated=cell(0,0);
  untreatment_specific.DeltaCenter_formated=cell(0,0);
  untreatment_specific.Deltawidth_formated=cell(0,0);
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
      untreatment_specific.DeltaHeight_formated(untreatment_row_counter,1)=interaction_final.DeltaHeight_formated(treatment_counter,1);
      untreatment_specific.DeltaCenter_formated(untreatment_row_counter,1)=interaction_final.DeltaCenter_formated(treatment_counter,1);
      untreatment_specific.Deltawidth_formated(untreatment_row_counter,1)=interaction_final.Deltawidth_formated(treatment_counter,1);
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
  
  
  
  tic
  fprintf('       8e. Write output files')
  writeOutput_ROC
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  tic
  fprintf('       8f. Make figures')
  makeFigures_ROC
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
end






