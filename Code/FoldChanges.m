%FOLDCHANGES Calculate changes in protein amount.
%   For datasets with multiple conditions, FOLDCHANGES calculates the ratio
%   of protein amount in condition A to protein amount in condition B, i.e. 
%   a fold-change in protein amount between conditions. If there are more 
%   than two conditions, condition pairs can be specified in master script 
%   PRINCE. For each protein, fold-changes are calculated using data
%   centered in a window around peaks identified by GAUSSBUILD.
%
%   For datasets with multiple replicates, FOLDCHANGES calculates the
%   significance of any fold changes via a paired t-test.
%
%   FOLDCHANGES is an optional module and can be turned on or off in the
%   master script PRINCE. If run, FOLDCHANGES must be run after PRINCE and 
%   GAUSSBUILD.
%
%   FOLDCHANGES produces two output folders: Output/Data/FoldChanges, which 
%   contains all output csv tables, and Output/Figures/FoldChanges, which
%   contains all figures.
%
%   Workflow:
%   1. Read co-fractionation profiles and GaussBuild output
%   2. Calculate fold changes
%   3. Collate variables, house clean
%   4. Calculate significance of fold changes
%   5. Write output tables
%   6. Make figures
%
%   See also PRINCE, GAUSSBUILD.


diary([user.maindir 'logfile.txt'])
disp('FoldChanges.m')

skipflag = 0;
if length(user.silacratios)==1
  disp('* NB: User set number of SILAC ratios to 1. Skipping FoldChanges...')
  user.skipcomparison = 1;
  skipflag = 1;
end
if user.skipcomparison==1 && skipflag==0
  disp('* NB: User set skipcomparison to True. Skipping FoldChanges...')
  skipflag = 1;
end

% check toolboxes
checktoolbox;

if ~skipflag
  
  %% 0. Initialize
  tic
  fprintf('\n    0. Initialize')
  
  
  % Load user settings
  maindir = user.maindir;
  User_Window = user.userwindow;
  Experimental_channels = user.silacratios;
  User_alignment_window1 = user.userwindow; % User defined window to consider Guassian the same
  Nchannels = length(Experimental_channels);
  %Diltuion_factor_master_mix = user.Dilution_factor;
  fraction_to_plot = user.Nfraction;
  
  
  % Define folders, i.e. define where everything lives.
  codedir = [maindir 'Code/']; % where this script lives
  funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
  datadir = [maindir 'Output/Data/FoldChanges/']; % where data files live
  figdir = [maindir 'Output/Figures/FoldChanges/']; % where figures live
  % Make folders if necessary
  if ~exist([maindir '/Output'], 'dir'); mkdir([maindir '/Output']); end
  if ~exist([maindir '/Output/Data'], 'dir'); mkdir([maindir '/Output/Data']); end
  if ~exist([maindir '/Output/Figures'], 'dir'); mkdir([maindir '/Output/Figures']); end
  if ~exist([maindir '/Output/tmp'], 'dir'); mkdir([maindir '/Output/tmp']); end
  if ~exist(figdir, 'dir'); mkdir(figdir); end
  if ~exist(datadir, 'dir'); mkdir(datadir); end
  if ~exist([figdir '/IndividualProteins'], 'dir'); mkdir([figdir '/IndividualProteins']); end
  
  % Find input files
  ChromatogramIn = cell(Nchannels,1);
  GaussIn = cell(Nchannels,1);
  if user.nickflag==1
    for ii = 1:Nchannels
      ChromatogramIn{ii} = [maindir '/Output/tmp/' 'Adjusted_' user.silacratios{ii} '_Raw_data_maxquant_modified.xlsx']; % from Alignment
      GaussIn{ii} = [maindir '/Output/tmp/' 'Adjusted_' user.silacratios{ii} '_Combined_OutputGaus.csv']; % from Alignment
    end
  else
    dd = dir([maindir '/Output/tmp/Adjusted_*_Combined_OutputGaus.csv']);
    if user.skipalignment==1 || isempty(dd)
      % If Alignment was skipped, use raw data + Gauss_Build output
      
      for di = 1:Nchannels
        ChromatogramIn{di} = user.MQfiles{di};
      end
      
      dd = dir([maindir '/Output/Data/GaussBuild/*_Combined_OutputGaus.csv']);
      for di = 1:length(dd)
        GaussIn{di} = [maindir '/Output/Data/GaussBuild/' dd(di).name];
      end
    else
      % If Alignment was not skipped, use Alignment output
      
      dd = dir([maindir '/Output/tmp/Adjusted_*_Raw_data_maxquant.csv']);
      for di = 1:length(dd)
        ChromatogramIn{di} = [maindir '/Output/tmp/' dd(di).name];
      end
      
      dd = dir([maindir '/Output/tmp/Adjusted_*_Combined_OutputGaus.csv']);
      for di = 1:length(dd)
        GaussIn{di} = [maindir '/Output/tmp/' dd(di).name];
      end
    end
  end
  
  
  % Plotting variables
  myC= [30/255 144/255 255/255
    255/255 215/255 0/255
    178/255 34/255 34/255
    193/255 205/255 193/255];
  colour_to_use=[0.254 0.411 0.882 %Colour: Royal Blue
    135/255 206/255 250/255 %Colour: Sky Blue
    0.28 0.82 0.8 %Colour: Medium Turquoise
    205/255 92/255 92/255 %Colour: Indian Red
    178/255 34/255 34/255 %Colour: Firebrick
    255/255 69/255 0 %Colour: Orange Red
    244/255 238/255 224/255 %Colour: Honeydew 2
    106/255 90/255 205/255 %Colour: Slate Blue
    34/255 139/255 34/255 %Colour: Forest Green
    222/255 184/255 135/255 %Colour: Burlywood
    186/255 85/255 211/255 %Colour: Medium Orchid
    219/255 112/255 147/255]; %Colour: Pale Violet Red
  
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 1. Read input
  tic
  fprintf('    1. Read input')
  
  % Import chromatogram data
  num_val = cell(Nchannels,1);
  txt_val = cell(Nchannels,1);
  for ii = 1:Nchannels
    %     tmp = importdata(ChromatogramIn{ii}); % from Alignment
    %     if isfield(tmp.data,'Sheet1')
    %       num_val{ii} = tmp.data.Sheet1;
    %       txt_val{ii} = tmp.textdata.Sheet1;
    %     else
    %       num_val{ii} = tmp.data;
    %       txt_val{ii} = tmp.textdata;
    %     end
    %
    
    fid = fopen(ChromatogramIn{ii});
    cc = 0;
    txt_val{ii} = cell(10^6,1);
    % put header in the first row of txt_val
    t = fgetl(fid);
    txt_val{ii}{1} = t;
    while ~feof(fid)
      cc = cc+1;
      t = fgetl(fid);
      t1 = strsplit(t,',');
      if cc ==1
        num_val{ii} = nan(10^6,length(t1)-1);
      end
      
      txt_val{ii}{cc+1} = t1{1};
      for kk = 1:size(num_val{ii},2)
        num_val{ii}(cc,kk) = str2double(t1{kk+1});
      end
    end
    txt_val{ii} = txt_val{ii}(1:cc+1);
    num_val{ii} = num_val{ii}(1:cc,:);
    
    % if txt_val includes protein groups, reduce it to the first protein in each group
    for jj = 1:size(txt_val{ii},1)
      tmp = strsplit(txt_val{ii}{jj},';');
      txt_val{ii}{jj} = tmp{1};
    end

    % Confirm that the 'Replicate' column is in the header
    Ihead = strfind(lower(txt_val{ii}(1,:)),'replicate');
    %replicate_in_header = find(~cellfun('isempty', Ihead));
    if isempty(Ihead)
      error('Comparison: Chromatogram numerical and text data mismatch.')
    end
    
    % Confirm that the first column of num_val is the replicate
    Nempty = sum(isempty(num_val{ii}(:,1)));
    if length(unique(num_val{ii}(:,1)))==user.Nreplicate && Nempty<2
      % num_val{ii}(:,1) is the replicate number
    else
      I = find(ismember(lower(txt_val{ii}(1,:)),'replicate'));
      replicate_cell = txt_val{ii}(:,I);
      replicate_vector = nan(size(replicate_cell));
      for jj = 1:length(replicate_cell)
        rep = str2num(replicate_cell{jj});
        if ~isempty(rep)
          replicate_vector(jj) = rep;
        end
      end
      
      if length(replicate_vector) == size(num_val{ii},1)
        num_val{ii} = [replicate_vector num_val{ii}];
      elseif length(replicate_vector) == size(num_val{ii},1)+1
        num_val{ii} = [replicate_vector(2:end) num_val{ii}];
      end
    end

    
    % Remove header from txt_val if necessary
    if size(txt_val{ii},1) == size(num_val{ii},1)+1
      txt_val{ii} = txt_val{ii}(2:end,:);
    elseif size(txt_val{ii},1) == size(num_val{ii},1)
      % do nothing
    else
    end
    
    %Remove data point with low values
    num_val{ii}(num_val{ii}<0.2) = nan;
    
    % Summary numbers: proteins, fractions, replicates
    Nproteins = length(num_val{1});
    Nfraction = size(num_val{1},2);
    Nfraction = Nfraction-1;
    a = num_val{1}(:,1);
    a(isnan(a)) = [];
    replicate_num = length(unique(a));
    
    % The data is nan-padded. Find where the real data starts and stops.
    nan_max = size(num_val{ii},1);
    tmp = find(sum(isnan(num_val{ii}))==size(num_val{ii},1));
    if isempty(tmp)
      tmp = -1;
    end
    frac1 = max([2 tmp(find(tmp<Nfraction/2,1,'last'))]); % start of real data
    frac2 = min([size(num_val{1},2) tmp(find(tmp>Nfraction/2,1,'first'))]); % end of real data
    
    % Ensure txt_val is a single column of protein names
    txt_val{ii} = txt_val{ii}(:,1);
    
    % Add unique identifiers to the chromatograms
    Unique_indentifer_maxqaunt = cell(Nproteins,1);
    Unique_indentifer_maxqaunt{1} = {'Unique_identifier'}; % replicate_ProteinName
    for jj = 1:Nproteins
      Replicate_number = num_val{1}(jj,1);
      %convert to number from string
      Unique_indentifer_maxqaunt{jj,1}=strcat(mat2str(Replicate_number),'_',txt_val{ii}{jj,1});
    end
    txt_val{ii} = [Unique_indentifer_maxqaunt, txt_val{ii}];
  end
  clear tmp2
  
  % Import Gaussian data
  GaussData = cell(Nchannels,1);
  Ngauss = zeros(Nchannels,1);
  Nreplicates = zeros(Nchannels,1);
  for ii = 1:Nchannels
    tmp1 = fopen(GaussIn{ii}); %This corresponds to the output from the Gaus script
    tmp2 = textscan(tmp1, '%s', 'Delimiter',',');
    %Reshape data for use in analysis
    tmp3=size(tmp2{:});
    GaussData{ii}=reshape(tmp2{:},10,(tmp3(1)/10))';
    
    % Summary number: Gaussians
    Ngauss(ii) = length(GaussData{ii})-1; %minus one to remove header, corresponds to the gaussian fitted
    
    % Add unique identifiers to the Gaussian lists
    Unique_indentifer=cell(1,1);
    Unique_indentifer{1,1} = ''; % replicate_GaussianNumber_channel_ProteinName
    Unique_indentifer{1,2} = 'Replicate_Protein_identifier'; % replicate_ProteinName
    for jj = 2:size(GaussData{ii},1)
      Replicate_number=str2double(GaussData{ii}{jj,3});
      %convert to number from string
      Index_Gaussian_number=str2double(GaussData{ii}{jj,1});
      Unique_indentifer{jj,1}=[num2str(Replicate_number) '_' num2str(Index_Gaussian_number) '_' Experimental_channels{ii} '_' GaussData{ii}{jj,4}];
      Unique_indentifer{jj,2}=[num2str(Replicate_number) '_' GaussData{ii}{jj,4}];
    end
    GaussData{ii} = [Unique_indentifer, GaussData{ii}]; %Add unique protein_replicate identifier
  end
  
  %expected amount of proteins, note 1 is the height
  %Standard_area=Nfraction*1*(1/Diltuion_factor_master_mix);
  
  
  % Make GaussSummary structure
  % This puts Gaussians for the same protein in one line  
  clear GaussSummary
  
  %Varible related to gaussians
  All_guassians_identified_redunent_list = [];
  for ii = 1:Nchannels
    tmp = unique(GaussData{ii}(2:Ngauss(ii)+1,2));
    All_guassians_identified_redunent_list = [All_guassians_identified_redunent_list;tmp];
  end
  All_guassians_identified_nonredunent_list=unique(All_guassians_identified_redunent_list);
  NuniqueGauss=length(All_guassians_identified_nonredunent_list);
  
  % pre-allocate?
  GaussSummary(2).Center(1) = -1;
  for ii = 1:Nchannels
    GaussSummary(ii).Unique_identifier = cell(NuniqueGauss,5);
    GaussSummary(ii).Guassian_index_number = zeros(NuniqueGauss,5);
    GaussSummary(ii).Protein_number = zeros(NuniqueGauss,5);
    GaussSummary(ii).Replicate = zeros(NuniqueGauss,5);
    GaussSummary(ii).Protein = cell(NuniqueGauss,5);
    GaussSummary(ii).Height = zeros(NuniqueGauss,5);
    GaussSummary(ii).Center = zeros(NuniqueGauss,5);
    GaussSummary(ii).Width = zeros(NuniqueGauss,5);
    GaussSummary(ii).SSE = zeros(NuniqueGauss,5);
    GaussSummary(ii).adjrsquare = zeros(NuniqueGauss,5);
    GaussSummary(ii).Complex_size = zeros(NuniqueGauss,5);
  end
  
  Number_guassian_pre_protein=zeros(Nchannels,NuniqueGauss);
  Location_of_guassians=zeros(Nchannels,NuniqueGauss,5);
  % Create array of values grouped by Protein and replicate
  % This loop goes through all unique replicate-protein combinations.
  % e.g. 1_A0AVT1, 1_A3KN83, 1_A8MV29, ...
  for ii = 1:Nchannels
    for gg = 1:NuniqueGauss
      
      %define protein to look for
      Protein_name = All_guassians_identified_nonredunent_list(gg);
      %find the position in MvsL
      internal_location_guassians = find(strcmp(Protein_name, GaussData{ii}(:,2)));
      number_counter1=length(internal_location_guassians);
      gcount=0;
      for kk = 1:number_counter1
        if str2double(GaussData{ii}(internal_location_guassians(kk),8))>fraction_to_plot
        else
          gcount = gcount+1;
          GaussSummary(ii).Unique_identifier{gg,gcount} = GaussData{ii}(internal_location_guassians(kk),1);
          GaussSummary(ii).Guassian_index_number(gg,gcount) = str2double(GaussData{ii}{internal_location_guassians(kk),3});
          GaussSummary(ii).Protein_number(gg,gcount) = str2double(GaussData{ii}{internal_location_guassians(kk),4});
          GaussSummary(ii).Replicate(gg,gcount) = str2double(GaussData{ii}{internal_location_guassians(kk),5});
          GaussSummary(ii).Protein{gg,gcount} = GaussData{ii}{internal_location_guassians(kk),6};
          GaussSummary(ii).Height(gg,gcount) = str2double(GaussData{ii}{internal_location_guassians(kk),7});
          GaussSummary(ii).Center(gg,gcount) = str2double(GaussData{ii}{internal_location_guassians(kk),8});
          GaussSummary(ii).Width(gg,gcount) = str2double(GaussData{ii}{internal_location_guassians(kk),9});
          GaussSummary(ii).SSE(gg,gcount) = str2double(GaussData{ii}{internal_location_guassians(kk),10});
          GaussSummary(ii).adjrsquare(gg,gcount) = str2double(GaussData{ii}{internal_location_guassians(kk),11});
          GaussSummary(ii).Complex_size(gg,gcount) = str2double(GaussData{ii}{internal_location_guassians(kk),12});
          
          %Calculate Area of Gaussian
          guassian_area_fun =  @(x)( GaussSummary(ii).Height(gg,gcount)...
            *exp(-((x- GaussSummary(ii).Center(gg,gcount))...
            /GaussSummary(ii).Width(gg,gcount)).^2));
          guassian_area= integral(guassian_area_fun,1,Nfraction);
          
          %Save Gaussian area for MvsL
          GaussSummary(ii).Gaussian_area(gg,gcount)= guassian_area;
        end
        GaussSummary(ii).Protein_name(gg) = GaussData{ii}(internal_location_guassians(1),6);
        GaussSummary(ii).Replicate_Protein_identifier(gg) = GaussData{ii}(internal_location_guassians(1),2);
      end
    end
    GaussSummary(ii).Protein_name(cellfun('isempty',GaussSummary(ii).Protein_name)) = {''};
  end
  
  %Define the size of the MvsL and HvsL arrays
  Dimension_of_array=size(GaussSummary(1).Center);
  
    
  % Determine which Gaussians are unique to each channel
  
  clear Combined_Gaussians
  
  %Determine which gaussians are unique and which are shared between MvsL and
  %HvsL of the same replicate
  gausscount = 0;
  for gg = 1:NuniqueGauss
    
    % 1. Pick the channel with the best R^2
    R2 = zeros(Nchannels,1);
    for ii = 1:Nchannels
      R2(ii) = GaussSummary(ii).adjrsquare(gg,1);
    end
    [~, Ibest] = max(R2); % the channel with the best R^2
    Inotbest = find(~ismember(1:Nchannels,Ibest)); % the other channels
    Cbest = GaussSummary(Ibest).Center(gg,:);
    Cbest(isnan(Cbest) | Cbest==0) = [];
    
    % 2. For each not-best Gaussian, check if it's within User_Window of a best Gaussian
    Gausskeep_notbest = cell(length(Inotbest),1);
    Cbest_multireplicates = zeros(length(Cbest),length(Experimental_channels));
    Cbest_multireplicates(:,Ibest) = 1;
    for ii = 1:length(Inotbest)
      % Get Gaussian Centers in this not-best channel
      tmp = GaussSummary(Inotbest(ii)).Center(gg,:);
      tmp = tmp(tmp>0 & ~isnan(tmp));
      
      Gausskeep_notbest{ii} = zeros(size(tmp));
      for jj = 1:length(tmp)
        D = find(abs(Cbest - tmp(jj)) < User_Window);
        
        if isempty(D) % Gaussian is unique, so keep it
          Gausskeep_notbest{ii}(jj) = 1;
        else % Gaussian is not unique, so 'add it' to the best Gaussian
          Cbest_multireplicates(D,Inotbest(ii)) = 1;
        end
      end
    end
    
    % 3. Return all Gaussians in the best channel.
    for ii = 1:length(Cbest)
      gausscount = gausscount+1;
      Combined_Gaussians.Unique_identifier(gausscount) = GaussSummary(Ibest).Replicate_Protein_identifier(gg);
      Combined_Gaussians.Protein_name{gausscount} = GaussSummary(Ibest).Protein{gg,ii};
      Combined_Gaussians.Replicate(gausscount) = GaussSummary(Ibest).Replicate(gg,ii);
      Combined_Gaussians.Center(gausscount) = GaussSummary(Ibest).Center(gg,ii);
      Combined_Gaussians.Height(gausscount) = GaussSummary(Ibest).Height(gg,ii);
      Combined_Gaussians.Width(gausscount) = GaussSummary(Ibest).Width(gg,ii);
      Combined_Gaussians.adjrsquare(gausscount) = GaussSummary(Ibest).adjrsquare(gg,ii);
      Combined_Gaussians.Channels{gausscount} = strjoin(Experimental_channels(Cbest_multireplicates(ii,:)==1), '+');
      Combined_Gaussians.SSE(gausscount) = GaussSummary(Ibest).SSE(gg,ii);
      Combined_Gaussians.Complex_size(gausscount) = GaussSummary(Ibest).Complex_size(gg,ii);
    end
    
    % 4. Return all Gaussians in other channels that have NOT been matched.
    for jj = 1:length(Inotbest)
      g2add = find(Gausskeep_notbest{jj});
      for kk = 1:length(g2add)
        gausscount = gausscount+1;
        Combined_Gaussians.Unique_identifier(gausscount) = GaussSummary(Inotbest(jj)).Replicate_Protein_identifier(gg);
        Combined_Gaussians.Protein_name{gausscount} = GaussSummary(Inotbest(jj)).Protein{gg,g2add(kk)};
        Combined_Gaussians.Replicate(gausscount) = GaussSummary(Inotbest(jj)).Replicate(gg,g2add(kk));
        Combined_Gaussians.Center(gausscount) = GaussSummary(Inotbest(jj)).Center(gg,g2add(kk));
        Combined_Gaussians.Height(gausscount) = GaussSummary(Inotbest(jj)).Height(gg,g2add(kk));
        Combined_Gaussians.Width(gausscount) = GaussSummary(Inotbest(jj)).Width(gg,g2add(kk));
        Combined_Gaussians.adjrsquare(gausscount) = GaussSummary(Inotbest(jj)).adjrsquare(gg,g2add(kk));
        Combined_Gaussians.Channels{gausscount} = strjoin(Experimental_channels(Inotbest(jj)), '+');
        Combined_Gaussians.SSE(gausscount) = GaussSummary(Inotbest(jj)).SSE(gg,g2add(kk));
        Combined_Gaussians.Complex_size(gausscount) = GaussSummary(Inotbest(jj)).Complex_size(gg,g2add(kk));
      end
    end
  end
  
  % Total number of Gaussians
  number_of_proteins = gausscount;
  Unique_gaussians_in_eachchannel = zeros(Nchannels,1);
  for ii = 1:Nchannels
    Unique_gaussians_in_eachchannel(ii) = sum(~cellfun(@isempty,regexp(Combined_Gaussians.Channels(:),Experimental_channels{ii})));
  end
  sjoin = strjoin(Experimental_channels(:), '+');
  shared_gaussian_counter = sum(strcmp(Combined_Gaussians.Channels(:),sjoin));
  
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 2. Detect 2-fold changes
  % This is done with Gaussian Center
  tic
  fprintf('    2. Detect 2-fold changes (WITHIN REPLICATE)')
  
  %%%%%%% you can insert a comparisonpairs loop here
  
  % for this comparison pair
  % find which of the user.silacratios is the no-treatment, i.e. denominator in the fold change, and
  % which is the treatment, i.e. numerator.
%   SR2use(1) = find(ismember(user.silacratios,user.comparisonpairs{1}));
%   SR2use(2) = find(ismember(user.silacratios,user.comparisonpairs{2}));
%   I1 = ismember(user.silacratios(SR2use),user.treatmentcondition); % treatment, numerator
%   I2 = SR2use(~I1); % no-treatment, denominator
%   I1 = find(I1);
  I1 = find(ismember(user.silacratios,user.comparisonpairs{1})); % treatment, numerator
  I2 = find(ismember(user.silacratios,user.comparisonpairs{2})); % non-treatment, denominator
  
  Combined_Gaussians.log2_of_gaussians = nan(number_of_proteins,1);
  for gg = 1:number_of_proteins
    % Record which channels are the denominator/numerator in this fold change
    Combined_Gaussians.numeratorChannnel{gg} = user.silacratios{I1};
    Combined_Gaussians.denominatorChannnel{gg} = user.silacratios{I2};
    
    % location of Protein in raw textdata
    Idata = find(strcmp(Combined_Gaussians.Unique_identifier{gg},txt_val{1}(:,1)));
    Center_to_test = round(Combined_Gaussians.Center(gg));
    
    % If Center is out of range, abort
    %if Center_to_test-2+frac1-1<1 || Center_to_test+2+frac1-1>Nfraction+1
    if Center_to_test-2<1 || Center_to_test+2>frac2-frac1+1
      continue;
    end
    
    % retrieve raw data
    rawratio = nan(Nchannels,5);
    baddata = 0;
    for ii = 1:Nchannels
      I = Center_to_test+frac1-1 -2:Center_to_test+frac1-1 +2;
      rawratio(ii,:) = (num_val{ii}(Idata,I));
      baddata = baddata | nnz(rawratio(ii,:))==0;
    end
    
    % Check if the fold change is measurable
    a = rawratio(:,2:4);
    baddata = baddata | sum(isnan(a(:)))/Nchannels>=2;
    if baddata;continue;end
    
    % If data is good, get fold changes between the comparison pair
    Combined_Gaussians.log2_of_gaussians(gg) = log2(nanmean(rawratio(I1,:))/nanmean(rawratio(I2,:)));
    

  end
  Combined_Gaussians.log2_normalised_gaussians = Combined_Gaussians.log2_of_gaussians - nanmean(Combined_Gaussians.log2_of_gaussians);
  
  Ibad = isnan(Combined_Gaussians.log2_normalised_gaussians);
  Iinc = Combined_Gaussians.log2_normalised_gaussians>1;
  Idec = Combined_Gaussians.log2_normalised_gaussians<-1;
  Combined_Gaussians.Observed_change = cell(number_of_proteins,1);
  Combined_Gaussians.Observed_change(Iinc & ~Ibad) = {'Increase'};
  Combined_Gaussians.Observed_change(Idec & ~Ibad) = {'Decrease'};
  Combined_Gaussians.Observed_change(~Idec & ~Iinc & ~Ibad) = {'No change'};
  Combined_Gaussians.Observed_change(Ibad) = {'Unquantifiable'};
  
  
  % Determine the trends of guassian curves
  % Aim: look at the determined gaussians within single replicate and compare
  % if the all the gaussians change, some change or none of them change
  
  %Define Variables
  No_change_consistent_across_gaus=0;
  Increase_consistent_across_gaus=0;
  Decrease_consistent_across_gaus=0;
  Increase_inconsistent_across_gaus=0;
  Decrease_inconsistent_across_gaus=0;
  Increase_and_decrease_across_gaus=0;
  
  %Loop over replicate-protein pairs to assess guassian information
  trendString = cell(NuniqueGauss,2);
  Ngaussperrep_prot = zeros(NuniqueGauss,1);
  for gg = 1:NuniqueGauss
    rep_prot = GaussSummary(1).Replicate_Protein_identifier{gg};
    
    % Find this replicate-protein pair in Combined_Gaussians
    I = find(strcmp(rep_prot,Combined_Gaussians.Unique_identifier));
    Ngaussperrep_prot(gg) = length(I);
    
    % Find the increase/decrease for this replicate-protein's Gaussians
    Ninc = sum(strcmp('Increase',Combined_Gaussians.Observed_change(I)));
    Ndec = sum(strcmp('Decrease',Combined_Gaussians.Observed_change(I)));
    Nnc = sum(strcmp('No change',Combined_Gaussians.Observed_change(I)));
    
    if Nnc>0 && Ninc== 0 && Ndec==0
      trendString(gg,:) = {'No change','Consistent across gaussians'};
      No_change_consistent_across_gaus=1+ No_change_consistent_across_gaus;
    elseif Nnc==0 && Ninc>0 && Ndec==0
      trendString(gg,:) = {'Increase','Consistent across gaussians'};
      Increase_consistent_across_gaus=1+ Increase_consistent_across_gaus;
    elseif Nnc==0 && Ninc==0 && Ndec>0
      trendString(gg,:) = {'Decrease','Consistent across gaussians'};
      Decrease_consistent_across_gaus=1+ Decrease_consistent_across_gaus;
    elseif Ninc>0 && Nnc>0
      trendString(gg,:) = {'Increase','Inconsistent across gaussians'};
      Increase_inconsistent_across_gaus=1+ Increase_inconsistent_across_gaus;
    elseif Ndec>0 && Nnc>0
      trendString(gg,:) = {'Decrease','Inconsistent across gaussians'};
      Decrease_inconsistent_across_gaus=1+ Decrease_inconsistent_across_gaus;
    elseif Ndec>0 && Ninc>0
      trendString(gg,:) = {'+/-','Inconsistent across gaussians'};
      Increase_and_decrease_across_gaus=1+Increase_and_decrease_across_gaus;
    else
      trendString(gg,:) = {'NaN','Unquantifiable'};
    end
  end
  
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 3. Create Gaussian Master List, combinining replicates and channels.
  % Comparsion 3 Part A:  Determine which Gaussians are shared across replicates
  % Aim: determined gaussians which are shared
  tic
  fprintf('    3. Create Gaussian Master List, combinining replicates.')
  
  % i) loop through unique protein names
  % ii) find all fitted Gaussians for that protein
  % iii) find the best replicate-channel for that protein using R^2
  % iv) merge all remining Gaussians with that best list
  % v) report the best list
  % vi) report the unmerged proteins
  
  %Unique gene name
  All_protein_names = [];
  for ii = 1:Nchannels
    All_protein_names=[All_protein_names GaussSummary(ii).Protein_name];
  end
  Ibad = cellfun('isempty',All_protein_names);
  All_protein_names = All_protein_names(~Ibad);
  Unique_protein_names=unique(All_protein_names);
  %remove NaN from list
  test_nan=ismember(Unique_protein_names, {'NaN'});
  Unique_protein_names(test_nan==1)=[];
  %number of unique names in list
  
  % i) loop through unique protein names
  countMaster = 0;
  clear Master_Gaussian_list
  for ii = 1:length(Unique_protein_names)
    protName = Unique_protein_names{ii};
    
    % ii) find all fitted Gaussians for that protein
    I = cell(Nchannels,1);
    for jj = 1:Nchannels
      %tmp = strfind(GaussSummary(jj).Protein_name,protName);
      %I{jj} = find(~cellfun('isempty',tmp));
      I{jj} = find(ismember(GaussSummary(jj).Protein_name, protName));
    end
    
    % iii) find the best replicate-channel for this protein using R^2
    count = 0;
    Gauss_thisprotein = zeros(Nchannels*replicate_num*5,5);
    for jj = 1:Nchannels
      for kk = 1:length(I{jj}) % replicates
        gausstmp = GaussSummary(jj).Center(I{jj}(kk),:);
        gausstmp = gausstmp(gausstmp>0 & ~isnan(gausstmp));
        for gg = 1:length(gausstmp)
          count = count+1;
          Gauss_thisprotein(count,1) = GaussSummary(jj).Center(I{jj}(kk),gg); % center
          Gauss_thisprotein(count,2) = jj; % channel
          Gauss_thisprotein(count,3) = I{jj}(kk); % replicate index
          Gauss_thisprotein(count,4) = gg; % gaussian number index
          Gauss_thisprotein(count,5) = GaussSummary(jj).adjrsquare(I{jj}(kk),gg);
        end
      end
    end
    Gauss_thisprotein = Gauss_thisprotein(1:count,:);
    if isempty(Gauss_thisprotein);continue;end
    [~,I2] = max(Gauss_thisprotein(:,5));
    bestChannel = Gauss_thisprotein(I2,2);
    bestReplicate = Gauss_thisprotein(I2,3);
    I3 = Gauss_thisprotein(:,2)==bestChannel & Gauss_thisprotein(:,3)==bestReplicate;
    bestCenters = Gauss_thisprotein(I3,1);
    
    % iv) merge Gaussians with the Gaussians from the best replicate-channel as well as new Gaussians
    Gauss2keep = false(size(Gauss_thisprotein,1),1);
    for jj = 1:length(Gauss2keep)
      Gauss2keep(jj) = (sum(abs(Gauss_thisprotein(jj,1)-bestCenters)<User_Window)==0 | I3(jj)) &  Gauss_thisprotein(jj,1)>1;
      if Gauss2keep(jj)
        bestCenters = [bestCenters;Gauss_thisprotein(jj,1)];
      end
    end
    Gauss_thisprotein = Gauss_thisprotein(Gauss2keep,:);
    
    % v) return Gaussians
    % note: need channel, replicate, gaussian-number-for-this-protein
    for jj = 1:size(Gauss_thisprotein,1)
      countMaster = countMaster+1;
      chan = Gauss_thisprotein(jj,2);
      repi = Gauss_thisprotein(jj,3);
      ggi = Gauss_thisprotein(jj,4);
      
      Master_Gaussian_list.Protein_name{ii,jj} = Unique_protein_names{ii};
      Master_Gaussian_list.Unique_identifier(ii,jj) = GaussSummary(chan).Unique_identifier{repi};
      Master_Gaussian_list.Protein_number(ii,jj) = GaussSummary(chan).Protein_number(repi,ggi);
      Master_Gaussian_list.Replicate(ii,jj) = GaussSummary(chan).Replicate(repi,ggi);
      Master_Gaussian_list.Channel{ii,jj} = Experimental_channels{chan};
      Master_Gaussian_list.Guassian_index_number(ii,jj) = GaussSummary(chan).Guassian_index_number(repi,ggi);
      Master_Gaussian_list.Center(ii,jj) = GaussSummary(chan).Center(repi,ggi);
      Master_Gaussian_list.Height(ii,jj) = GaussSummary(chan).Height(repi,ggi);
      Master_Gaussian_list.Width(ii,jj) = GaussSummary(chan).Width(repi,ggi);
      Master_Gaussian_list.SSE(ii,jj) = GaussSummary(chan).SSE(repi,ggi);
      Master_Gaussian_list.adjrsquare(ii,jj) = GaussSummary(chan).adjrsquare(repi,ggi);
      Master_Gaussian_list.Complex_size(ii,jj) = GaussSummary(chan).Complex_size(repi,ggi);
      Master_Gaussian_list.Replicate_Protein_identifier{ii,jj} = GaussSummary(chan).Replicate_Protein_identifier{repi};
      Master_Gaussian_list.Gaussian_area(ii,jj) = GaussSummary(chan).Gaussian_area(repi,ggi);
    end
  end
    
  %Create Chromatogram file containing the raw isotoplogue data for each guassian that has been kept
  Chromatograms=zeros(countMaster,Nfraction);
  countMaster = 0;
  for ii = 1:length(Unique_protein_names)
    for jj = 1:nnz(Master_Gaussian_list.Center(ii,:))
      countMaster = countMaster+1;
      
      chan = find(strcmp(Master_Gaussian_list.Channel{ii,jj},Experimental_channels));
      I = find(strcmp(Master_Gaussian_list.Replicate_Protein_identifier{ii,jj},txt_val{chan}(:,1)));
      if length(I)~=1
        error('Comparison: ummm!! oh dear')
      end
      Chromatograms(countMaster,:) = num_val{chan}(I,2:end);
    end
  end
  
  Finalised_Master_Gaussian_list = Master_Gaussian_list;
  Total_number_of_unique_gaussians_master_list = countMaster;
  Finalised_Master_Gaussian_list.Total_number_of_unique_gaussians_master_list = Total_number_of_unique_gaussians_master_list;
  
  
  
  % Use fold-change > 2 to determine if complexes change
  
  %%%%%%% you can insert a comparisonpairs loop here
  
  % for this comparison pair
  % find which of the user.silacratios is the no-treatment, i.e. denominator in the fold change, and
  % which is the treatment, i.e. numerator.
  I1 = find(ismember(user.silacratios,user.comparisonpairs{1})); % treatment, numerator
  I2 = find(ismember(user.silacratios,user.comparisonpairs{2})); % non-treatment, denominator
  
  foldLabel = cell(length(Unique_protein_names),30);
  foldChange = nan(length(Unique_protein_names),30);
  foldChange_byreplicate = nan(length(Unique_protein_names),10,replicate_num);
  foldLabel_byreplicate = cell(length(Unique_protein_names),10,replicate_num);
  for ii = 1:length(Unique_protein_names)
    protName = Finalised_Master_Gaussian_list.Protein_name{ii,1};
    if isempty(protName)
      continue
    end
    tmp = strfind(Combined_Gaussians.Protein_name,protName);
    for hh = 1:nnz(Finalised_Master_Gaussian_list.Center(ii,:))
      for rep = 1:replicate_num
        %rep = Finalised_Master_Gaussian_list.Replicate(ii,hh);
        
        %locate protein to compare
        Iraw = find(strcmp(protName,txt_val{1}(:,2)) & num_val{1}(:,1)==rep);
        rounded_center = floor(Finalised_Master_Gaussian_list.Center(ii,hh)); % determine the center of the Gaussian
        
        %Reset rounded_center if Gaussian found in replicate
        %Check if a Gaussian peak was detected within two fractions of this master Gaussian value?
        %Find protein in Combined_Gaussian analysis
        
        location_Protein_in_combined_Gaus = find(~cellfun('isempty',tmp));
        for find_rep_counter=1:length(location_Protein_in_combined_Gaus) %As Combined Gaussian is used to create Final will also have atleast one value
          if rep==Combined_Gaussians.Replicate(location_Protein_in_combined_Gaus(find_rep_counter))
            %check if Center is within 2 fractions of master list center
            temp_value=abs(Combined_Gaussians.Center(location_Protein_in_combined_Gaus(find_rep_counter))- rounded_center);
            if  temp_value<=(User_Window+0.49)
              rounded_center= floor(Combined_Gaussians.Center(location_Protein_in_combined_Gaus(find_rep_counter)));
            end
          end
        end
        
        % Is Center out of bounds?
        if rounded_center-2+frac1-1<2+1 || rounded_center+2+frac1-1>Nfraction+1
          foldLabel{ii,hh} = 'Unquantifiable_Not_observed_in_replicate';
          continue;
        end
        
        % retrieve raw data
        rawratio = nan(Nchannels,5);
        baddata = false(Nchannels,1);
        for jj = 1:Nchannels
          I = rounded_center-2+frac1-1:rounded_center+2+frac1-1;
          rawratio(jj,:) = num_val{jj}(Iraw,I);
          baddata(jj) = sum(isnan(rawratio(jj,:)))==5 | sum(isnan(rawratio(jj,2:4)))>=2;
        end
        
        % Is the raw data around the Center bad?
        for jj = 1:Nchannels
          if baddata(jj)
            foldLabel{ii,hh} = ['Unquantifiable_Not_observed_in_' Experimental_channels{jj} ' channel'];
            continue;
          end
        end
        
        % Calculate log2-fold-change
        foldChange(ii,hh) = log2(nanmean(rawratio(I1,:)./rawratio(I2,:)));
        foldChange_byreplicate(ii,hh,rep) = foldChange(ii,hh);
        if foldChange(ii,hh)>2
          foldLabel{ii,hh} = 'Increase';
        elseif foldChange(ii,hh)<0.5
          foldLabel{ii,hh} = 'Decrease';
        else
          foldLabel{ii,hh} = 'No change';
        end
        foldLabel_byreplicate{ii,hh,rep} = foldLabel{ii,hh};
        
      end
    end
  end
  
  % Determine whether changes were consistent across replicates
  Decrease_protein_persaus_input = cell(size(Finalised_Master_Gaussian_list.Center));
  Increase_protein_persaus_input = cell(size(Finalised_Master_Gaussian_list.Center));
  Increase_and_Decrease_protein_persaus_input = cell(size(Finalised_Master_Gaussian_list.Center));
  Decrease_gauss_persaus_input = cell(size(Finalised_Master_Gaussian_list.Center));
  Increase_gauss_persaus_input = cell(size(Finalised_Master_Gaussian_list.Center));
  Increase_and_Decrease_gauss_persaus_input = cell(size(Finalised_Master_Gaussian_list.Center));
  
  for ii= 1:length(Unique_protein_names)
    for hh= 1:nnz(Finalised_Master_Gaussian_list.Center(ii,:))
      incyn = sum(foldChange_byreplicate(ii,hh,:)>1)>=2; % do two more replicates increase?
      decyn = sum(foldChange_byreplicate(ii,hh,:)<-1)>=2; % do two more replicates decrease
      
      if incyn && ~decyn
        Increase_protein_persaus_input{ii} = '+';
        Increase_gauss_persaus_input{ii} = '+';
      elseif ~incyn && decyn
        Decrease_protein_persaus_input{ii} = '+';
        Decrease_gauss_persaus_input{ii} = '+';
      elseif incyn && decyn
        Increase_and_Decrease_protein_persaus_input{ii} = '+';
        Increase_and_Decrease_gauss_persaus_input{ii} = '+';
      end
    end
  end
  
  Finalised_Master_Gaussian_list.foldLabel = foldLabel;
  Finalised_Master_Gaussian_list.foldChange = foldChange;
  Finalised_Master_Gaussian_list.foldChange_normalized = foldChange - nanmedian(foldChange(:));
  Finalised_Master_Gaussian_list.foldChange_byreplicate = foldChange_byreplicate;
  Finalised_Master_Gaussian_list.foldLabel_byreplicate = foldLabel_byreplicate;
  Finalised_Master_Gaussian_list.Increase_and_Decrease_protein_persaus_input = Increase_and_Decrease_protein_persaus_input;
  Finalised_Master_Gaussian_list.Increase_protein_persaus_input = Increase_protein_persaus_input;
  Finalised_Master_Gaussian_list.Decrease_protein_persaus_input = Decrease_protein_persaus_input;
  Finalised_Master_Gaussian_list.Increase_and_Decrease_gauss_persaus_input = Increase_and_Decrease_gauss_persaus_input;
  Finalised_Master_Gaussian_list.Increase_gauss_persaus_input = Increase_gauss_persaus_input;
  Finalised_Master_Gaussian_list.Decrease_gauss_persaus_input = Decrease_gauss_persaus_input;
  
  %%%%%%% you can end the comparisonpairs loop here
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  
  %% 4. Use mwwtest and ttest to determine if complexes change
  tic
  fprintf('    4. MWW Test and T-test across replicates')
  
  countMaster = 0;
  p_tt = nan(length(Unique_protein_names),30);
  df_tt = nan(length(Unique_protein_names),30);
  p_mww = nan(length(Unique_protein_names),30);
  u_mww = nan(length(Unique_protein_names),30);
  nobs = nan(length(Unique_protein_names),30);
  for ii = 1:length(Unique_protein_names) % Protein
    protName = Finalised_Master_Gaussian_list.Protein_name{ii};
    
    % Find this protein in the raw data
    location_Protein_in_raw_textdata = find(strcmp(protName,txt_val{1}(:,2)));
    
    for jj = 1:nnz(Finalised_Master_Gaussian_list.Center(ii,:)) % Gaussian
      
      % The Center of this Gaussian
      Center_to_test = round(Finalised_Master_Gaussian_list.Center(ii,jj)) + 6;
      
      % Is Center out of bounds?
      if Center_to_test<1 || Center_to_test>Nfraction
        foldLabel{ii,hh} = 'Unquantifiable_Not_observed_in_replicate';
        continue;
      end
      
      % Get data from all channels and replicates
      data = cell(Nchannels,1);
      nobs(ii,jj) = 10^6; % keep track of how many data points are in each channel
      Ibad = zeros(length(location_Protein_in_raw_textdata),1);
      for kk = 1:Nchannels % Channel
          tmp = log(num_val{kk}(location_Protein_in_raw_textdata,Center_to_test));
          data{kk} = tmp;
          Ibad = Ibad | isnan(data{kk});
          
      end
      for kk = 1:Nchannels
          data{kk} = data{kk}(not(Ibad));
          nobs(ii,jj) = min(nobs(ii,jj),length(data{kk}));
      end
      
      % Does the data pass Nick's good-data checks?
      % Check 1: Is channel 1 data all nans?
      % Check 2: Is channel 2 data all nans?
      % Check 3: Are there more than 2 nans around the center in channel 1?
      % Check 4: Are there more than 2 nans around the center in channel 1?
      % NB this ^ is irrelevant, since nans are just removed.
      if nobs(ii,jj)<2;continue;end
      
      % MWWTEST
      stats = mwwtest(data{1},data{2});
      p_mww(ii,jj) = stats.p;
      u_mww(ii,jj) = stats.U;
      
      % TTEST
      [~,p_tt(ii,jj),~,stats] = ttest(data{1},data{2});
      df_tt(ii,jj) = stats.df;
      
    end
  end
  
  %Add ttest output
  Finalised_Master_Gaussian_list.Ttest_p_values = p_tt;
  Finalised_Master_Gaussian_list.degree_freedom = df_tt;
  Finalised_Master_Gaussian_list.Number_observation = nobs;
  Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue = p_tt < 0.05/Total_number_of_unique_gaussians_master_list;
  
  %add MWW test output
  Finalised_Master_Gaussian_list.MWW_p_values = p_mww;
  Finalised_Master_Gaussian_list.MWW_U = u_mww;
  Finalised_Master_Gaussian_list.MWW_Satifies_Adjusted_pvalue = p_mww < 0.05/Total_number_of_unique_gaussians_master_list;
  
  Singficant_change_ttest_persaus_protein = cell(length(Unique_protein_names),1);
  Singficant_change_mww_persaus_protein = cell(length(Unique_protein_names),1);
  for ii= 1:length(Unique_protein_names)
    for jj= 1:nnz(Finalised_Master_Gaussian_list.foldChange(ii,:))
      %Check MWW test
      if Finalised_Master_Gaussian_list.MWW_Satifies_Adjusted_pvalue(ii,jj)
        Singficant_change_mww_persaus_protein(ii)= {'+'};
      end
      %Check ttest
      if Finalised_Master_Gaussian_list.Satifies_Adjusted_pvalue
        Singficant_change_ttest_persaus_protein(ii)= {'+'};
      end
    end
  end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 5. Write output
  tic
  fprintf('    5. Write output\n')
  
  writeOutput_foldchanges
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 6. Make figures
  tic
  fprintf('    6. Make figures\n')
  
  try
    makeFigures_foldchanges
  catch
    warning('makeFigures_foldchanges failed.')
  end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
end

diary('off')

