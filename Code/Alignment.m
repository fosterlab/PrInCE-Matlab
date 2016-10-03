%%%%%%%%%%%%%%% Logic:
% Since replicates can be run on different columns, etc., it's necessary to align replicates. Each
% replicate will be offset by a constant amount, so it's sufficient to find offset between a few
% landmarks
% i. Read chromatogram data (MaxQuant verbose output) and Gaussian fitting data.
% ii. Impute missing chromatogram values.
% iii. Find all proteins with a single Gaussian in all replicates.
% iv. Choose the best replicate, i.e. the one with the most overlap of single-Gaussian proteins with
%     the other two. Align to this replicate.
% v. For every other replicate, fit the line: y_mu = m * x_mu + b, where y_mu is the centers of the
%     best replicate, x_mu is the centers of the other replicate.
% vi. Scale the other replicates, x, by the above function. Use linear interpolation.
%
% NB1 - Best replicate doesn't seem to be chosen using SSE of fits.
% NB2 - x and y likely need to be rescaled to reflect fraction spacing (e.g. 20 in 10 minutes, etc.)
%     - That is, the spacing b/w fractions is not constant. That seems major, right?


%%%%%%%%%%%%%%% Big changes:
% 1. Impute multiple missing values (more than just one, as is currently being done)
% 2. Don't write csv files (these seem useful for debugging, but otherwise aren't needed).
% 3. Include an option not to plot (separate plotting script?).
% 4. Why does Gauss_Build clean the chromatograms more thoroughly? Here we just impute single
%    values. Should we do more?




diary([user.maindir 'logfile.txt'])
disp('Alignment.m')

skipflag = 0;
if user.Nreplicate==1
  disp('* NB: User set Number of Replicates to 1. Skipping Alignment...')
  user.skipalignment = 1;
  skipflag = 1;
end
if user.skipalignment==1 && skipflag==0
  disp('* NB: User set skipalignment to True. Skipping Alignment...')
  skipflag = 1;
end

if ~skipflag
  
  %% 0. Initialize
  tic
  fprintf('\n    0. Initialize')
  
  
  % Load user settings
  maindir = user.maindir;
  Experimental_channels = user.silacratios;
  User_alignment_window1 = user.User_alignment_window1;
  Nreplicates = user.Nreplicate;
  Nchannels = length(Experimental_channels);  % Defines the number of experiments to be compared
  
  
  % Define folders, i.e. define where everything lives.
  codedir = [maindir 'Code/']; % where this script lives
  funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
  datadir = [maindir 'Output/Data/Alignment/']; % where data files live
  figdir = [maindir 'Output/Figures/Alignment/']; % where figures live
  % Make folders if necessary
  if ~exist(datadir, 'dir'); mkdir(datadir); end
  if ~exist([maindir '/Output'], 'dir'); mkdir([maindir '/Output']); end
  if ~exist([maindir '/Output/Data'], 'dir'); mkdir([maindir '/Output/Data']); end
  if ~exist([maindir '/Output/Figures'], 'dir'); mkdir([maindir '/Output/Figures']); end
  if ~exist([maindir '/Output/tmp'], 'dir'); error('Data from Gauss_Build not found'); end
  if ~exist([maindir '/Output/Data/GaussBuild'], 'dir'); error('Data from Gauss_Build not found'); end
  if ~exist(datadir, 'dir'); mkdir(datadir); end
  if ~exist(figdir, 'dir'); mkdir(figdir); end
 
  
  flag3 = 0;
  
  
  % for now, read files in from Nick's data
  % in the future read them in from my data
  if user.nickflag %nick's output
    tmpdir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/LegacyData/';
    GaussInputFile = cell(Nchannels, Nreplicates);
    GassSumInputFile = cell(Nchannels, Nreplicates);
    for ei=1:Nchannels
      tmp = Experimental_channels{ei};
      for replicates= 1:Nreplicates
        GaussInputFile{ei,replicates} = [tmpdir tmp '_Combined_OutputGaus_rep' num2str(replicates) '.csv'];
        GassSumInputFile{ei,replicates} = [tmpdir tmp '_Summary_Gausians_for_individual_proteins_rep' num2str(replicates) '.csv'];
      end
    end
  else %my output
    dd = dir([maindir '/Output/tmp/' '*Combined_OutputGaus_rep*csv']);
    tmp = length(dd) / length(Experimental_channels);
    
    % Check that all replicate files were made in Gauss_Build
    if Nreplicates ~= tmp
      disp(['Error: Alignment: Number of replicates set to ' num2str(Nreplicates) ', ' num2str(tmp) ' detected'])
    end
    
    GaussInputFile = cell(Nchannels, Nreplicates);
    GassSumInputFile = cell(Nchannels, Nreplicates);
    for ei=1:Nchannels
      for replicates= 1:Nreplicates
        GaussInputFile{ei,replicates} = [maindir '/Output/tmp/' Experimental_channels{ei} '_Combined_OutputGaus_rep' num2str(replicates) '.csv'];
        GassSumInputFile{ei,replicates} = [maindir '/Output/tmp/' Experimental_channels{ei} '_Summary_Gausians_for_individual_proteins_rep' num2str(replicates) '.csv'];
      end
    end
  end
  % *vsL_Combined_OutputGaus.csv: data on all the fitted Gaussians
  % *vsL_Summary_Gausians_for_individual_proteins.csv: how many Gaussians were fitted per protein
  
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 1. Read input
  tic
  fprintf('    1. Read input')
  
  
  % Import MaxQuant data files
  num_val = cell(Nchannels,1);
  txt_val = cell(Nchannels,1);
  for ii = 1:Nchannels
    tmp = importdata(user.MQfiles{ii});
    
    % remove 'sheet1' fields
    if isfield(tmp,'Sheet1')
      tmp = tmp.Sheet1;
    end
    fn = fieldnames(tmp);
    for jj = 1:length(fn)
      tmp1 = tmp.(fn{jj});
      if isfield(tmp1,'Sheet1')
        tmp.(fn{jj}) = tmp.(fn{jj}).Sheet1;
      end
    end
    
    num_val{ii} = tmp.data;
    txt_val{ii} = tmp.textdata(:,1);
    
    % if txt_val includes protein groups, reduce it to the first protein in each group
    for jj = 1:size(txt_val{ii},1)
      tmp = strsplit(txt_val{ii}{jj},';');
      txt_val{ii}{jj} = tmp{1};
    end
  end
  
  
  % Import Gauss fits for each replicate
  %   Gaus_import: mx6, where m is the number of proteins with a fitted Gaussian
  %         columns: Height,Center,Width,SSE,adjrsquare,Complex Size
  Gaus_import = cell(Nchannels, Nreplicates);
  letters = 'abcdefghijklmnopqrstuvwxyz';
  LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
  Nfit = cell(Nchannels,Nreplicates);
  for ci = 1:Nchannels
    for replicates= 1:Nreplicates
      Gaus_import{ci,replicates} = importdata(GaussInputFile{ci,replicates});
      %Summary_gausian_infomration{ci,replicates} = importdata(GassSumInputFile{ci,replicates});
      
      % Ensure Gaus_import.textdata is a single column of protein names
      % simple rule: protein names are the column with the most letters
      a = (cellfun(@(x) ismember(x,LETTERS),Gaus_import{ci,replicates}.textdata(2:end,:),'uniformoutput',0));
      nLETTERS = sum(cell2mat(cellfun(@(x) sum(x),a,'uniformoutput',0)));
      a = (cellfun(@(x) ismember(x,letters),Gaus_import{ci,replicates}.textdata(2:end,:),'uniformoutput',0));
      nletters = sum(cell2mat(cellfun(@(x) sum(x),a,'uniformoutput',0)));
      [~,I] = max(nLETTERS + nletters);
      Gaus_import{ci,replicates}.textdata = Gaus_import{ci,replicates}.textdata(:,I);
      if size(Gaus_import{ci,replicates}.textdata,1) == size(Gaus_import{ci,replicates}.data,1)+1
        Gaus_import{ci,replicates}.header = Gaus_import{ci,replicates}.textdata{1};
        Gaus_import{ci,replicates}.textdata = Gaus_import{ci,replicates}.textdata(2:end);
      end
      
      % Count how many Gaussians were fit to each protein name
      % Needed to identify which proteins have a single Gaussian fit
      Nfit{ci,replicates} = size(Gaus_import{ci,replicates}.textdata,1);
      for ii = 1:size(Gaus_import{ci,replicates}.textdata,1)
        protName = Gaus_import{ci,replicates}.textdata{ii};
        Nfit{ci,replicates}(ii) = sum(ismember(Gaus_import{ci,replicates}.textdata,protName));
      end
      
    end
  end  
  
  % Calibration
  %SEC_fit=polyfit(SEC_size_alignment(1,:),SEC_size_alignment(2,:),1);
  
  %Number of fractions
  [~, fraction_number]=size(num_val{1});
  fraction_number = fraction_number-1;
  if fraction_number~= user.Nfraction
    disp('Alignment: user.Nfraction does not equal detected number of fractions')
  end
  
  % replicates
  replicates =  num_val{1}(:,1);
  
  % Clean chromatograms
  cleandata = cell(Nchannels,1);
  for ii = 1:Nchannels
    cleandata{ii} = nan(size(num_val{ii},1),size(num_val{ii},2)+10);
    %if flag3; cleandata{3} = nan(size(num_val_HvsM,1),size(num_val_HvsM,2)+10);end
    cleandata{ii}(:,1) = num_val{ii}(:,1);
    %if flag3; cleandata{3}(:,3) = num_val_HvsM(:,1);end
    for ri = 1:size(num_val{1},1) % loop over proteins
      cleandata{ii}(ri,2:end) = cleanChromatogram(num_val{ii}(ri,2:end),[1 3 5 7]);
      %if flag3; cleandata{3}(ri,2:end) = cleanChromatogram(num_val_HvsM(ri,2:end),[1 3 5 7]);end
    end
  end
  
  % The data is zero-padded. Find where the real data starts and stops.
  tmp = find(sum(cleandata{ii}==0)==size(cleandata{ii},1));
  frac1 = max([2 tmp(find(tmp<user.Nfraction/2,1,'last'))])+1; % start of real data
  frac2 = tmp(find(tmp>user.Nfraction/2,1,'first'))-1; % end of real data
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 2. Find the best replicates to align to
  % i) Find names of proteins with a single Gaussian
  % ii) Find the overlap in single Gaussians b/w replicates
  % iii) Choose the replicate with maximum overlap
  tic
  fprintf('    2. Find the best replicates to align to')
  
  replicate_to_align_against = nan(Nchannels,1);
  summerised_names_G1 = cell(Nchannels,Nreplicates); % proteins with 1 Gaussian
  for ci = 1:Nchannels
    
    % i) Find names of proteins with a single Gaussian in each replicate
    %summerised_protein_number_G1 = cell(Nreplicates,1);
    for rr= 1:Nreplicates
      %Ngauss = size(Summary_gausian_infomration{ci,rr},1);
      
      Isingle = find(Nfit{ci,rr} == 1 & ~cellfun('isempty',Gaus_import{ci,rr}.textdata)');
      summerised_names_G1{ci,rr} = Gaus_import{ci,rr}.textdata(Isingle,1);
      
      % Using the same syntax as section 3, confirm that these occur exactly once in replicate.
      overlap = find(ismember(Gaus_import{ci,rr}.textdata(:,1),summerised_names_G1{ci,rr}));
      if length(overlap) ~= length(summerised_names_G1{ci,rr})
        error('Alignment: Problem finding proteins with a single Gaussian fit')
      end
    end
    
    % ii) Find the overlap in single Gaussians b/w replicates
    Nintersect = zeros(Nreplicates,Nreplicates);
    for rr1 = 1:Nreplicates
      for rr2 = 1:Nreplicates
        Nintersect(rr1,rr2) = length(intersect(summerised_names_G1{ci,rr1},summerised_names_G1{ci,rr2}));
      end
    end
    Nintersect(eye(Nreplicates)==1) = 0; % remove diagonal, i.e. self-overlap
    Nintersect = sum(Nintersect); % collapse to one dimension
    
    [~,replicate_to_align_against(ci)] = max(Nintersect); % ding ding ding!
    
  end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 3. Calculate best fit lines for adjustment
  % i) Find the overlapping proteins
  % ii) Find the Centers of these proteins
  % iii) Fit a line
  tic
  fprintf('    3. Calculate best fit lines for adjustment')
  
  pfit = nan(Nchannels,Nreplicates,2);
  
  for ci = 1:Nchannels
    align_rep = replicate_to_align_against(ci);
    for rr = 1:Nreplicates
      
      rep_to_align = rr;
      
      % i) Find the overlapping proteins
      overlap = intersect(summerised_names_G1{ci,align_rep},summerised_names_G1{ci,rep_to_align});
      %overlap([1 2]) = [];
      
      % ii) Find their centers
      Ia = find(ismember(Gaus_import{ci,align_rep}.textdata(:,1),overlap));
      Ib = find(ismember(Gaus_import{ci,rep_to_align}.textdata(:,1),overlap));
      Ca = Gaus_import{ci,align_rep}.data(Ia,2); % align to this replicate, x
      Cb = Gaus_import{ci,rep_to_align}.data(Ib,2);% align this replicate, y
      
      % iii) Fit a line
      I = abs(Ca - Cb)<User_alignment_window1;
      %pfit(ci,rr,:) = robustfit(Ca(I), Cb(I));
      pfit(ci,rr,:) = robustfit(Cb(I), Ca(I));
      
    end
  end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 4. Using fitted curves, adjust replicate data
  % i) Gaussian fits: shift Center
  % ii) Chromatograms: shift data points
  tic
  fprintf('    4. Using fitted curves, adjust replicate data')
  
  Adjusted_Gaus_import = Gaus_import;
  
  for ci = 1:Nchannels
    align_rep = replicate_to_align_against(ci);
    for rr = 1:Nreplicates
      
      b = pfit(ci,rr,1); % intercept
      m = pfit(ci,rr,2); % slope
      
      % i) Gaussian fits: shift Center
      adjustedCenters = Gaus_import{ci,rr}.data(:,2)*m + b;
      Adjusted_Gaus_import{ci,rr}.data(:,2) = adjustedCenters;
    end
  end
  
  % ii) Chromatograms: shift data points
  x = -4:fraction_number+5;
  adjusted_raw_data = cell(Nchannels,1);
  for ii = 1:Nchannels
    adjusted_raw_data{ii} = nan(size(cleandata{ii},1),size(cleandata{ii},2)-1);
    for ri=1:size(num_val{1})
      y = cleandata{ii}(ri,2:end);
      y(isnan(y)) = 0;
      rr = replicates(ri);
      b = pfit(ci,rr,1); % intercept
      m = pfit(ci,rr,2); % slope
      x2 = x*m + b;
      y2 = interp1(x,y,x2);
      y2(y2==0) = nan;
      adjusted_raw_data{ii}(ri,:) = y2;
    end
  end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 5. Make some summary statistics
  % Delta_Center
  % Delta_Height
  % Delta_Width
  % EuDist
  % Divergence
  
  tic
  fprintf('    5. Make some summary statistics')
  
  ci = 1;
  
  Delta_center = cell(Nchannels,Nreplicates,Nreplicates);
  Delta_height = cell(Nchannels,Nreplicates,Nreplicates);
  Delta_width = cell(Nchannels,Nreplicates,Nreplicates);
  EuDis = cell(Nchannels,Nreplicates,Nreplicates);
  
  x = 1:user.Nfraction+10;
  
  for ci = 1:Nchannels,
    for rr1 = 1:Nreplicates
      for rr2 = 1:Nreplicates
        
        % i) Find the overlapping proteins
        overlap = intersect(summerised_names_G1{ci,rr1},summerised_names_G1{ci,rr2});
        overlap([1 2]) = [];
        Ia = find(ismember(Gaus_import{ci,rr1}.textdata(:,1),overlap));
        Ib = find(ismember(Gaus_import{ci,rr2}.textdata(:,1),overlap));
        
        % ii) Calculate Gaussian curves
        G1 = zeros(user.Nfraction+10,length(overlap));
        G2 = zeros(user.Nfraction+10,length(overlap));
        for ri = 1:length(overlap)
          c1 = Gaus_import{ci,rr1}.data(Ia(ri),1:3);
          c2 = Gaus_import{ci,rr2}.data(Ib(ri),1:3);
          G1(:,ri) = c1(1)*exp( -(x-(c1(2)+5)).^2 /c1(3).^2 /2);
          G2(:,ri) = c2(1)*exp( -(x-(c2(2)+5)).^2 /c2(3).^2 /2);
        end
        
        % ii3) Calculate statistics
        Delta_center{ci,rr1,rr2} = abs(Gaus_import{ci,rr1}.data(Ia,2) - Gaus_import{ci,rr2}.data(Ib,2));
        Delta_height{ci,rr1,rr2} = abs(Gaus_import{ci,rr1}.data(Ia,1) - Gaus_import{ci,rr2}.data(Ib,1));
        Delta_width{ci,rr1,rr2} = abs(Gaus_import{ci,rr1}.data(Ia,3) - Gaus_import{ci,rr2}.data(Ib,3));
        EuDis{ci,rr1,rr2} = sqrt(sum((G1 - G2) .^ 2));
        
      end
    end
  end
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 6. Write output
  tic
  fprintf('    6. Write output')
  
  writeOutput_alignment
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
  %% 7. Make figures
  
  tic
  fprintf('    7. Make figures')
  
  makeFigures_alignment
  
  tt = toc;
  fprintf('  ...  %.2f seconds\n',tt)
  
  
  
end

diary('off')

