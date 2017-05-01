%GAUSSBUILD Find peaks in co-fractionation profiles.
%   GAUSSBUILD identifies peaks in co-fractionation profiles by
%   deconvolving profiles into Gaussian mixture models. Models are mixtures 
%   of between 1 and 5 Gaussians. Profiles are pre-processed (single data
%   point imputation, smoothing) before model fitting. Model fitting and
%   model selection are performed by CHOOSEMODEL_AIC, which uses the 'fit'
%   command.
% 
%   Master script PRINCE must be run before GAUSSBUILD.
%
%   GAUSSBUILD produces two output folders: Output/Data/GaussBuild, which
%   contains all output csv tables, and Output/Figures/GausBuild, which
%   contains all figures.
%
%   Workflow:
%   1. Read co-fractionation profiles
%   2. Pre-process profiles
%   3. For all profiles, fit all five models and choose the best model
%   4. Write output tables
%   5. Make figures
%
%   See also PRINCE, CLEANPROFILE, CHOOSEMODEL_AIC.

diary([user.maindir 'logfile.txt'])
disp('Gauss_Build.m')


%% 0. Initialize
tic
fprintf('\n    0. Initialize')

% Load user settings
maindir = user.maindir;
experimental_channels = user.silacratios;

Nchannels = length(experimental_channels);

minGoodValues = 5;

% Define folders, i.e. define where everything lives.
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir = [maindir 'Output/Data/GaussBuild/']; % where data files live
figdir = [maindir 'Output/Figures/GaussBuild/']; % where figures live
% Make folders if necessary
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist([maindir '/Output'], 'dir'); mkdir([maindir '/Output']); end
if ~exist([maindir '/Output/Data'], 'dir'); mkdir([maindir '/Output/Data']); end
if ~exist([maindir '/Output/Figures'], 'dir'); mkdir([maindir '/Output/Figures']); end
if ~exist([maindir '/Output/tmp'], 'dir'); mkdir([maindir '/Output/tmp']); end
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist(figdir, 'dir'); mkdir(figdir); end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)


%% 1. Read input data
tic
fprintf('\n    1. Read input data')

% Import data from all input files
rawdata = cell(size(experimental_channels));
txt_val = cell(size(experimental_channels));
for ii = 1:Nchannels
  %[rawdata{ii},txt_val{ii}] = xlsread(user.MQfiles{ii});
  tmp = readchromatogramfile2(user.MQfiles{ii});
  
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
  
  rawdata{ii} = tmp.data;
  txt_val{ii} = tmp.textdata;
  
  % if txt_val includes protein groups, reduce it to the first protein in each group
  for jj = 1:size(txt_val{ii},1)
    tmp = strsplit(txt_val{ii}{jj},';');
    txt_val{ii}{jj} = tmp{1};
  end
  
  % if rawdata & txt_val are the same length, assume they both have headers, remove rawdata header
  if size(rawdata{ii},1)==size(txt_val{ii},1)
    rawdata{ii} = rawdata{ii}(2:end,:);
  end
  
  % Remove first column of rawdata as the replicate
  if sum(mod(rawdata{ii}(:,1),1)~=0)>0
    disp('Warning: Gauss_Build: Replicate column in chromatogram tables is badly formatted.')
    disp('Warning: Gauss_Build: Assuming all chromatograms are from a single replicate...')
    replicate = ones(size(rawdata{ii},1),1);
  else
    replicate = rawdata{ii}(:,1);
    rawdata{ii} = rawdata{ii}(:,2:end);
  end
  
end

% How many proteins and fractions are there?
[Nproteins, Nfractions] = size(rawdata{1});

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 2. Clean the chromatograms
%   1. Impute (fill in) single missing values (nans) by linear interpolation.
%   2. Remove "lone" singletons and doubletons.
%   3. Replace all missing values with near-zero noise.
%   4. Add 5 near-zeros to either side.
tic
fprintf('\n    2. Clean the chromatograms')

% store all clean chromatograms in cleandata
cleandata = cell(size(rawdata));
Xraw = (1:Nfractions)';
Xclean = (1:Nfractions+10)' - 5;

% a couple of housekeeping variables
tmp1 = cell(size(rawdata));
tmp2 = cell(size(rawdata));

for ci = 1:Nchannels % loop over channels
  cleandata{ci} = zeros(size(rawdata{ci},1),size(rawdata{ci},2)+10);
  for ri = 1:Nproteins % loop over proteins
    
    raw_chromatogram = rawdata{ci}(ri,1:Nfractions);
    clean_chromatogram = cleanProfile(raw_chromatogram,1);
    
    % store clean_chromatogram in cleandata
    cleandata{ci}(ri,:) = clean_chromatogram;
    
  end
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 3. Fit 1-5 Gaussians on each cleaned chromatogram
t1=tic;
fprintf('\n    3. Fit 1-5 Gaussians on each cleaned chromatogram')

% pre-allocate some variables
Coef = cell(Nchannels, Nproteins);
SSE = nan(size(Coef));
adjrsquare = nan(size(Coef));
Try_Fit = zeros(size(Coef));

gausscount = 0;
t2=tic;
for ci = 1:Nchannels % loop over channels
  cleandata_nonbc = cleandata{ci};
  rawdata_nonbc = rawdata{ci};
  txt_val_nonbc = txt_val{1};
  parfor ri = 1:Nproteins % loop over proteins
    
    % get a single clean chromatogram
    clean_chromatogram = cleandata_nonbc(ri,:);
    
    % don't fit Gaussians if there are less than 5 good data points in the chromatogram
    if sum(clean_chromatogram > 0.05)<minGoodValues
      continue
    end
    Try_Fit(ci,ri)=1;
    
    % Throw out any models where the number of non-imputed data points, i.e. the real data, is less than
    % the number of parameters.
    Ngaussmax = floor(sum(rawdata_nonbc(ri,1:Nfractions)>0.01)/3);
    Ngaussmax = min([5 Ngaussmax]);
    
    % fit 5 models and chose the best one
    %model = choosemodel_holdout(clean_chromatogram,MaxIter);
    try
      model = choosemodel_AIC(clean_chromatogram,Xclean,'AICc',Ngaussmax);
      Coef{ci,ri} = model.coeffs;
      SSE(ci,ri) = model.SSE;
      adjrsquare(ci,ri) = model.adjrsquare;
      
      % if bad fit, try again
      maxIter = 2;
      iter = 0;
      while adjrsquare(ci,ri)<0.85 && iter<=maxIter
        iter = iter+1;
        fprintf(['\n    re-fitting ' txt_val_nonbc{ri+1} '...'])
        model = choosemodel_AIC(clean_chromatogram,Xclean,'AICc',Ngaussmax);
        Coef{ci,ri} = model.coeffs;
        SSE(ci,ri) = model.SSE;
        adjrsquare(ci,ri) = model.adjrsquare;
      end
      
      fprintf(['\n    fit ' txt_val_nonbc{ri+1} ' with ' num2str(model.Ngauss) ' Gaussians, R^2=' num2str(round(adjrsquare(ci,ri)*100)/100)])
      
    catch
      Try_Fit(ci,ri)=0;
    end
        
  end
end

% Make protgausI{ci}, where each row is a Gaussian: [protein number, gaussian number, replicate number]
% i.e. an indexing variable between gaussians and proteins.
protgausI = cell(Nchannels,1); 
Ngauss = zeros(Nchannels,1);
for ci = 1:Nchannels
  protgausI{ci} = zeros(100,3);
  gausscount = 0;
  for ri = 1:Nproteins % loop over proteins
    % don't fit Gaussians if there are less than 5 good data points in the chromatogram
    if sum(cleandata{ci}(ri,:) > 0.05)>5
      Ng = length(Coef{ci,ri})/3;
      rep = replicate(ri);
      for gi = 1:Ng
        gausscount = gausscount+1;
        protgausI{ci}(gausscount,:) = [ri, gi, rep];
      end
    end
  end
  protgausI{ci} = protgausI{ci}(1:gausscount,:);
  Ngauss(ci) = gausscount;
end

tt = toc(t1);
fprintf('  ...  %.2f seconds\n',tt)



%% 4. Write output
tic
fprintf('\n    4. Write output')

writeOutput_gaussbuild

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 5. Make figures
tic
fprintf('\n    5. Make figures')

makeFigures_gaussbuild

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)

diary('off')
