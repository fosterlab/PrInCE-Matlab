%%%%%%%%%%%%%%% Instructions for running Gauss_Build.m:
% - Rename maindir for your computer. (The command 'pwd' will give you the full path to a folder.)
% - Add maindir/Code and maindir/Code/Functions to Matlab path.
%
%%%%%%%%%%%%%%% Logic:
% 0. Initalize
% 1. Read input data (MaxQuant output)
% 2. Clean the chromatograms
% 3. Fit 1-5 Gaussians on each cleaned chromatogram
% 4. Write output
%
%%%%%%%%%%%%%%% Custom functions called by Gauss_build.m:
% - cleanChromatogram.m
% - choosemodel_holdout.m
% - gaussfitICs.m
% - fitgaussmodel.m
%
%%%%%%%%%%%%%%% To do:
% - Include the parfoor loops
% - Right now the logic is a bit roundabout:
%       First select a model and then fit that model. This is dangerous because the final fit
%       could be bad! Better to keep the parameters of the best fit from model selection.
% - Add 'user defined settings' at the very top of initialize
% - Fix InputFile, MainOutputFile, and DebugOutputFile which are not currently used.

diary([user.maindir 'logfile.txt'])
disp('Gauss_Build.m')



%% 0. Initialize
tic
fprintf('\n    0. Initialize')


% Load user settings
maindir = user.maindir;
experimental_channels = user.silacratios;


Nchannels = length(experimental_channels);

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
%[num_val_MvsL,txt_MvsL] = xlsread(InputFile{1}); %Import file MvsL
%[num_val_HvsL,txt_HvsL] = xlsread(InputFile{2}); %Import file HvsL
rawdata = cell(size(experimental_channels));
txt_val = cell(size(experimental_channels));
for ii = 1:Nchannels
  %[rawdata{ii},txt_val{ii}] = xlsread(user.MQfiles{ii});
  %tmp = importdata(user.MQfiles{ii});
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
  
  % turn txt_val into a list of protein names
  %txt_val{ii} = txt_val{ii}(:,1);
end
%SEC_size_alignment = xlsread(user.calfile);


% How many proteins and fractions are there?
[Nproteins, Nfractions] = size(rawdata{1});


tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 2. Clean the chromatograms
%   1. Impute (fill in) single missing values (nans) by linear interpolation;
%   2. Fill in first and last fractions with a 0 if they're missing;
%   3. Replace values <0.2 with nan;
%   4. If a value is not part of 5 consecutive values, replace it with 0.05;
%   5. Add 5 zeros to either side of the chromatogram;
%   6. Smooth whole chromatogram with a boxcar filter (moving average).
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
    clean_chromatogram = cleanChromatogram2(raw_chromatogram);
    
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
%fit_flag = nan(size(Coef));
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
    if sum(clean_chromatogram > 0.05)<5
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
% It's just an indexing variable, that matches up gaussians to protein number.
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
