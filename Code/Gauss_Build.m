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
InputFile{1} = user.MQfiles{1};
InputFile{2} = user.MQfiles{2};


Nchannels = length(experimental_channels);

% Define folders, i.e. define where everything lives.
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir0 = [maindir 'Data/']; % where data files live
datadir = [maindir 'Data/GaussBuild/']; % where data files live
datadir1 = [maindir 'Data/GaussBuild/Output_Chromatograms/'];
datadir2 = [maindir 'Data/GaussBuild/Output_Chromatograms_filtered_out/'];
datadir3 = [maindir 'Data/GaussBuild/OutputGaus/'];
datadir4 = [maindir 'Data/GaussBuild/OutputGaus_filtered_out/'];
figdir = [maindir 'Figures/']; % where figures live
% Make folders if necessary
if ~exist(datadir0, 'dir'); mkdir(datadir0); end
if ~exist(datadir, 'dir'); mkdir(datadir); end
if ~exist(datadir1, 'dir'); mkdir(datadir1); end
if ~exist(datadir2, 'dir'); mkdir(datadir2); end
if ~exist(datadir3, 'dir'); mkdir(datadir3); end
if ~exist(datadir4, 'dir'); mkdir(datadir4); end
if ~exist(figdir, 'dir'); mkdir(figdir); end
if ~exist([figdir '/GaussBuild'], 'dir'); mkdir([figdir '/GaussBuild']); end



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
  [rawdata{ii},txt_val{ii}] = xlsread(InputFile{ii});
  
  % Remove first column as the replicate
  replicate = rawdata{ii}(:,1);
  rawdata{ii} = rawdata{ii}(:,2:end);
  
  % turn txt_val into a list of protein names
  %txt_val{ii} = txt_val{ii}(:,1);
end
SEC_size_alignment = xlsread(user.calfile);


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
    [clean_chromatogram, x1, x2] = cleanChromatogram(raw_chromatogram);
    
    % store clean_chromatogram in cleandata
    cleandata{ci}(ri,:) = clean_chromatogram;
    
    % store housekeeping variables
    tmp1{ci}(ri,:) = x1;
    tmp2{ci}(ri,:) = x2;
  end
end

tt = toc;
fprintf('  ...  %.2f seconds\n',tt)



%% 3. Fit 1-5 Gaussians on each cleaned chromatogram
tic
fprintf('\n    3. Fit 1-5 Gaussians on each cleaned chromatogram')

% pre-allocate some variables
Coef = cell(Nchannels, Nproteins);
SSE = nan(size(Coef));
adjrsquare = nan(size(Coef));
%fit_flag = nan(size(Coef));
Try_Fit = zeros(size(Coef));

gausscount = 0;
for ci = 1:Nchannels % loop over channels
  parfor ri = 1:Nproteins % loop over proteins
    
    % get a single clean chromatogram
    clean_chromatogram = cleandata{ci}(ri,:);
    
    % don't fit Gaussians if there are less than 5 good data points in the chromatogram
    if sum(clean_chromatogram > 0.05)<5
      continue
    end
    Try_Fit(ci,ri)=1;
    
    % choose 5 models and fit the best one
    %model = choosemodel_holdout(clean_chromatogram,MaxIter);
    model = choosemodel_AIC(clean_chromatogram,Xclean,'AICc');
    disp(['    fit protein number ' num2str(ri) ' (' txt_val{1}{ri+1} ') with ' num2str(model.Ngauss) ' Gaussians...'])
    %[Coef{ci,ri},SSE(ci,ri),adjrsquare(ci,ri),fit_flag(ci,ri)] = fitgaussmodel(clean_chromatogram,model);
    Coef{ci,ri} = model.coeffs;
    SSE(ci,ri) = model.SSE;
    adjrsquare(ci,ri) = model.adjrsquare;
    
%     figure,hold on
%     x = -4:60;
%     plot(x,clean_chromatogram)
%     yhat = feval(model.curveFit,x);
%     plot(x,yhat,'r')
%     model.BIC
%     pause
%     close all
    
  end
end

% Make protgausI{ci}, where each row is for a Gaussian: [protein number, gaussian number, replicate number]
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
  Ngauss(ci) = gausscount;
end

tt = toc;
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
