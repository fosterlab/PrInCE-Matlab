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


%%%%%%%%%%%%%%% Small changes:
% 1. Make all file references and directories absolute. Remove 'cd'. Make this directory-independent!
% 2. Replace all '\' with '/'
% 3. Replace all 'if test==1' with 'if test'
% 4. 'for replicate_to_align_against= y_maximum' does not need a for loop.


%%%%%%%%%%%%%%% Big changes:
% 1. Impute multiple missing values (more than just one, as is currently being done)
% 2. Don't write csv files (these seem useful for debugging, but otherwise aren't needed).
% 3. Include an option not to plot (separate plotting script?).
% 4. Why does Gauss_Build clean the chromatograms more thoroughly? Here we just impute single
%    values. Should we do more?


%%%%%%%%%%%%%%% Questions:
% 1. Why do the chromatograms have 5 extra points on either side?
%
%


%% 0. Initialization

% User defined variables
User_alignment_window1=20; %for initial alignment with proteins fitted to a single Gaussian
User_alignment_window2=8; %for Second round of alignment
Number_of_experimental_channels=2;  % Defines the number of experiments to be compared
Alignment_to_user= 'MvsL'; %Define the experimental channel to use for alignment
User_defined_zero_value = 0.2; %lowest value to be shown in Adjusted Chromatograms

% Define folders, i.e. define where everything lives.
maindir = '/Users/Mercy/Academics/Foster/NickCodeData/GregPCP-SILAC/'; % where everything lives
codedir = [maindir 'Code/']; % where this script lives
funcdir = [maindir 'Code/Functions/']; % where small pieces of code live
datadir = [maindir 'DataFiles/']; % where data files live
figdir = [maindir 'Figures/']; % where figures live
% Make folders if necessary
if ~exist(datadir1, 'dir'); mkdir(datadir1); end
if ~exist(datadir2, 'dir'); mkdir(datadir2); end
if ~exist(datadir3, 'dir'); mkdir(datadir3); end
if ~exist(datadir4, 'dir'); mkdir(datadir4); end


% List all input files. These contain data that will be read by this script.
InputFile{1} = [datadir 'Input/Combined_replicates_2014_04_22_contaminates_removed_for_MvsL_scripts.xlsx'];
InputFile{2} = [datadir 'Input/Combined_replicates_2014_04_22_contaminates_removed_for_HvsL_scripts.xlsx'];
InputFile{3} = [datadir 'Input/Combined_replicates_2014_04_22_contaminates_removed_for_HvsM_scripts.xlsx'];
InputFile{4} = [datadir 'Input/SEC_alignment.xlsx'];



%% 1. Read in data

% Import data from input files
[num_val_MvsL,txt_MvsL] = xlsread(InputFile{1}); %Import file raw Maxqaunt output
[num_val_HvsL,txt_HvsL] = xlsread(InputFile{2}); %Import file raw Maxqaunt output
[num_val_HvsM,txt_HvsM] = xlsread(InputFile{3}); %Import file raw Maxqaunt output
[SEC_size_alignment]=xlsread(InputFile{4});
SEC_fit=polyfit(SEC_size_alignment(1,:),SEC_size_alignment(2,:),1);

%Number of fractions
fraction_number=size(num_val_MvsL);
fraction_number(2) = fraction_number(2)-1;

% number of fractions,equal to # fractions + 10
Chromatograms_axis = fraction_number(2)+10;

Image_couter=1; %Start Image counter



%% 2. Clean the chromatograms


