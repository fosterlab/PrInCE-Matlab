%PRINCE Master script for the PrInCE bioinformatics pipeline
%   PRINCE.m collates information about the co-fractionation experiment,
%   such as the number of fractions, replicates, and biological conditions,
%   as well the location of input data files and the reference database of
%   gold-standard protein complexes. These experimental parameters are
%   combined in the user structure, which controls the behaviour of the 5
%   PrInCE modules: GAUSSBUILD, ALIGNMENT, FOLDCHANGES, INTERACTIONS, and
%   COMPLEXES.
%
%   Correct formatting of experiment information is confirmed by the
%   function CLEANUSER.
%
%   See https://github.com/fosterlab/PRInCE for use instructions.
%
%   See also GAUSSBUILD, ALIGNMENT, FOLDCHANGES, INTERACTIONS, COMPLEXES, CLEANUSER.


clear user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data files

user.maindir = '/path/to/your/analysis/folder/';
user.MQfiles{1} = 'name_of_cofractionation_file_1.csv';
user.MQfiles{2} = 'name_of_cofractionation_file_2.csv';
user.majorproteingroupsfile = 'name_of_majorproteingroups_file.csv';
user.corumfile = 'allComplexes.csv';


% Experimental design
% Names of SILAC ratios
user.silacratios = {'MvsL' 'HvsL'};
% Name of treated and non-treated ratios. Leave blank if no treatment.
user.treatmentcondition = 'HvsL';
user.notreatment = 'MvsL';
% Number of fractions
user.Nfraction = 55;
% Number of replicates
user.Nreplicate = 3;
% Model organism. This must match corum file.
user.organism = 'human';


% Script parameters

% GaussBuild
user.fastgaussbuild = 0;        % Set this to 1 to skip making plots for individual proteins.
                                % Set this to 0 to make the plots, which takes a lot of time.
% Alignment
user.skipalignment = 1;
user.User_alignment_window1 = 20; % The maximum distance between Gaussian centers, in fractions, that
                                % will be considered for alignment. Larger value means more data 
                                % and a looser alignment

% FoldChanges
user.skipcomparison = 0;
user.userwindow = 2;            % Only compare Gaussians that are within this many fractions of each other.
user.comparisonpairs = {'MvsL' 'HvsL'}; % Which pairs of conditions should be used?
user.fastcomparison = 1;        % Set this to 1 to skip making plots for individual proteins.
                                % Set this to 0 to make the plots, which takes a lot of time.

% Interactions
user.desiredPrecision = 0.5;   % At what precision do you want to evaluate protein interactions? This
                               % should be a number between 0 and 1.

% Complexes
user.separateByReplicate = 0;
user.separateByChannel = 1;
user.fdr = 0.01; % FDR cutoff for enrichment analysis


% Confirm 'user' is properly formatted
if ~exist(user.maindir,'dir')
  user.maindir = pwd;
  user.maindir = [user.maindir '/'];
end
f1 = [user.maindir 'Code'];
f2 = [user.maindir 'Code/Functions'];
addpath(f1,f2)
user = cleanuser(user);

% Confirm input files are properly formatted
standardinput(user);

user.nickflag = 0;

% Save 'user'
save([user.maindir 'userconfig.mat'],'user')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GaussBuild
%Alignment
%FoldChanges
%Interactions
Complexes
