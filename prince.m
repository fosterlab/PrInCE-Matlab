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

user.maindir = [pwd '/'];
% Add Code and Functions to path
if ~isdeployed
  f1 = [user.maindir 'Code'];
  f2 = [user.maindir 'Code/Functions'];
  addpath(f1,f2)
end

% Get user input
user = princeGUI;
if isempty(user)
  return;
end
user.maindir = [pwd '/'];

% Confirm 'user' is properly formatted
user = cleanuser(user);

% Confirm input files are properly formatted
user = standardinput(user);

% Save 'user'
save([user.maindir 'userconfig.mat'],'user')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GaussBuild
Alignment
FoldChanges
Interactions
Complexes
