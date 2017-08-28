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

% Display feedback figure
figure(101)
axis([0 1 0 1])
text(.25,.6,'Starting PrInCE...', 'fontsize',20)
text(.15,.5,'Check logfile.txt for progress.','fontsize',18)
set(gca,'xtick',[],'ytick',[])
set(gcf,'units','normalized','position',[.3 .3 .3 .3])

user.maindir = [pwd '/'];
% Add Code and Functions to path
if ~isdeployed
  f1 = [user.maindir 'Code'];
  f2 = [user.maindir 'Code/Functions'];
  addpath(f1,f2)
end
user = readExperimentFile(user);

% Confirm 'user' is properly formatted
user = cleanuser(user);

% Confirm input files are properly formatted
standardinput(user);

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
