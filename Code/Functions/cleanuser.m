function user_new = cleanuser(user)
% This little function does a couple things.
% First, it checks that the user structure in pcpsilac_master.m is formatted correctly.
%   1. Does user.maindir end in a forward slash?
%   2. Do the data files exist? If not, are there similarly-named data files?
% Second, it adds user.maindir/Code and user.maindr/Code/Functions to path

user_new = user;


%% Check that files exist

% check that user.maindir ends in a forward slash
if user.maindir(end) ~= '/'
  user_new.maindir = [user.maindir '/'];
end

% check that nickflag exists, if not set it to zero
if ~isfield(user,'nickflag')
  user_new.nickflag = 0;
end

% do the data files exist?
% user.MQfiles
for ii = 1:length(user.MQfiles)
  if ~exist([user.maindir '/Input/' user.MQfiles{ii}],'file')
    error('\n The following file could not be found: \n %s', user.MQfiles{ii})
  else
    user_new.MQfiles{ii} = [user.maindir '/Input/' user.MQfiles{ii}];
  end
end
% user.majorproteingroupsfile
if ~exist([user.maindir '/Input/' user.majorproteingroupsfile],'file')
  error('\n The following file could not be found: \n %s', user.majorproteingroupsfile)
else
  user_new.majorproteingroupsfile = [user.maindir '/Input/' user.majorproteingroupsfile];
end
if ~exist([user.maindir '/Input/' user.corumfile],'file')
  error('\n The following file could not be found: \n %s', user.corumfile)
else
  user_new.corumfile = [user.maindir '/Input/' user.corumfile];
end

% Make user.corumpairwisefile
fn_corumpair = [user.maindir '/Output/tmp/Corum_pairwise.csv'];
if ~exist(fn_corumpair,'file')
  disp('cleanuser: Making Corum_pairwise.csv.')
  corum2pairwise(user_new)
end
user_new.corumpairwisefile = fn_corumpair;

% Make user.corumcomplexfile
fn_corumcomplex = [user.maindir '/Output/tmp/Corum_complex.csv'];
try
  corumextractcomplex(user_new)
end
user_new.corumcomplexfile = fn_corumcomplex;



%% Check that user is formatted correctly

% ensure that user.silacratios is a cell, not a string. this is a problem when Nchannels=1.
if ischar(user.silacratios)
  user_new.silacratios = {user.silacratios};
end

% ensure that fractions are between 0 and 1, not 0% and 100%

% ensure that user.Dilution_factor exists

% ensure that the silac ratios in user.MQfiles and user.silacratios are the same order

% ensure that treatmentcondition is part of every comparisonpairs

% ensure that 



%% Add user.maindir/Code and user.maindr/Code/Functions to path

f1 = [user.maindir 'Code'];
f2 = [user.maindir 'Code/Functions'];
if ~exist(f1, 'dir'); mkdir(f1); end
if ~exist(f2, 'dir'); mkdir(f2); end

addpath(f1,f2)

