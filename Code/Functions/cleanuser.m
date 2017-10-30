function user_new = cleanuser(user)
% CLEANUSER Formats 'user' structure for pcpsilac.
%
%   CLEANUSER(user) ensures that fields within the user structure are 
%   formatted correctly. Catches and corrects errors in file names and 
%   experiment parameters such as wrong type, length, or value. user is a 
%   structure made by master script PRINCE.
%
%   See also PRINCE


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
if isempty(user.majorproteingroupsfile)
  disp('')
  disp(' Continuing without a major protein groups file...')
else
  if ~exist([user.maindir '/Input/' user.majorproteingroupsfile],'file')
    warning('\n The following file could not be found: \n %s', user.majorproteingroupsfile)
    disp('\n Continuing without a major protein groups file...')
    user_new.majorproteingroupsfile = '';
  else
    user_new.majorproteingroupsfile = [user.maindir '/Input/' user.majorproteingroupsfile];
  end
end
% user.corumfile
if exist([user.maindir '/Input/' user.corumfile],'file')~=2
  if isempty(user.corumfile)
    
  else
    error('\n The following file could not be found: \n %s', user.corumfile)
  end
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

% ensure that desiredPrecision is a scalar
if ~isnumeric(user.desiredPrecision) || length(user.desiredPrecision)~=1
  error('\n user.desiredPrecision must be a scalar, and it is not.')
end

% ensure that fractions are between 0 and 1, not 0% and 100%

% ensure that user.Dilution_factor exists

% ensure that the silac ratios in user.MQfiles and user.silacratios are the same order

% ensure that treatmentcondition is part of every comparisonpairs

% ensure that user.organism is right

% ensure that user.silacratios corresponds with user.comparisonpairs
if ~isempty(user.comparisonpairs) && isempty(intersect(user.silacratios,user.comparisonpairs))
  warning('user.silacratios does not correspond with user.comparisonpairs.')
  warning('Setting user.comparisonpairs to first two entries in user.silacratios...')
  user.comparisonpairs = [];
  for ii = 1:min(2,length(user.silacratios))
    user.comparisonpairs{ii} = user.silacratios{ii};
  end
end



%% Add user.maindir/Code and user.maindr/Code/Functions to path

if ~isdeployed
  f1 = [user.maindir 'Code'];
  f2 = [user.maindir 'Code/Functions'];
  if ~exist(f1, 'dir'); mkdir(f1); end
  if ~exist(f2, 'dir'); mkdir(f2); end
  addpath(f1,f2)
end

