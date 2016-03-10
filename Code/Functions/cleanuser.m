function user_new = cleanuser(user)
% This little function does a couple things.
% First, it checks that the user structure in pcpsilac_master.m is formatted correctly.
%   1. Does user.maindir end in a forward slash?
%   2. Do the data files exist? If not, are there similarly-named data files?
% Second, it adds user.maindir/Code and user.maindr/Code/Functions to path

user_new = user;

%% Check that user is formatted correctly

% check that user.maindir ends in a forward slash
if user.maindir(end) ~= '/'
  user_new.maindir = [user.maindir '/'];
end

% do the data files exist?
% user.MQfiles
% user.calfile
% user.majorproteingroupsfile
% user.corumfile
% user.mastergaussian
% user.fastafile
% user.omimfile

% ensure that user.silacratios is a cell, not a string. this is a problem when Nchannels=1.
if ischar(user.silacratios)
  user_new.silacratios = {user.silacratios};
end

% ensure that fractions are between 0 and 1, not 0% and 100%


%% Add user.maindir/Code and user.maindr/Code/Functions to path

f1 = [user.maindir 'Code'];
f2 = [user.maindir 'Code/Functions'];
if ~exist(f1, 'dir'); mkdir(f1); end
if ~exist(f2, 'dir'); mkdir(f2); end

addpath(f1,f2)

