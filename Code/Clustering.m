%%%%%%%%%%%%%%% Instructions for running:
%
%
%%%%%%%%%%%%%%% Logic:
% Identifies protein complexes, as opposed to just binary interactions. Uses a markov-chain
% clustering algorithm, but Nick says similar results could be obtained from just looking for
% Gaussians with approximately equal Centers.
% for each precision level (70%, 60%, 40%)
%   for each Gaussian
%       if this Gaussian is in an interaction
%           - add it to the interaction matrix
%   - perform MCL on the interaction matrix
%
%
%%%%%%%%%%%%%%% Custom functions called:
%
%
%%%%%%%%%%%%%%% To do:
% - According to the mcl doc, it's rare to find overlapping clusters. But Nick's code finds a lot of
%   overlapping clusters. What's going on?
%
%
%%%%%%%%%%%%%%% Fix the logic:
%
%
%%%%%%%%%%%%%%% Bugs in my (this) code:
%
%
%%%%%%%%%%%%%%% Questions:
% 1. How does the normalization (for the interaction matrix) work?
%