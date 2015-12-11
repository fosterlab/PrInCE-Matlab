%%%%%%%%%%%%%%% Instructions for running Gauss_Build.m:
%
%
%%%%%%%%%%%%%%% Logic:
% This script just makes figures (at least for now). It fills in the gap between raw instrument data
% and identified proteins; it looks at "what's going on with the peptides?"
% for each isotope channel
%   for each protein group
%       - find every peptide in that group
%       for each replicate
%           - clean the chromatogram
%           - get a "chromatogram" for each peptide
%           for each peptide
%               make a bunch of summary stuff
%           - create figures
%
%
%%%%%%%%%%%%%%% Custom functions called by Gauss_build.m:
%
%
%%%%%%%%%%%%%%% To do:
%
%
%%%%%%%%%%%%%%% Fix the logic:
%
%
%%%%%%%%%%%%%%% Bugs in my (this) code:
%
%
%%%%%%%%%%%%%%% Questions:
% - Where do the input files come from, e.g. 'UP_human_canonical_isoforms_TrEMBL_Oct2013.fasta' and
%   'modificationSpecificPeptides.txt'?
%   - fasta file comes from uniprot
%   - the txt file might be output from maxquant (?)
% - What is 'proteinGroups.txt'?
%   - MQ output
% - What are "peptide results", e.g. from "read in peptide results"?
%   - a huge MQ output file. This is split up into chunks to make for easier reading.
% - Why do you need to process HvsL and MvsL if you're processing HvsM?
% - Why are most of the figures 'HvsM'?
% - How does the fasta file come in?
%