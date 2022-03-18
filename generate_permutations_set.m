% USE PALM WITH FIXELCFESTATS (FIXEL-BASED ANALYSIS)
% refer to: https://mrtrix.readthedocs.io/en/latest/reference/commands/fixelcfestats.html
% prereq:
% (1) download PALM and add to Matlab path; save new dir in local user Documents dir. 
% (2) exchangeability block file exchangeability_blocks.csv table that has the structure:
%       4 columns (x1, x2, x3, x4) by N rows (where N=number of subjects)
%           x1: all negative one (-1) so that blocks of data are exchangeable within a block only
%           x2: number of siblings the subject has in the study including that subject (i.e. minimum = 1)
%           x3: index number for each family with a specific number of siblings
%           x4: index number for each sibling in each family
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make m x n matrix in matlab/octave %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refer to: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/ExchangeabilityBlocks%Generating_the_set_of_permutations_only
% Do not need to run PALM on the diffusion data--rather, use PALM to Generate the set of permutations only (i.e. text file defining the siblings)
% load exchangeability block file
EB=csvread('P:/3022028.01/fixel_based_analysis/fba_long_w1_w2/fixelcfestats/palm/exchangeability_blocks.csv');
% generate the set of permutation only 
[Pset,VG]=palm_quickperms([],EB,5000);
% save permutations output as a text file that can be imported elsewhere
dlmwrite('P:/3022028.01/fixel_based_analysis/fba_long_w1_w2/fixelcfestats/palm/permutations_set_5000.csv',Pset);