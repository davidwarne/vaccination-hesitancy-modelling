function S = smry(Data)
%Creates summary statistics from provided data. 
%This version simply returns all the data.
%
% Parameters:
%    Data - data structure as produced by simulation function
% Returns:
%    S - a vector of summary statistics to be used in discrepancy metric
%
% Author:
%   David J. Warne (david.warne@qut.edu.au)
%         School of Mathematical Sciences
%         Science and Engineering Faculty
%         Queensland University of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = [Data.C(:);Data.D(:);Data.V1(:);Data.V2(:)]';
