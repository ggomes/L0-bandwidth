function []=add_dependencies()

root = fullfile(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(root,'src')))