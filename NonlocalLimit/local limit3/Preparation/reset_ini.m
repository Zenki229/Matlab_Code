clc;clear;close all;

% Potential
q = @(x) 2+ sin (5*pi.*x);
save('pot','q');
% Source
f = @(x) 10 + 0.*x;
save('soce','f');
