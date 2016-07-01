close all
zslice = 15;
load('mat\bI.mat');
T = load('mat\hmatoof.mat');
T = T.T;
plot_DTI_slice(T, zslice, I);
load('mat\diffussion.mat');
plot_DTI_slice(T, zslice, I);
T = load('mat\hmatvess.mat');
T = T.T;
plot_DTI_slice(T, zslice, I);
