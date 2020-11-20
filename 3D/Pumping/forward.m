function [OBS,A_hard] = forward(TRUE)
nx = 40;                                                                   % number of blocks in each direction
ny = 120;
nz = 40;                                       % size of blocks in each direction (m)
trueplot = flip(permute(reshape(TRUE,[nx,ny,nz]),[3,2,1]),3);
[P_true] = pressure_calculation(trueplot);

load wellloc wellloc
IND = wellloc;
OBS_P = P_true(IND);

A_hard = zeros(length(IND),length(TRUE(:)));
for i=1:length(IND)
    A_hard(i,IND(i)) = 1;
end
OBS_HD = A_hard*TRUE(:);

OBS = [OBS_P;OBS_HD];
end

