function F_ij = convert_F_to_F_ij(F, i, j)
%Converts matrix of fitness effects F into evolutionarily equivalent
% matrix F_ij, where entry F_ij becomes zero
k = 1/(1+F(i,j));
F_ij = (k-1) + k * F;

end