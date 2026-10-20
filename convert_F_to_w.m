function w = convert_F_to_w(F)
%Converts a matrix of fitness effects F into a matrix of relative fitnesses
% w, where the genotype with the highest fitness has rel fitness of 1.0

max_fitness = max(max(F+1.0));
w = (F+1.0) ./ max_fitness;

end