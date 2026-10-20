% Loops over (part of) drosophila genome and analyzes F for each locus

clear all
close all
clc

% load F matrixes
F_results_unfiltered = load("F_results_hcf_all_new.mat").F_results;

% filter:
% Exclude loci where optimization failed %unfiltered only for hcf
keep_indices = ~isnan([F_results_unfiltered.cost]);
F_results = F_results_unfiltered(keep_indices);
opt_failed = isnan([F_results_unfiltered.cost]);
n_opt_failed = length(F_results_unfiltered(opt_failed));
p_opt_failed = 100 * (n_opt_failed / length(F_results_unfiltered));

% Exclude loci where freqs where all 0.25, F only zeros and cost = 0
F_results = F_results([F_results.cost] ~= 0);

% extract and save number of loci
n_wins = length(F_results);
disp(n_wins);

% save analysis results
temp_struct = struct('chr', [], ...
    'pos', [], ...
    'F', [], ...
    'w', [], ...
    'hav1', [], ...
    'type_top_genotype', [], ...
    'cond', false, ...
    'comment', []);
F_analysis_results = repmat(temp_struct, n_wins, 1);

% initialise with 0
n_het_top = 0;
n_condtrue = 0;

%loop through all loci
for i = 1:n_wins
    chr = F_results(i).chr; % extract chromosome name for this locus
    pos = F_results(i).pos; % extract position for this locus

    F = F_results(i).F; % extract F for this locus
    
    w = convert_F_to_w(F); % convert F to w matrix (relative fitness)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate "strength" of heterozygous advantage
    fhet = 0;
    fhom = 0;
    fhet_vektor = [];
 
    [m, n] = size(w); % m = number of lines, n = number of columns; in our case always m = n
    for j = 1:m
        fhom = fhom + w(j, j); % Accumulate the diagonal elements of w
        for k = j+1:m % only consider heterozygote values once, since matrix is symmetrical (k > j)
            fhet = fhet + w(j, k); % accumulate values
            fhet_vektor(end+1) = w(j, k);
        end
    end
    nhet = n*(n - 1)/4;
    dhet = fhet/nhet; % average fitness of heterozygotes
    dhom = fhom/n; % average fitness of homozygotes

    hav1 = log10(dhet/dhom); % heterozygous advantage value 1
    F_analysis_results(i).hav1 = hav1; % Store the heterozygous advantage value 1
    
   
    % condition from Lewontin paper
    s = 1-(fhom/fhet);
    sigma = std(fhet_vektor);
    condv1 = s/sigma;
    condv2 = 2*sqrt(m);

    if condv1 > condv2
        F_analysis_results(i).cond = true;
        n_condtrue = n_condtrue +1;
    else 
        F_analysis_results(i).cond = false;  
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find type of top genotype for every matrix

    [max_w, linear_index] = max(w(:));

    % Convert linear index to row and column indices
    [max_row, max_col] = ind2sub(size(w), linear_index);

    F_analysis_results(i).chr = chr;
    F_analysis_results(i).pos = pos;
    F_analysis_results(i).F = F;
    F_analysis_results(i).w = w;

    if max_row == max_col
        type_top_genotype = "hom";
    else
        type_top_genotype = "het";
        n_het_top = n_het_top + 1;
    end
    F_analysis_results(i).type_top_genotype = type_top_genotype;

    
end

%calculate percentages
percent_het_top = 100*(n_het_top/n_wins);
percent_n_condtrue = 100*(n_condtrue/n_wins);

fprintf("In %g percent of original loci " + ...
    "optimization failed\n", p_opt_failed);

fprintf("In %g percent of loci the most fit genotype is " + ...
    "a heterozygote\n", percent_het_top);

fprintf("In %g percent of loci " + ...
    "the condition is true\n", percent_n_condtrue);


% filename
file_str = "F_analysis_results_hcf_all_new" + ".mat";

% save matrix as .mat file
save(file_str, 'F_analysis_results');





