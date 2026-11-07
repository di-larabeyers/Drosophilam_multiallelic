
clear all
close all
clc

% load F matrixes
F_analysis_results = load("F_analysis_results_all.mat");

data_struct = F_analysis_results.F_analysis_results;

data_table = struct2table(data_struct);


% Convert the 'type_top_genotype' column to categorical
chroms = categorical(data_table.chr);
type_top_genotype = categorical(data_table.type_top_genotype);

% Extract the 'pos' column (numerical data)
positions = data_table.pos;

% exract hav1
hav1_values = data_table.hav1;

% keep only valid values (nan and inf come if you try dividing by zero)
valid_indices = ~isinf(hav1_values);

n_wins = length(F_analysis_results.F_analysis_results);

inv_perc = (n_wins - sum(valid_indices))/n_wins * 100;

fprintf("In %g percent of original windows " + ...
    "Inf values occured for hav1.\n", inv_perc);

% filter positions, chromosomes and hav1
positions_cleaned = positions(valid_indices);
chroms_cleaned = chroms(valid_indices);
hav1_cleaned = hav1_values(valid_indices);


%%%%%%%%%%%%%%%%%%%%
%quantiles
% Calculate quantiles for hav1_log
quantiles_hav1 = quantile(hav1_cleaned, [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]);

disp(quantiles_hav1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save csv-files for subgroups

%%%
% type_top_genotype == "hom"

homtg_chr = [];
homtg_pos = [];


% loop over all windows
for i = 1:height(data_table)
    % check if type_top_genotype == "hom"
    if type_top_genotype(i) == 'hom'
        % add chr and pos if true
        homtg_chr = [homtg_chr; data_table.chr(i)];
        homtg_pos = [homtg_pos; data_table.pos(i)];
    end
end

% make new table with filtered values
homtg_table = table(homtg_chr, homtg_pos, ...
    'VariableNames', {'chr', 'pos'});

% save table as csv-file
writetable(homtg_table, 'homtg_hf.csv');

%%%%

% homozygotes are better than heterozygotes on average

indices_homadv = hav1_cleaned < 0;

% extract 'chr' and 'pos' based on indices_homadv
chr_homadv = data_table.chr(indices_homadv);
pos_homadv = data_table.pos(indices_homadv);

% make new table with filtered values
homadv_table = table(chr_homadv, pos_homadv, ...
    'VariableNames', {'chr', 'pos'});

% save table as csv-file
writetable(homadv_table, 'homadv_hf.csv');

n_homadv = height(homadv_table);
disp(['Anzahl homadv: ', num2str(n_homadv)]);

n_homtg = height(homtg_table);
disp(['Anzahl homtg: ', num2str(n_homtg)]);


% overlap
[~, idx] = ismember(homtg_table, homadv_table, 'rows');

% amount of overlap
amount_overlap = sum(idx > 0);

disp(['Amount of overlap: ', num2str(amount_overlap)]);

%%%%
% loci with average heterozygous advantage

indices_mid = ((quantile(hav1_cleaned, 0.499)) <= hav1_cleaned) & (hav1_cleaned <= (quantile(hav1_cleaned, 0.501)));

% extract 'chr' and 'pos' based on indices_mid
chr_mid = data_table.chr(indices_mid);
pos_mid = data_table.pos(indices_mid);

% make new table with filtered values
mid_table = table(chr_mid, pos_mid, ...
    'VariableNames', {'chr', 'pos'});

% save table as csv-file
writetable(mid_table, 'mid_hf.csv');

%%%
% outliers

indices_high_hav = hav1_cleaned >= 13;

% extract 'chr' and 'pos' based on indices_high_hav
chr_high_hav = data_table.chr(indices_high_hav);
pos_high_hav = data_table.pos(indices_high_hav);

% make new table with filtered values
high_hav_table = table(chr_high_hav, pos_high_hav, ...
    'VariableNames', {'chr', 'pos'});

% save table as csv-file
writetable(high_hav_table, 'high_hav_hf.csv');


%%%%
% loci with low heterozygous advantage

indices_five = ((quantile(hav1_cleaned, 0.01)) <= hav1_cleaned) & (hav1_cleaned <= (quantile(hav1_cleaned, 0.05)));

% extract 'chr' and 'pos' based on indices_mid
chr_five = data_table.chr(indices_five);
pos_five = data_table.pos(indices_five);

% make new table with filtered values

five_table = table(chr_five, pos_five, ...
    'VariableNames', {'chr', 'pos'});

% save table as csv-file
writetable(five_table, 'five_hf.csv');



%%%%
% loci with strong heterozygous advantage

indices_ninefive = ((quantile(hav1_cleaned, 0.95)) <= hav1_cleaned) & (hav1_cleaned <= (quantile(hav1_cleaned, 0.99)));

% extract 'chr' and 'pos' based on indices_mid
chr_ninefive = data_table.chr(indices_ninefive);
pos_ninefive = data_table.pos(indices_ninefive);

% make new table with filtered values
ninefive_table = table(chr_ninefive, pos_ninefive, ...
    'VariableNames', {'chr', 'pos'});

% save table as csv-file
writetable(ninefive_table, 'ninefive_hf.csv');


%%%
% Lewontin multiallelic condition

% Extract indices based on 'cond' being true
indices_cond_true = data_table.cond;

% Extract 'chr' and 'pos' based on indices_cond_true
chr_cond_true = data_table.chr(indices_cond_true);
pos_cond_true = data_table.pos(indices_cond_true);

% Create a new table with filtered values
cond_true_table = table(chr_cond_true, pos_cond_true, ...
    'VariableNames', {'chr', 'pos'});

% Save the table as a CSV file
writetable(cond_true_table, 'cond_true_hf.csv');