%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtert F_analysis_results_hcf_all.mat (output von scan_dros_analysis) nach einem genomischen Bereich
%
% Bereich:
%   3L:9145838..9182149
%
% Alle Felder der passenden F_analysis_results-Einträge bleiben erhalten:
%   chr
%   pos
%   F
%   w
%   hav1
%   type_top_genotype
%   cond
%   comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc



%% F_analysis_results laden


input_file = "F_analysis_results_nice_all.mat";

load(input_file, "F_analysis_results");

fprintf("F_analysis_results geladen.\n");
fprintf("Anzahl ursprünglicher Einträge: %d\n\n", ...
    numel(F_analysis_results));



%% Genomischen Bereich festlegen


chromosome = "3L";

min_pos = 9145838;

max_pos = 9182149;



%% Chromosomen und Positionen auslesen


all_chr = string({F_analysis_results.chr});

all_pos = [F_analysis_results.pos];



%% Passende Einträge bestimmen


keep = ...
    all_chr == chromosome & ...
    all_pos >= min_pos & ...
    all_pos <= max_pos;



%% F_analysis_results filtern


F_analysis_results_subset = F_analysis_results(keep);



%% Prüfen, ob Daten gefunden wurden


if isempty(F_analysis_results_subset)

    error(...
        "Keine Einträge für %s:%d..%d gefunden.",...
        chromosome,...
        min_pos,...
        max_pos);

end



%% Anzahl gefundener Einträge ausgeben


fprintf(...
    "Gesuchter Bereich: %s:%d..%d\n",...
    chromosome,...
    min_pos,...
    max_pos);

fprintf(...
    "Gefundene F_analysis_results-Einträge: %d\n",...
    numel(F_analysis_results_subset));



%% Unterschiedliche Positionen bestimmen


unique_positions = unique([F_analysis_results_subset.pos]);

n_positions = numel(unique_positions);

fprintf(...
    "Unterschiedliche Positionen/Sites: %d\n\n",...
    n_positions);



%% Positionen anzeigen


disp("Positionen im gefilterten Datensatz:");

disp(unique_positions');



%% Gefilterte Ergebnisse wieder als F_analysis_results speichern


F_analysis_results = F_analysis_results_subset;



%% Neue MAT-Datei speichern


output_file = sprintf(...
    "F_analysis_results_nice_hf_%s_%d_%d_new.mat",...
    chromosome,...
    min_pos,...
    max_pos);

save(...
    output_file,...
    "F_analysis_results",...
    "-v7.3");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Zusammenfassung


disp(" ");
disp("------------------------------------------------------");
disp("Filter erfolgreich abgeschlossen.");
disp(" ");

fprintf(...
    "Bereich: %s:%d..%d\n",...
    chromosome,...
    min_pos,...
    max_pos);

fprintf(...
    "Gespeicherte F_analysis_results-Einträge: %d\n",...
    numel(F_analysis_results));

fprintf(...
    "Unterschiedliche Positionen: %d\n",...
    n_positions);

disp(" ");

disp("Neue Datei:");
disp(output_file);

disp("------------------------------------------------------");