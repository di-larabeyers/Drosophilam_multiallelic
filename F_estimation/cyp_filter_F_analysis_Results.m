%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtert F_analysis_results_hcf_all.mat nach mehreren genomischen Bereichen
%
% Alle F_analysis_results-Einträge, die in mindestens einem der
% angegebenen Bereiche liegen, werden beibehalten.
%
% Beibehaltene Felder:
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



%% Genomische Bereiche festlegen


% Spalten:
%   Chromosom | Startposition | Endposition

regions = {
    "2R", 11117225, 11121884;
    "2R",  6779749,  6781491;
    "2R", 14887171, 14889010;
    "2R", 14879356, 14881372;
    "2R", 14875829, 14877654;
    "3R", 19135158, 19137136;
    "2R", 12566185, 12569061;
    "2R", 19164143, 19178223;
    "2R",  6173684,  6176195
};



%% Chromosomen und Positionen auslesen


all_chr = string({F_analysis_results.chr});
all_pos = [F_analysis_results.pos];



%% Passende Einträge bestimmen


% Zunächst wird kein Eintrag ausgewählt

keep = false(size(all_pos));


% Alle Bereiche einzeln durchsuchen

for i = 1:size(regions,1)

    chromosome = string(regions{i,1});
    min_pos    = regions{i,2};
    max_pos    = regions{i,3};


    % Einträge im aktuellen Bereich suchen

    current_keep = ...
        all_chr == chromosome & ...
        all_pos >= min_pos & ...
        all_pos <= max_pos;


    % Zum Gesamtfilter hinzufügen

    keep = keep | current_keep;


    % Anzahl gefundener Einträge ausgeben

    n_found = sum(current_keep);

    fprintf(...
        "%s:%d..%d -> %d Einträge gefunden\n",...
        chromosome,...
        min_pos,...
        max_pos,...
        n_found);

end



%% F_analysis_results filtern


F_analysis_results_subset = F_analysis_results(keep);



%% Prüfen, ob überhaupt Daten gefunden wurden


if isempty(F_analysis_results_subset)

    error("In keinem der angegebenen Bereiche wurden Einträge gefunden.");

end



%% Unterschiedliche Sites bestimmen


% Da Positionen auf verschiedenen Chromosomen identisch sein könnten,
% werden chr + pos gemeinsam als Site betrachtet.

subset_chr = string({F_analysis_results_subset.chr});
subset_pos = [F_analysis_results_subset.pos]';

site_ID = subset_chr' + ":" + string(subset_pos);

unique_sites = unique(site_ID);

n_sites = numel(unique_sites);



%% Ergebnisse ausgeben


fprintf("\n--------------------------------------------\n");

fprintf(...
    "Insgesamt gefundene F_analysis_results-Einträge: %d\n",...
    numel(F_analysis_results_subset));

fprintf(...
    "Anzahl unterschiedlicher Sites (chr + pos): %d\n",...
    n_sites);



%% Gefilterte Ergebnisse wieder als F_analysis_results speichern


F_analysis_results = F_analysis_results_subset;



%% Neue MAT-Datei speichern


output_file = "F_analysis_results_nice_hf_cyp_new.mat";

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
    "Anzahl untersuchter Bereiche: %d\n",...
    size(regions,1));

fprintf(...
    "Gespeicherte F_analysis_results-Einträge: %d\n",...
    numel(F_analysis_results));

fprintf(...
    "Unterschiedliche Sites (chr + pos): %d\n",...
    n_sites);

disp(" ");

disp("Neue Datei:");
disp(output_file);

disp("------------------------------------------------------");