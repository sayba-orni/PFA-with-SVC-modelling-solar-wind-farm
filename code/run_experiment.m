function run_experiment
clc; clear; close all;

% Paths 
thisFile = mfilename('fullpath');
codeDir  = fileparts(thisFile);
projDir  = fileparts(codeDir);          % project root (parent of code/)
addpath(genpath(codeDir))

% Inputs / config 
rng(42,'twister');
C  = config();                          % assumes C.data.* are 'data/pvdata2.csv' etc (under proj root)
Nmc = 200;

% Scenarios 
S = build_scenarios(C, Nmc, 'p_outage', 0.05, 'timeMode', "random");

% Candidates (4 bus + 7 midline) 
k = 0; candidates = struct('type',{},'k',{},'a',{},'b',{});
for b = C.svc_buses
    k = k + 1; candidates(k) = struct('type',"bus",'k',b,'a',[],'b',[]);
end
for e = 1:size(C.mid_list,1)
    k = k + 1; candidates(k) = struct('type',"mid",'k',[],'a',C.mid_list(e,1),'b',C.mid_list(e,2));
end

%  Rank 
lambda = 100;            % MW per unit voltage penalty
alpha  = 0.90;           % CVaR tail
T = rank_candidates(C, S, candidates, lambda, alpha);

% Output folders (UNDER PROJECT ROOT) 
resultsDir = fullfile(projDir,'results');
figDir     = fullfile(projDir,'figures');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
if ~exist(figDir,'dir'),     mkdir(figDir);     end

%  Save table & show top rows 
outCsv = fullfile(resultsDir,'ranking.csv');
writetable(T, outCsv);
disp(T(1:min(10,height(T)),:))
fprintf('Saved ranking table to: %s\n', outCsv);

% Plot CVaR bar chart 
set(0,'DefaultFigureVisible','on');
f = figure('Name','SVC candidate ranking','NumberTitle','off','Color','w','Visible','on');
bar(T.cvarScore);
ax = gca;
ax.XTick = 1:height(T);
ax.XTickLabel = cellstr(T.name);
ax.XTickLabelRotation = 45;
grid on
ylabel(sprintf('CVaR_{%.2f}(Loss + \\lambda·VoltPen)', alpha));
title('Lower is better');
drawnow;

outPng = fullfile(figDir,'ranking.png');
exportgraphics(f, outPng, 'Resolution', 200);
fprintf('Saved figure to: %s\n', outPng);
end
