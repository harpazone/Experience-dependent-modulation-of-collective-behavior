function [output] = collectTrajReplaySimulationData(varargin)
%description: collects all simulated trajectory replay data and saves it in
%format matching the real experiments so we can use the exact same codes
%for analysis

%Input: 



%Output: 

default_Folder  = 'G:\simulations\density_habituation/ab_traj_replay_7_neigh_wall1/';

default_Tframes  = 100*60;


% parse
vars = inputParser;
addParameter(vars,'Folder',default_Folder);
addParameter(vars,'Tframes',default_Tframes);

parse(vars,varargin{:})

Folder = vars.Results.Folder;
Tframes = vars.Results.Tframes;

%% ...........Local Variable definitions..........
% close all; clear all;

% Path = 'G:\simulations\density_habituation/';

% Folder = 'ab_traj_replay_7_neigh_wall1/';

arena_diam = 6.5;
simFR=50;
cd([Folder]);


file_list = dir('sim_larva*');
file_list = {file_list.name};

Nfish = length(file_list);
Tframes = Tframes*simFR;

% make variables
x = zeros(Tframes,Nfish);
y = zeros(Tframes,Nfish);
time = zeros(Tframes,Nfish);
ori_smoothed = zeros(Tframes,Nfish);
bout_start = zeros(Tframes,Nfish);
bout_end = zeros(Tframes,Nfish);
bout_start_smoothed_ori = zeros(Tframes,Nfish);
bout_end_smoothed_ori = zeros(Tframes,Nfish);
bout_start_path = zeros(Tframes,Nfish);
bout_end_path = zeros(Tframes,Nfish);
movie_fr = [];
camera_fr = [];
fish_path = zeros(Tframes,Nfish);
ori = zeros(Tframes,Nfish);
fish_integrator_potential = zeros(Tframes,Nfish);

for i = 1: length(file_list)
    cd([Folder,file_list{i}]);
    load('coor','xfish','yfish');
    to_lose = xfish==1;
    
    xfish(to_lose) = nan;
    yfish(to_lose) = nan;
    
    load('Header','Fs','arena_diam')
    load('group_prop');
    
    angle(to_lose) = nan;
    Speed(to_lose) = nan;
    
    T = Tframes/Fs;
    x(:,i) = xfish/(arena_diam/2);
    y(:,i) = yfish/(arena_diam/2);
    tt = 1/Fs:1/Fs:T;
    time(:,i) = tt;
    
    ori_smoothed(:,i) = angle;
    ori(:,i) = mod(angle,360);
    
    Speed = Speed./(arena_diam/2);
    curr_fish_path = cumsum(Speed);
    fish_path(:,i) = curr_fish_path;
    if exist('integrator_potential','var')
        fish_integrator_potential(:,i) = integrator_potential;
    end
    
    % bout_vec = Speed>0;
    % [st,nd] = findNonZeroSeq(bout_vec);
    % st = st-1;
    % nd = nd+1;

    st = find(Bouts)-1;
    nd = find(Bouts_end)+1;

    
    st(st<1) =1;
    nd(nd>length(Bouts)) = length(Bouts);
    
    bout_start(1:length(st),i) = st/Fs;
    bout_end(1:length(nd),i) = nd/Fs;
    
    bout_start_smoothed_ori(1:length(st),i) = angle(st);
    bout_end_smoothed_ori(1:length(nd),i) = angle(nd);
    
    bout_start_path(1:length(st),i) = curr_fish_path(st);
    bout_end_path(1:length(nd),i) = curr_fish_path(nd);
    
    

    i
end


cd([Folder]);


save('data','x','y','time','ori_smoothed','ori','bout_start','bout_end',...
    'bout_start_smoothed_ori','bout_end_smoothed_ori','bout_start_path',...
    'bout_end_path','movie_fr','camera_fr','fish_path','fish_integrator_potential');

disp(['simulation saved from folder: ',Folder])
% cd(Path)
% close all; clear all;
%.................Main Function.................



%............Call for local functions...........




