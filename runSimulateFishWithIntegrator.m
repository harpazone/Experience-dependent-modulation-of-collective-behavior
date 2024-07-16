%% run simulaiton of a fish with integrator responding to a replayed trajectory of neighbors


% This demo explains how  to run 'SimulateLarvaFishTrajReplayForGitHub.m' code
% %used in Harpaz et al, 'Experience-dependent modulation of collective behavior'. 
% 
%  To run the code, the user needs supply the following input variables:
%
% folderpath – exact path of where the .mat files containing the response functions 
%              and trajectory replay stimuli can be found. [files can be downloaded from repository]
% 
% replay_file_load - name of the trajectory replay file to load. The file
%                   defines the structure of the experiment. For example 5_20_5.mat 
%                   will present the agent with simulates trajectories of 4, 19 and then 4
%                   neighbors. Trajectory replay experiment should match the
%                   frame rate of the simulations
% 
% Optional inputs [default values are in parenthesis match values used in the paper]:
%
% Fs - samples per sec (50)
% bout_rate - in Hz sec (1.65Hz) 
% Speed - BL/s (1.5)
% arena_diam - in cm (6.5)
% T - total simulation time [frames]
% BL - fish body length [in cm] (0.4)
% One_side_blindang (15) degrees s.t. 2*ang is the full blind angle behind the fish
% Wall_th (2) - distance for wall avoidance reaction
% num_rep (24) - number of simulate repetitions
% USEINTEGRATOR (1); 1 simulate with /0 without integrator

% r_m (13); % input scaling factor
% tau_m = (360000*3); leak time constant in seconds
% S0 = (0.6); integrator starting point
% Srest = (0); integrator resting value
% DeltaS = (0.19); Scaling factor that relates integrator values S(t) to slope 
%                  of occupancy responses of fish p1' such that p1' =
%                  p1(1-[(S(t)-S0)/DeltaS])


% PLOT - 1-yes(default), 0-no
% PLOTfs (50) % frame rate of plotted data (compare to Fs of sim)
% SAVE = 1; %save data yes/no
% MOVIE = (false); % either false or name of the file to save


% Saved Output: x,y - position [NxT matrices]
%         Vx,Vy - velocitiy [NxT matrices]
%         speed - norm of velocity [NxT]
%         angle - heading angle [in degrees]
%         State - indicates which behavior was executed ....
%         WallD - distance to closes wall [in cm];
%         relAng, relOri, relDist - NxTxNneighbors of relative angle,
%         relative orientaiton(heading) and relative distance of neighbors
%         from the focal fish
%         visualAngle -NxTxNneighbors - visual angle (in deg) of neighbors
%         internal_state_vlaues - NxT arbitrary units
%         Bouts,Bout_end - NxT mid and end times of bouts
%         all_p - NxT calculated probability to turn right
%         all_ret_ang - NxTx2xNneighbors - the visual angle of each
%                       neighboring fish. Dimension 3 determines if neighbor
%                       was on the left or right of the focal fish.


close all;clear all; 
folderpath = 'put your folder here';
addpath(folderpath);


SAVE = 1; 
PLOT = 1;
MOVIE = 0;
% MOVIE = 0;
%.................Main Function.................



% name of file to load (must be in the path or folderpath);
replay_file_to_load = '5_20_5for_sim';

SimulateLarvaFishTrajReplayForGitHub(folderpath,replay_file_to_load,'PLOT',...
    PLOT,'PLOTfs',50,'SAVE',SAVE,'MOVIE',MOVIE)

%% run simulation of a group of interacting agents with integrator

% This is a demo explains how  to run the 'SimulateLarvaFishWithDensitySwitchesGitHub.m' code
% %used in Harpaz et al, 'Experience-dependent modulation of collective behavior'. 
% 
%  To run the code, the user needs to set the following parameters:
%
% folderpath – exact path where the .mat files containing the response functions 
%              [files can be downloaded from repository]
% 
% Norder - groups size by phase of experiment simulation. vector of 1xd. example [5,20,5] will
%          switch between group of 5 to 20 and back to 5
% 
% switch_times - number of time points for each phase of the simulation.
%                vector of 1xd+1. example [0,2000,4000,6000]
%               will run phase 1 for the first 2000 times points, then
%               switch to the second phase for another 2000 time points
%               etc'

% Optional inputs [default in parenthesis match values used in the paper]:
%
% Fs - samples per sec (50)
% bout_rate - in Hz sec (1.65Hz) 
% Speed - BL/s (1.5)
% arena_diam - in cm (6.5)
% T - total simulation time [frames]
% BL - fish body length [in cm] (0.4)
% One_side_blindang (15) degrees s.t. 2*ang is the full blind angle behind the fish
% Wall_th (2) - distance for wall avoidance reaction
% num_rep (24) - number of simulate repetitions
% USEINTEGRATOR (1); 1 simulate with /0 without integrator

% r_m (13); % input scaling factor
% tau_m = (360000*3); leak time constant in seconds
% S0 = (0.6); integrator starting point
% Srest = (0); integrator resting value
% DeltaS = (0.19); Scaling factor of that relates integrator values S(t) to slope 
%                  of occupency responses of fish p1' such that S(t) = S0+DeltaS 
%                  cause the response slope p1' = 0 (no response to neighbor)

% PLOT - 1-yes(default), 0-no
% PLOTfs (50) % frame rate of ploted data (compare to Fs of sim)
% SAVE = 1 (default); %save data yes/no
% MOVIE = false; %make a movie - false or text of movie name


% Saved Output: x,y - position [NxT matrices]
%         Vx,Vy - velocity [NxT matrices]
%         speed - norm of velocity [NxT]
%         angle - heading angle [in degrees]
%         State - indicates which behavior was axecuted ....
%         WallD - distance to closes wall [in cm];
%         relAng, relOri, relDist - NxTxNneighbors of relative angle,
%         relative orientaiton(heading) and relative distance of neighobrs
%         from the focal fish
%         visualAngle -NxTxNneighbors - visual angle (in deg) of neighbors
%         internal_state_vlaues - NxT arbitrary units
%         Bouts,Bout_end - NxT mid and end times of bouts
%         all_p - NxT calculated probability to turn right
%         all_ret_ang - NxTx2xNneighbors - the visual angle of each
%                       neighboring fish. Dimension 3 determines if neighbor
%                       was on the left or right of the focal fish.


close all; clear all; 


folderpath = 'add your folder here';
addpath(folderpath);

sampling_rate = 50; % should match simulation code (default is 50Hz)
Nmin = 5; % low density
Nmax = 20; % high density

% set experiment switch times
switch_times = [0,20*60*sampling_rate 40*60*sampling_rate];
            
% set group order
Norder = [Nmin Nmax];

MOVIE = false; 

SimulateLarvaFishWithDensitySwitchesGitHub(folderpath,Norder,switch_times,...
    'PLOT',1,'PLOTfs',50,'MOVIE',MOVIE);





