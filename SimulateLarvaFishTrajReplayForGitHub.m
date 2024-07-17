function SimulateLarvaFishTrajReplayForGitHub(loadpath,replay_file_to_load,varargin)
% description: 

% set defaults
defaultFs = 50; % sampling rate Hz
defaultN = 1; % % number of fish
defaultbout_rate = 1.65; %
defaultSpeed = 1.5; % BL/s
defaultBL = 0.4; % fish body length cm;
defaultArena_diam = 6.5; % arena diameter in cm
defaultOne_side_blindang = 15;  % in degrees s.t. 2*ang is the full blind angle
defaultWall_th = 2; % distance to repond to wall in BL
default_num_rep = 1; % number of repetitions
defaultPLOT = 1; % plot data yes/no
defaultPLOTfs = 50; % frame rate of ploted data (compare to Fs of sim)
defualtSAVE = 1; %save data yes/no
defualtUSEINTEGRATOR = 1; %simulate with/without integrator
defualtMOVIE = false; % either false or name of the file to save



defualt_r_m = 13; % input scaling factor
defualt_tau_m = 360000*3; % leak time constant in seconds
defualt_S0 = 0.6; % starting point
defualt_Srest = 0; %resting value
defualt_DeltaS = 0.19; % S(t) = S0+DeltaS causes p1' = 0; (no repsone to neighbors)

% parse
vars = inputParser;
addParameter(vars,'Fs',defaultFs);
addParameter(vars,'N',defaultN);
addParameter(vars,'PLOT',defaultPLOT);
addParameter(vars,'PLOTfs',defaultPLOTfs);

addParameter(vars,'SAVE',defualtSAVE);
addParameter(vars,'all_bout_rate',defaultbout_rate);
addParameter(vars,'all_avg_speed',defaultSpeed);
addParameter(vars,'all_BL',defaultBL);
addParameter(vars,'all_arean_rad',defaultArena_diam);
addParameter(vars,'blindang',defaultOne_side_blindang);
addParameter(vars,'wall_th',defaultWall_th);
addParameter(vars,'USEINTEGRATOR',defualtUSEINTEGRATOR);
addParameter(vars,'num_rep',default_num_rep);
addParameter(vars,'MOVIE',defualtMOVIE);


addParameter(vars,'r_m',defualt_r_m);
addParameter(vars,'tau_m',defualt_tau_m);
addParameter(vars,'S0',defualt_S0);
addParameter(vars,'Srest',defualt_Srest);
addParameter(vars,'DeltaS',defualt_DeltaS);

parse(vars,varargin{:})

% create variables
all_bout_rates = vars.Results.all_bout_rate;
avg_speed = vars.Results.all_avg_speed;
BL = vars.Results.all_BL;
arena_diam = vars.Results.all_arean_rad;
blindang = vars.Results.blindang;
WD_th = vars.Results.wall_th; %wall distance threshold
PLOT = vars.Results.PLOT;
PLOTfs = vars.Results.PLOTfs;
MOVIE = vars.Results.MOVIE;

SAVE = vars.Results.SAVE;
Fs = vars.Results.Fs;
dt = 1/Fs;
N = vars.Results.N;
USEINTEGRATOR = vars.Results.USEINTEGRATOR;
Nreps = vars.Results.num_rep;

r_m = vars.Results.r_m;
tau_m = vars.Results.tau_m;
S0 = vars.Results.S0;
Srest = vars.Results.Srest;
DeltaS = vars.Results.DeltaS;

% for saving
suffix_file_name = ['_with_integrator_',replay_file_to_load(1:end-4)];

load([loadpath,replay_file_to_load])
disp(['loaded ',replay_file_to_load]);

Nstim = size(x,1); % max number of stimuli fish
x = x.*arena_diam/2; % normalize stim coordinates
y = y.*arena_diam/2;



% bout rates from data:
% Translate to bout probability, take into account down time that fish can't bout
est_dec_in_rate = 0.55; %compensate for 'down times'
basic_bout_prob = ones(1,N)*all_bout_rates/Fs/est_dec_in_rate;  % here all fish have the same bout rates


ModelType = 'traj_replay';

vert_size = BL; % fish vertical size

% load recorded response function
load([loadpath,'5 fitted_response.mat']);
p1 = mdlp.Coefficients.Estimate;
p0 = 0.5;


% load wall interactions from data:
load([loadpath,'wall_response']);
wall_fun = @(x) aexp.*exp(bexp.*x); % note - the function is for BL not cm

% turning parameters:
Wall_sig = 45; % % sigma of normal distribution for turning angles from wall
Neigh_sig = 30; % sigma for normal distribution turning angles from neighbros


% create bout shape - using a  sigmoid function
% parameters for sigmoidal bouts
xsig = (0:7)/Fs;
x0 = 3.5/Fs;
maxV = avg_speed; % taken from data of 7 dpf (we may need to change for the others)
slope = 41; % set the slope (or add a distribution)

% speed in body length per second
Sup = maxV./(1+exp(-slope*(xsig-x0)));
Sdown = maxV./(1+exp(slope*(xsig-x0)));

% trnasform to cm per frame (so we can sum it)
S = ([Sup Sdown]*BL)/Fs;

Ttotal = length(T);  % define the simulation time

% define arena size and bounderies:
rad = arena_diam/2; % arena radius in cm
rad_start = rad; % radius to start agents in
boundx = sind(0:360)*rad;
boundy = cosd(0:360)*rad;



savepath = [loadpath,ModelType,suffix_file_name,'/'];

% define a 'tail' for plotting trajectories
if PLOT
    Tail_time = 2; % in sec
    Step = 15; % in frames
    Tail_length = round(Tail_time*Fs/(Step)); % in plotted points
    Tail = ceil(linspace(0,Tail_time,Tail_length)*Fs);
    % alpha values of the dots
    Alpha = linspace(0,0.5,length(Tail));Alpha = Alpha(end:-1:1);
    Cmap =LineMap;
    Cmap = [Cmap(1,:);repmat(Cmap(2,:),63,1)];
end


if MOVIE
    % make a movie file if needed
    obj = VideoWriter(MOVIE,'MPEG-4');
    obj.Quality = 100;
    obj.FrameRate = 30;
    open(obj);
end

for r = 1:Nreps %curr_rep
    a =tic;
    rng(r); % set random number for reproducability

    % file name for saving
    filename = ['sim_larva_rep',num2str(r),'_',ModelType];

    % create directory if needed
    if ~exist([savepath,filename],'dir')
        mkdir([savepath,filename]);
        cd([savepath,filename]);
    else
        cd([savepath,filename]);
    end

    %%
    % define state of each fish (0 - idle, 1 - in bout)
    State = zeros(N,Ttotal);

    % variables for position
    xfish = ones(N,Ttotal); % 1,1 values are to be ignored in analysis
    yfish = ones(N,Ttotal);

    % random starting positions:
    [x_strt,y_strt] = dotsInAcircle(1,0.9*(rad_start),(rad_start)/3,[0 0]);

    xfish(:,N) = x_strt;
    yfish(:,N) = y_strt;

    % variables for neighbor relative angle, orientation and distance
    relAng = zeros(N,Ttotal,Nstim);
    relOri = zeros(N,Ttotal,Nstim);
    relDist = zeros(N,Ttotal,Nstim);
    visualAngle = zeros(N,Ttotal,Nstim);
    internal_state_value = ones(N,Ttotal)*S0;
    input_to_integrator = ones(N,Ttotal);

    % fish body angle
    angle = zeros(N,Ttotal);

    % startig angle
    angle(:,1) = rand(N,1)*360; % set random angles for start (angle from north)

    % velocity vectors
    Vx = zeros(N,Ttotal);
    Vy = zeros(N,Ttotal);

    % speed
    Speed = zeros(N,Ttotal);

    % distance to wall
    wallD = zeros(N,Ttotal);


    % vector from center of arena to fish
    Norm = calculateNorm([xfish(:,1) yfish(:,1)]);

    % distance from wall:
    wallD(:,1) = (rad - Norm)/BL; % transform to body length since wall distance funciton in in BL

     % save the angle between fish orientation and the norm vec from arena
    % center to fish  values < 180 means wall is to the left (90 is perpendicular),
    % values > 180 means wall is to the right (270 is perpendicular)
    angFromWall = zeros(N,Ttotal);

    % normalized vector from center towards fish:
    vec_2_fish = [xfish(:,1)./Norm yfish(:,1)./Norm];

    % starting angle from wall 
    [~,ang_from_norm] = angOfVectors([sind(angle(:,1)) cosd(angle(:,1))],vec_2_fish);
    angFromWall(:,1) = ang_from_norm;

    Bouts = zeros(N,Ttotal); % var for bout times
    Bouts_end = zeros(N,Ttotal); % bout ends
    all_predicted_pright = zeros(N,Ttotal); % probabilites to turn right

    % variable for the obtained retina angles:
    all_ret_angles = zeros(N,Ttotal,2,N-1);


    % claculate relative properties for all fish
    t = 1;
    for f = 1:N
        [temp_dist,temp_angle,~] = relativeNeighborProp(f,[xfish(:,t);x(:,t)],[yfish(:,t);y(:,t)],angle(:,t));
        relDist(f,t,:) = temp_dist(2:end);
        relAng(f,t,:) = temp_angle(2:end);
        relOri(f,t,:) = ones(Nstim,1)*90;
        visualAngle(f,t,:) = 2*atand((vert_size/2)./relDist(f,t,:));
    end
    
    if PLOT
        figure('Units','Normalized','Position',[0.1 0.1 0.7 0.7])
    end

    for t = 2:Ttotal-15 % run simulation (lose time points at the end)


        % find all idle fish (currently not moving)
        idle = find(State(:,t)==0); % fish that are not moving

        for f = 1:length(idle) % loop over idle fish

            % get wall distance:
            tempWD = wallD(idle(f),t-1);

            % get relative distance, angle, and orientation of all neighboring fish
            tempD = squeeze(relDist(idle(f),t-1,:));
            neighD = tempD; % neighbor distnace
            tempA = squeeze(relAng(idle(f),t-1,:));
            neighA = tempA; % angle to neighbor
            neighA(neighA > 180) = neighA(neighA>180) - 360;
            tempO = squeeze(relOri(idle(f),t-1,:));
            neighO = tempO; % neighbor orientations

            % retinal size
            ret_angle = 2 * atand((vert_size/2)./neighD);

            % define prob of bout:
            p_bout = basic_bout_prob(idle(f));


            % flip a coin to decide on a bout
            bout_flag =  rand(1) < p_bout;
            Bouts(idle(f),t) = bout_flag; % save to variable

            % if a bout is starting - choose kinematics
            if bout_flag

                % Kinematics phase:
                % for fish switching to bouts, choose angle, and speed and duration
                % according to neighbors and walls.

                % turn away from wall as first priority

                if tempWD < WD_th

                    % check if fish will respond to the wall
                    pwall = rand(1) < wall_fun(tempWD)+0.5;

                    if pwall % if fish is turning from wall

                        % if wall is to the left (<180) turn right) and vice
                        p_right = angFromWall(idle(f),t-1) < 180;
                        ang = randn(1)*Wall_sig;

                        % choose turn direction according to p_right:
                        ang = abs(ang)*p_right-abs(ang)*(1-p_right);

                        ang = angle(idle(f),t-1) + ang;

                    end
                else
                    pwall = 0;
                end

                % if fish is responding to neighbors:
                if pwall==0 && abs(p1) > 0 


                    % calculate predicted response for each fish
                    % according to its size and side
                    left_right_side_fish = ones(size(neighA));
                    left_right_side_fish(neighA < 0) = -1; % neighbors to the left
                    left_side_ii = left_right_side_fish<0;

                    % remove angles of fish in the blind spot:
                    blind_spot_ii = neighA > (180-blindang) | neighA < (-180+blindang);
                    if sum(blind_spot_ii) > 0
                        left_right_side_fish(blind_spot_ii) = [];
                        left_side_ii(blind_spot_ii) = [];
                        ret_angle(blind_spot_ii) = [];
                    end

                    % save to var:
                    all_ret_angles(idle(f),t,1,left_side_ii) = ret_angle(left_side_ii);
                    all_ret_angles(idle(f),t,2,~left_side_ii) = ret_angle(~left_side_ii);
                    temp_right =  ret_angle(~left_side_ii);
                    temp_left = ret_angle(left_side_ii);
                    
                    % use the avergae squared veritcal size
                    left_response = nansum(temp_left.*...
                        temp_left)./nansum(temp_left);

                    right_response = nansum(temp_right.*...
                        temp_right)./nansum(temp_right);

                    % take the singed sum between the eyes
                    Diff = nansum([-right_response,left_response]);


                    % use the integrator potenatial to change
                    % the fish response function.

                    % we liearly scale the state value such that for S0 p1' = p1.
                    % and for S0+DeltaS will result in p1' = 0;
                    curr_state_value = -(internal_state_value(t-1)-S0)/DeltaS; 
      

                    % update slope according to integrator
                    if USEINTEGRATOR
                        curr_slope = p1+p1*curr_state_value;
                    else
                        curr_slope = p1;
                    end
                    % make sure it non negative (no attraction)

                    % curr_slope = p1;    % for the case of no modulation
                    curr_slope(curr_slope<0) = 0;


                    predicted_pright = Diff*curr_slope+p0;
                    all_predicted_pright(idle(f),t) = predicted_pright;


                    % calculate the prob to turn left or right:
                    p_right = rand(1) < predicted_pright;

                    ang = randn(1)*Neigh_sig; % angle from a dist with sigma = 30;

                    % choose turn direction according to p_right:
                    % positive is a right turn negative is a left turn
                    ang = abs(ang)*p_right-abs(ang)*(1-p_right);

                   
                    % change fish angle of swimming
                    ang = angle(idle(f),t-1) + ang;

                elseif pwall==0 % if there is no social or wall interactions

                    new_turn = abs(randn(1)*Wall_sig); % angles no sign
                    
                    % sign of turn, -1 or 1 randomly
                    new_turn_side = (double((rand < 0.5))-0.5)/0.5;
                   
                    % set the turn and its direction
                    ang = angle(idle(f),t-1)+new_turn*new_turn_side;

                end
                % set kinematics

                % speed in body length per second
                Sup = maxV./(1+exp(-slope*(xsig-x0)));
                Sdown = maxV./(1+exp(slope*(xsig-x0)));

                % trnasform to cm per frame (so we can sum it)
                S = ([Sup Sdown]*BL)/Fs;
                % velocity
                vx = S*sind(ang);
                vy = S*cosd(ang);

                % position
                xx = cumsum(vx);
                yy = cumsum(vy);

                % check if the expected trajectory will go outside the arena
                dd = calculateNorm([(xx+xfish(idle(f),t-1))' (yy+yfish(idle(f),t-1))']); % radius

                outside = dd > rad; % position outside
                % if needed correct to stop at the wall
                if sum(outside)>0
                    vx(outside) = 0;
                    vy(outside) = 0;
                    xx = cumsum(vx);
                    yy = cumsum(vy);
                    S(outside) = 0;
                end

                % Update phase:
                % update all positions and kinematic variables of moving fish
                State(idle(f),t:t+length(S)-2) = ~outside(1:end-1); % set bout times
                State(idle(f),t:t+length(S)-1) = 1;
                angle(idle(f),t:t+length(S)-1) = ang;
                Bouts_end(idle(f),t+sum(~outside)-1) = 1;

                Speed(idle(f),t:t+length(S)-1) = S;
                Vx(idle(f),t:t+length(vx)-1) = vx;
                Vy(idle(f),t:t+length(vy)-1) = vy;

                xfish(idle(f),t:t+length(xx)-1) = xx+xfish(idle(f),t-1);
                yfish(idle(f),t:t+length(xx)-1) = yy+yfish(idle(f),t-1);


            else % if no bout - keep same position and angle
                xfish(idle(f),t) = xfish(idle(f),t-1);
                yfish(idle(f),t) = yfish(idle(f),t-1);
                angle(idle(f),t) = angle(idle(f),t-1);
            end

        end

        % calculate relative position and angle for each fish
        for f = 1:N
            [temp_dist,temp_angle,~] = relativeNeighborProp(f,[xfish(:,t);x(:,t)],[yfish(:,t);y(:,t)],angle(:,t));
            relDist(f,t,:) = temp_dist(2:end);
            relAng(f,t,:) = temp_angle(2:end);
            relOri(f,t,:) = ones(Nstim,1)*90;
            visualAngle(f,t,:) = 2*atand((vert_size/2)./relDist(f,t,:));

            
            if t>2 % start updating intrnal state variable
                max_va = max(visualAngle(f,t-2:t,:),[],3); % maximum visual angles
                d_va = diff(max_va,1,2); % find positive changes
                d_va(d_va<0 | isnan(d_va)) = 0; 
                input_to_integrator(f,t) = d_va(2); % use as input

                %calc the change in internal state variable
                dv = (Srest-internal_state_value(t-1) + r_m*d_va(2))/tau_m; 
                
                % update internal state
                internal_state_value(f,t) = internal_state_value(f,t-1)+dv;
            end
        end


        % wall distance
        wallD(:,t) = (rad - calculateNorm([xfish(:,t) yfish(:,t)]))/BL;

        % vector from center of arena to fish
        Norm = calculateNorm([xfish(:,t) yfish(:,t)]);
        vec_2_fish = [xfish(:,t)./Norm yfish(:,t)./Norm];

        % angle from wall
        [~,angFromWall(:,t)] = angOfVectors([sind(angle(:,t)) cosd(angle(:,t))],vec_2_fish);

        % plotting
        if PLOT && t>max(Tail) && mod(t,PLOTfs)==0
            subplot(3,4,[1:3,5:7,9:11])
            
            % neighbors have a 'shadow' in time
            alpha_p = 1;
            trail_back = t - Tail;
            trail_back = trail_back(trail_back>0);
            tempx = [xfish(1,trail_back);x(:,trail_back)];
            tempy = [yfish(1,trail_back);y(:,trail_back)];
            for cc = 1:length(trail_back)
                scatter(tempx(:,cc),tempy(:,cc),15,Cmap(1:size(tempx,1),:),'filled',...
                    'MarkerEdgeColor','none','MarkerFaceAlpha',Alpha(alpha_p));
                hold on;
                alpha_p = alpha_p + 1;
            end
            
            % arena walls
            plot(boundx,boundy,'k');
            hold off
            axis image;
           
            title(round(t/Fs,3));
            axis off;

            % input to integrator
            subplot(3,4,4)
            plot(T(1:t),input_to_integrator(1:t));
          
            ylim([0 40]);
            xlim([0 T(t)+0.1])
            ylabel({'Max increase in',' occupancy [deg]'})
            xlabel('Time [s]')
            box off;
            set(gca,'FontSize',15)

            % internal state variable
            subplot(3,4,8)
 
            plot(T(1:t),internal_state_value(1:t)); hold on;
            plot([0, T(t)+0.1],[S0 S0],'k--','LineWidth',1)
            ylim([S0-0.2,S0+0.2]);
            xlim([0 T(t)+0.1])
            ylabel('S(t)')
            xlabel('Time [s]')
            set(gca,'FontSize',15)
            box off;
            hold off;
            if MOVIE
                f = getframe(gcf);
                writeVideo(obj,f);
            end
            pause(0.01);
        end


        if mod(t,10000)==0
            disp(t)
        end


    end


    groupname = filename;
    % save data
    toc(a)
    if SAVE


        save('coor','xfish','yfish','boundx','boundy');
        save('group_prop','relAng','relOri','relDist','visualAngle',...
            'angle','Vx','Vy','Speed','wallD','angFromWall',...
            'Bouts','Bouts_end','all_ret_angles','all_predicted_pright','internal_state_value',...
            'input_to_integrator');

        avg_BL = BL;
        save('Header','Fs','avg_BL','N','groupname','arena_diam');
        disp('data saved')



    else
        disp('data not saved')
    end


end

if MOVIE
    close(obj);
end

if SAVE
cd(savepath)
save('model_params','r_m','tau_m','S0','DeltaS')


% save a version of the code we run
FileNameAndLocation=[mfilename('fullpath')];
newbackup=[savepath,'simFile.m'];
currentfile=strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);

% consolidata all fish simulated
collectTrajReplaySimulationData('Tframes',length(xfish)*dt,'Folder',savepath);
disp('files consolidated together')
cd ..
end
