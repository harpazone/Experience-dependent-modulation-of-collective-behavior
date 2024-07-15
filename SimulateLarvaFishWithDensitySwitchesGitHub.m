function SimulateLarvaFishWithDensitySwitchesGitHub(loadpath,Norder,switch_times, varargin)
%Description - simulating a group of agents interacting based on visual
%occucpancy intractions. Group sizes switches according to Norder at switch
%times. Interaction strength is modulated by an internal state values. 


% set defaults
defaultFs = 50; % sampling rate Hz

defaultbout_rate = 1.65; %
defaultSpeed = 1.5; % BL/s
defaultBL = 0.4; % fish body length cm;
defaultArena_diam = 6.5; % arena diameter in cm
defaultOne_side_blindang = 15;  % in degrees s.t. 2*ang is the full blind angle
default_num_rep = 1; % number of repetitions
defaultPLOT = 1; % plot data yes/no
defaultPLOTfs = 50; % frame rate of ploted data (compare to Fs of sim)
defualtSAVE = 1; %save data yes/no
defualtUSEINTEGRATOR = 1; %simulate with/without integrator
defualtMOVIE = false; %make a movie - false or text of movie name


defualt_r_m = 13; % input scaling factor
defualt_tau_m = 360000*3; % leak time constant in seconds
defualt_S0 = 0.6; % starting point
defualt_Srest = 0; %resting value
defualt_DeltaS = 0.19; % S(t) = S0+DSresttaS gives p1' = 0; (no repsone to neighbors)

% parse
vars = inputParser;
addParameter(vars,'Fs',defaultFs);

addParameter(vars,'PLOT',defaultPLOT);
addParameter(vars,'PLOTfs',defaultPLOTfs);

addParameter(vars,'SAVE',defualtSAVE);
addParameter(vars,'all_bout_rate',defaultbout_rate);
addParameter(vars,'all_avg_speed',defaultSpeed);
addParameter(vars,'all_BL',defaultBL);
addParameter(vars,'all_arean_rad',defaultArena_diam);
addParameter(vars,'blindang',defaultOne_side_blindang);
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
PLOT = vars.Results.PLOT;
PLOTfs = vars.Results.PLOTfs;
MOVIE = vars.Results.MOVIE;

SAVE = vars.Results.SAVE;
Fs = vars.Results.Fs;
dt = 1/Fs;
Nmin = min(Norder); % low density
Nmax = max(Norder); % high density

USEINTEGRATOR = vars.Results.USEINTEGRATOR;
Nreps = vars.Results.num_rep;

r_m = vars.Results.r_m;
tau_m = vars.Results.tau_m;
S0 = vars.Results.S0;
Srest = vars.Results.Srest;
DeltaS = vars.Results.DeltaS;


% add order and times to group name
folder_name = 'inteacting_group';
for i = 1:length(Norder)
    folder_name = [folder_name,'_',num2str(Norder(i))];
end

folder_name = [folder_name,'_swith_times'];

for i = 1:length(Norder)
    folder_name = [folder_name,'_',num2str(switch_times(i)/60/Fs)];
end

suffix_file_name = ['_with_integrator_',folder_name];



% bout rates from data:
% Translate to bout probability, take into account down time that fish
% 'can't' bout
est_dec_in_rate = 0.55; %compensate for 'down times' when fish cant bout (found from data)
basic_bout_prob = all_bout_rates/Fs/est_dec_in_rate*ones(1,Nmax);

% model type:
ModelType = 'group_interactions';

% model parameters
vert_size = BL/2; % vertical height
load([loadpath,'turn_7dpf_VR.mat']);
p1 = FIT.p1;

p0 = FIT.p2;


load([loadpath,'wall_fun_groups']);
wall_fun_min = wall_fun5;
wall_fun_max = wall_fun20; % different wall functions for the group sizes
% as seen in data. Also note - the function is for BL not cm

Cmap = LineMap;


% turning parameters:
Wall_sig = 45; % sigma of normal distribution for turning angles from wall
Neigh_sig = 30; % sigma for normal distribution turning angles from neighbros


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

Ttotal = switch_times(end); % define the simulation time

% define arena size and bounderies:
rad = arena_diam/2; % arena radius in cm
rad_start = rad; % radius to start agents in
boundx = sind(0:360)*rad;
boundy = cosd(0:360)*rad;

switch_times(end) = switch_times(end)-15; % remove 15 frames of simulation not
% to exceed size of variables


% path to save simulation
savepath = [loadpath,ModelType,suffix_file_name,'/'];

% create a movie if needed
if MOVIE
obj = VideoWriter(MOVIE,'MPEG-4');
obj.Quality = 100;
obj.FrameRate = 30;
open(obj);
end


% define a 'tail' for plotting trajectories
if PLOT
    Tail_time = 2; % in sec
    Step = 10; % in frames
    Tail_length = Tail_time*Fs/(Step); % in plotted points
    Tail = ceil(linspace(0,Tail_time,Tail_length)*Fs);
    % alpha values of the dots
    Alpha = linspace(0,0.5,length(Tail));Alpha = Alpha(end:-1:1);
    Cmap =distinguishable_colors(Nmax);
    colors = [Cmap, ones(Nmax,1)*max(Alpha)];
end


for r = 1: Nreps
    a =tic;
    rng(r); % set random number for reproducability
    
    % filename for saving
    filename = ['sim_group_rep',num2str(r)];


    % create directory if needed
    if ~exist([savepath,filename],'dir')
        mkdir([savepath,filename]);
        cd([savepath,filename]);
    else
        cd([savepath,filename]);
    end

    %%
    % define state of each fish (0 - idle, 1 - in bout)
    State = zeros(Nmax,Ttotal);

    % variables for position
    xfish = ones(Nmax,Ttotal)*nan; % nan values are to be ignored in simulations
    yfish = ones(Nmax,Ttotal)*nan;


    % variables for neighbor relative angle, orientation and distance
    relAng = zeros(Nmax,Ttotal,Nmax-1);
    relOri = zeros(Nmax,Ttotal,Nmax-1);
    relDist = zeros(Nmax,Ttotal,Nmax-1);
    visualAngle = zeros(Nmax,Ttotal,Nmax-1);
    internal_state_value= ones(Nmax,Ttotal)*S0;
    input_to_integrator = ones(Nmax,Ttotal)*nan;

    % fish body angle
    angle = ones(Nmax,Ttotal)*nan;


    % velocity vectors
    Vx = zeros(Nmax,Ttotal);
    Vy = zeros(Nmax,Ttotal);

    % speed
    Speed = zeros(Nmax,Ttotal);

    % distance to wall
    wallD = zeros(Nmax,Ttotal);

    % save the angle between fish orientation and the norm vec from arena
    % center to fish  values < 180 means wall is to the left (90 is perpendicular),
    % values > 180 means wall is to the right (270 is perpendicular)
    angFromWall = zeros(Nmax,Ttotal);


    Bouts = zeros(Nmax,Ttotal); % center of bout times
    Bouts_end = zeros(Nmax,Ttotal); % end of bout times
    all_p = zeros(Nmax,Ttotal); % probabilities to turn right/left 
    
    % SHAM variable unused. Included so that random number generator will
    % explicitly match ones used in simulation in Harpaz et al 2024.
    bout_side_flag = ones(Nmax,1);
    bout_side_flag(rand(Nmax,1) < 0.5) = -1;

    % variable for the obtained retina angles:
    all_ret_angles = zeros(Nmax,Ttotal,2,Nmax-1);

    if PLOT
        figure('Units','normalized','Position',[0.1 0.1 0.7 0.7])
    end
    p_t = 1; % counter for actual positions

    % loop over experiment parts
    for phase_id = 1:length(Norder)
        Ncurr = Norder(phase_id);
        % random starting positions:
        if phase_id == 1 % experiment starts
            [x_strt,y_strt] = dotsInAcircle(Ncurr,0.9*(rad_start),(rad_start)/3,[0 0]);

            xfish(1:Ncurr,p_t) = x_strt;
            yfish(1:Ncurr,p_t) = y_strt;
            % startig angle
            angle(1:Ncurr,p_t) = rand(Ncurr,1)*360; % set random angles for start (angle is angle from north)

        else
            % if this is an advanced phases take previous positions
            internal_state_value(1:Ncurr,p_t) = internal_state_value(1:Ncurr,p_t-1);
            xfish(1:Ncurr,p_t) = xfish(1:Ncurr,p_t-1);
            yfish(1:Ncurr,p_t) = yfish(1:Ncurr,p_t-1);
            angle(1:Ncurr,p_t) = angle(1:Ncurr,p_t-1);

            if Ncurr == Nmax % if we moved to higher N and we are missing fish addd them
                % add missing fish
                [x_strt,y_strt] = dotsInAcircle(Nmax-Nmin,0.9*(rad_start),(rad_start)/3,[0 0]);

                xfish(Nmin+1:Nmax,p_t) = x_strt;
                yfish(Nmin+1:Nmax,p_t) = y_strt;
                angle(Nmin+1:Nmax,p_t) = rand(Nmax-Nmin,1)*360; %

                % use the average and std of the internal state value to
                % assign these value to newly added agents
                mean_p = mean(internal_state_value(1:Nmin,p_t-1));
                std_p = std(internal_state_value(1:Nmin,p_t-1));
                internal_state_value(Nmin+1:Nmax,p_t) = randn(Nmax-Nmin,1)*std_p+mean_p;
            else
                % if we moved to low N
                % remove uneeded fish
                xfish(Nmin+1:Nmax,p_t) = nan;
                yfish(Nmin+1:Nmax,p_t) = nan;
                internal_state_value(Nmin+1:Nmax,p_t) = nan;

            end
        end

        % vector from center of arena to fish
        Norm = calculateNorm([xfish(:,1) yfish(:,1)]);
        % distance from wall:
        wallD(:,p_t) = (rad - Norm)/BL; % transform to body length since wall distance funciton in in BL
        % normalized vector from center towards fish:
        vec_2_fish = [xfish(:,p_t)./Norm yfish(:,p_t)./Norm];

        [~,ang_from_norm] = angOfVectors([sind(angle(:,p_t)) cosd(angle(:,p_t))],vec_2_fish);
        angFromWall(:,p_t) = ang_from_norm;

        % claculate relative propertied for all fish
        for f = 1:Ncurr
            [temp_dist,temp_angle,~] = relativeNeighborProp(f,xfish(:,p_t),yfish(:,p_t),angle(:,p_t));
            temp_dist(f) = []; % remove self
            temp_angle(f) = [];
            relDist(f,p_t,:) = temp_dist;
            relAng(f,p_t,:) = temp_angle;
            relOri(f,p_t,:) = ones(Nmax-1,1)*90;
            visualAngle(f,p_t,:) = 2*atand((vert_size/2)./relDist(f,p_t,:));
        end

        p_t = p_t+1;

        % loop over times in the current experiment phase
        for t = 2: switch_times(phase_id+1)-switch_times(phase_id) % for all points


            % State phase:
            % flip a coin for the bout state

            idle = find(State(1:Ncurr,p_t)==0); % fish that are not moving

            for f = 1:length(idle) % loop over idle fish

                % get wall distance:
                tempWD = wallD(idle(f),p_t-1);

                % get relative distance, angle, and orientation of all neighboring
                % fish to estimate retina angle:

                tempD = squeeze(relDist(idle(f),p_t-1,:));
                neighD = tempD; % neighbor distnace
                % neighD(idle(f)) = [];
                tempA = squeeze(relAng(idle(f),p_t-1,:));
                neighA = tempA; % angle to neighbor
                % neighA(idle(f)) = [];
                neighA(neighA > 180) = neighA(neighA>180) - 360;

                tempO = squeeze(relOri(idle(f),p_t-1,:));
                neighO = tempO; % neighbor orientations
                % neighO(idle(f)) = [];
                ret_angle = 2 * atand((vert_size/2)./neighD);


                % get prob of bout:
                p_bout = basic_bout_prob(idle(f));

                % flip a coin to decide on a bout

                bout_flag =  rand(1) < p_bout;
                Bouts(idle(f),p_t) = bout_flag;
                % if a bout is starting - choose kinematics
                if bout_flag

                    % Kinematics phase:
                    % for fish switching to bouts, choose angle, and speed and duration
                    % according to neighbors and walls.


                    % turn away from wall as first priority
                    % check if fish will respond to the wall

                    if Ncurr==Nmin
                        pwall = rand(1) < wall_fun_min(tempWD); % low density

                    else
                        pwall = rand(1) < wall_fun_max(tempWD); % high density
                    end

                    if pwall % if fish is turning from wall

                        % if wall is to the left (<180) turn right) and vice
                        p_right = angFromWall(idle(f),p_t-1) < 180;
                        ang = randn(1)*Wall_sig;

                        % choose turn direction according to p_right:
                        ang = abs(ang)*p_right-abs(ang)*(1-p_right);

                        ang = angle(idle(f),p_t-1) + ang;

                    end


                    % if fish is responding to neighbors:
                    if p_t > 1 && pwall==0 && abs(p1)>0



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
                        all_ret_angles(idle(f),p_t,1,left_side_ii) = ret_angle(left_side_ii);
                        all_ret_angles(idle(f),p_t,2,~left_side_ii) = ret_angle(~left_side_ii);

                        % get normalized potential to adjust response slope
                        curr_potential = -(internal_state_value(idle(f),p_t-1)-S0)/DeltaS;

                        if USEINTEGRATOR
                            curr_slope = p1+p1*curr_potential; % adjust slope
                        else
                            curr_slope = p1; % no modulation
                        end
                        
                        % calculate response at each angle
                        temp_pright = (ret_angle.*left_right_side_fish*-1)*curr_slope+p0;


                        % take a weighted avg on each side (reponse
                        % tendency is the difference from 0.5)
                        left_response = nansum(temp_pright(left_side_ii).*...
                            ret_angle(left_side_ii))./nansum(ret_angle(left_side_ii))-0.5;

                        right_response = nansum(temp_pright(~left_side_ii).*...
                            ret_angle(~left_side_ii))./nansum(ret_angle(~left_side_ii))-0.5;


                        % sum up the signed responses on both sides:
                        predicted_pright = nansum([left_response,right_response])+0.5;
                        all_p(idle(f),p_t) = predicted_pright;



                        % calculate the prob to turn left or right:
                        p_right = rand(1) < predicted_pright;

                        ang = randn(1)*Neigh_sig; % angle from a normal distribution

                        % choose turn direction according to p_right:
                        % positive is a right turn negative is a left turn
                        ang = abs(ang)*p_right-abs(ang)*(1-p_right);


                        % change angle of swimming
                        ang = angle(idle(f),p_t-1) + ang;

                    elseif pwall==0 % if there is no interaction of any sort

                        new_turn = abs(randn(1)*Wall_sig); % angles no sign

                        % sign of turn, -1 or 1 randomly
                        new_turn_side = (double((rand < 0.5))-0.5)/0.5;

                        % set new swimming direction
                        ang = angle(idle(f),p_t-1)+new_turn*new_turn_side;

                    end
                    % set kinematics
                    S = ([Sup Sdown]*BL)/Fs;

                    % velocity
                    vx = S*sind(ang);
                    vy = S*cosd(ang);

                    % position
                    xx = cumsum(vx);
                    yy = cumsum(vy);

                    % check if the expected trajectory will go outside the arena
                    dd = calculateNorm([(xx+xfish(idle(f),p_t-1))' (yy+yfish(idle(f),p_t-1))']); % radius

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
                    State(idle(f),p_t:p_t+length(S)-2) = ~outside(1:end-1); % set bout times
                    State(idle(f),p_t:p_t+length(S)-1) = 1;
                    angle(idle(f),p_t:p_t+length(S)-1) = ang;
                    Bouts_end(idle(f),p_t+sum(~outside)-1) = 1;

                    Speed(idle(f),p_t:p_t+length(S)-1) = S;
                    Vx(idle(f),p_t:p_t+length(vx)-1) = vx;
                    Vy(idle(f),p_t:p_t+length(vy)-1) = vy;

                    xfish(idle(f),p_t:p_t+length(xx)-1) = xx+xfish(idle(f),p_t-1);
                    yfish(idle(f),p_t:p_t+length(xx)-1) = yy+yfish(idle(f),p_t-1);

                else
                    xfish(idle(f),p_t) = xfish(idle(f),p_t-1);
                    yfish(idle(f),p_t) = yfish(idle(f),p_t-1);
                    angle(idle(f),p_t) = angle(idle(f),p_t-1);
                end


            end

            % calculate relative position and angle for each fish
            for f = 1:Ncurr
                [temp_dist,temp_angle,~] = relativeNeighborProp(f,xfish(:,p_t),yfish(:,p_t),angle(:,p_t));
                temp_dist(f) = []; % remove self
                temp_angle(f) = [];
                relDist(f,p_t,:) = temp_dist;
                relAng(f,p_t,:) = temp_angle;
                relOri(f,p_t,:) = ones(Nmax-1,1)*90;
                visualAngle(f,p_t,:) = 2*atand((vert_size/2)./relDist(f,p_t,:));

                if  p_t>2 % update integrator module
                    max_va = max(visualAngle(f,p_t-2:p_t,:),[],3); % maximum visual angle
                    d_va = diff(max_va,1,2); % get positive changes in maximum visual angle
                    d_va(d_va<0 | isnan(d_va)) = 0;
                    
                    % calculate change in integrator states
                    input_to_integrator(f,p_t) = d_va(2); 
                    dv = (Srest-internal_state_value(f,p_t-1) + r_m*d_va(2))/tau_m; %calc the change in voltage
                    % update internal state of integrator
                    internal_state_value(f,p_t) = internal_state_value(f,p_t-1)+dv;
                end
            end

            % calc wall distance
            wallD(:,p_t) = (rad - calculateNorm([xfish(:,p_t) yfish(:,p_t)]))/BL;

            % vector from center of arena to fish
            Norm = calculateNorm([xfish(:,p_t) yfish(:,p_t)]);
            vec_2_fish = [xfish(:,p_t)./Norm yfish(:,p_t)./Norm];

            % angle from wall
            [~,angFromWall(:,p_t)] = angOfVectors([sind(angle(:,p_t)) cosd(angle(:,p_t))],vec_2_fish);

            % plot if needed
            if PLOT && p_t>max(Tail) && mod(p_t,PLOTfs)==0
                subplot(3,4,[1:3,5:7,9:11])
                % plot the dots 'shadow' in the past
                alpha_p = 1;
                trail_back = p_t - Tail;
                trail_back = trail_back(trail_back>0);
                for cc = trail_back
                    scatter(xfish(1:Ncurr,cc)',yfish(1:Ncurr,cc)',15,Cmap(1:Ncurr,:),'filled',...
                        'MarkerEdgeColor','none','MarkerFaceAlpha',Alpha(alpha_p));
                    hold on;
                    alpha_p = alpha_p + 1;
                end

                % add arena
                plot(boundx,boundy,'k');
                hold off
                axis image;

                title(round(p_t/Fs,3));
                axis off;

                % input to integrator
                subplot(3,4,4)
                h = plot((1:p_t)/Fs,input_to_integrator(1:Ncurr,1:p_t)');
                set(h,{'Color'},num2cell(colors(1:Ncurr,:),2));
                ylim([0 40]);
                xlim([0 p_t/Fs+0.1])
                ylabel({'Max increase in ','occupancy [deg]'})
                xlabel('Time [s]')
                box off;

                set(gca,'FontSize',15)

                % internal state variable
                subplot(3,4,8)

                h1 = plot((1:p_t)/Fs,internal_state_value(1:Ncurr,1:p_t));
                hold on;
                plot([0 p_t/Fs+0.1],[0.6 0.6],'k--','LineWidth',1);
                set(h1,{'Color'},num2cell(colors(1:Ncurr,:),2));
                ylim([S0-0.2,S0+0.2]);
                xlim([0 p_t/Fs+0.1])
                ylabel('S(t)')
                xlabel('Time [s]')
                box off;
                hold off;
                set(gca,'FontSize',15)
                
                if MOVIE % add to movie if needed

                f = getframe(gcf);
                writeVideo(obj,f);
                end
                pause(0.01);
            end


            if mod(p_t,10000)==0
                disp(p_t)
            end
            p_t = p_t+1;

        end
    end
end
groupname = filename;
% save data
toc(a)

if MOVIE
    close(obj);
end

if SAVE

    save('coor','xfish','yfish','boundx','boundy');
    save('group_prop','relAng','relOri','relDist','visualAngle',...
        'angle','Vx','Vy','Speed','wallD','angFromWall',...
        'Bouts','Bouts_end','all_ret_angles','all_p','internal_state_value',...
        'input_to_integrator');

    avg_BL = BL;
    save('Header','Fs','avg_BL','Nmin','Nmax','groupname','arena_diam');
    disp('data saved')


else
    disp('data not saved')
end


cd(savepath)
save('model_params','r_m','tau_m','S0')


% save a version of the code we run
FileNameAndLocation=[mfilename('fullpath')];
newbackup=[savepath,'simFile.m'];
currentfile=strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);

cd ..






