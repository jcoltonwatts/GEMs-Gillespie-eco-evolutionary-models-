function [x_stand, pop_stand, stand_times, pop_data_out, x_data_out, x_var_data_out] = GEMv4_bdCompModel_bdTradeoff(params, state_parameter_match, cv_vector, h2_vector, num_replicates, y0, t_max)

% this is a GEM for an alternative formulation of logistic growth
% call this function to run it

%% standardized time steps for storing time series
% import the time span (t_max) but decide on step lengths
stand_times = 0:0.25:t_max;
num_time_steps = length(stand_times); % calculate the number of standardized time steps

%% preallocate space in which to log data at standardized times
% need one for each population, one for each evolving trait, and one for
% the variance in each evolving trait
no_species = length(y0); % the number of populations in the model = number of starting sizes
no_params = length(params);
no_columns = length(params)+1; % user must update this

pop_stand = zeros(no_species,num_time_steps,num_replicates); % population size
x_stand = nan(no_params,num_time_steps,num_replicates); % trait
x_var_stand = nan(no_params,num_time_steps,num_replicates); % trait variance

%% run for loop for each replicate simulation
parfor i = 1:num_replicates % for parallel computing
%for i = 1:num_replicates % not parallel
    t = 0; % initial time
    rng('shuffle'); % change random number seed    
    i % display replicate in the Command Window
    % pre-allocate for the sliced standardized variables
        blank_comm_matrix = nan(sum(y0),no_columns); % set up a NaN matrix for all individuals and traits
        pop_slice = zeros(no_species,num_time_steps);
        x_slice = nan(length(params),num_time_steps);
        x_var_slice = nan(length(params),num_time_steps);    

    % assign initial population sizes to sliced abundance variables    
        pop_slice(:,1) = y0';
            
    % pull initial distribution of trait (repeat for all populations)
    end_row = cumsum(y0);
    starting_row = [1 1+end_row(1:length(end_row)-1)];
    % determine how many of each state to pull
    ind_to_assign = state_parameter_match.*y0';
    
    for qq = 1:length(y0) % loop through states
        blank_comm_matrix(starting_row(qq):end_row(qq),1) = qq;
        for zz = 1:length(params) % loop through parameters   
            temp  = V4_pick_individuals(params(zz),cv_vector(zz)*params(zz),ind_to_assign(qq,zz));
            if isempty(temp) == 0
                blank_comm_matrix(starting_row(qq):end_row(qq),1+zz) = temp;
%                 if qq == 3 && zz > 24 % special case of norm distributed parameters
%                     blank_comm_matrix(starting_row(qq):end_row(qq),1+zz) = normrnd(0,1,ind_to_assign(qq,zz),1);
%                 end     
            end
        end
        % now overwrite values for d_min1 (param 5) and d_min2 (param 7)
        % according to b-d tradeoff. That is, make them a function of b_max1
        % (param 1) and b_max2 (param 3), respectively, depending on which
        % state (qq) is currently being entered in blank_comm_matrix
        if qq == 1
            blank_comm_matrix(starting_row(qq):end_row(qq), 5 + 1) = 2 * (blank_comm_matrix(starting_row(qq):end_row(qq), 1 + 1)).^2
        end
        if qq == 2
            blank_comm_matrix(starting_row(qq):end_row(qq), 7 + 1) = 2 * (blank_comm_matrix(starting_row(qq):end_row(qq), 3 + 1)).^2
        end
        % Check with a plot
        % figure(999)
        % plot(blank_comm_matrix(blank_comm_matrix(:,1) == 1, 5+1), blank_comm_matrix(blank_comm_matrix(:,1) == 1, 1+1), '-', 'Color', [1 0.77 0.3], 'LineWidth', 0.5)
    end
    
    x_dist_init = blank_comm_matrix;
    x_slice(1:no_columns-1,1) = nanmean(x_dist_init(:,2:no_columns),1)'; % initial mean trait vector
    x_var_slice(:,1) = nanvar(x_dist_init(:,2:no_columns),1)'; % initial variances in trait

    %% Initiate core GEM algorithm
    time_step_index = 2; % start a counter for standard times
    time_step = stand_times(time_step_index); % assign first standard time step

    % assign initial pop sizes to reduction variables (updated abundances)
        R1 = y0(1);
        R2 = y0(2);        
        x_dist = x_dist_init;
        
    while t < t_max
        % loop through species to find individuals for each state
        params_next = nan(no_species,no_params);
        whosnext = nan(size(y0));
        for zz = 1:no_species
            inds_in_state = find(x_dist(:,1)==zz);
            which_params = find(state_parameter_match(zz,:));
            if isempty(inds_in_state) == 0
                which_row = randi(length(inds_in_state));
                whosnext(zz) = inds_in_state(which_row);
                params_next(zz,which_params) = x_dist(whosnext(zz),1+which_params);
            else
                params_next(zz,which_params) = 0; % if pop is gone, set parameters to 0
            end
        end
        
        % pull out and re-assign parameters with names
        % state 1 parameters (R1) [1,2,5,6,9,11]
        b_max1 = params_next(1,1);
        b11 = params_next(1,2);
        b12 = params_next(1,9);
        d_min1 = params_next(1,5);
        d11 = params_next(1,6);
        d12 = params_next(1, 11);
        % state 2 parameters (R2) [3,4,7,8,10,12]
        b_max2 = params_next(2,3);
        b22 = params_next(2,4);
        b21 = params_next(2, 10);
        d_min2 = params_next(2, 7);
        d22 = params_next(2, 8);
        d21 = params_next(2, 12);

        % set up rates of each possible event
        % 1 birth resource 1
            b_R1 = max((b_max1 - b11*R1 - b12*R2),0)*R1;
        % 2 death resource 1
            d_R1 = max((d_min1 + d11*R1 + d12*R2), 0)*R1;
        % 3 birth resource 2
            b_R2 = max((b_max2 - b22*R2 - b21*R1),0)*R2;
        % 4 death resource 2
            d_R2 = max((d_min2 + d22*R2 + d21*R1), 0)*R2;
            
    % sum the events to make wheel of fortune
        CS_vector = cumsum([b_R1 d_R1 b_R2 d_R2]);
        Slice_widths = CS_vector./CS_vector(end);
        LI = rand < Slice_widths;
        Event_index = find(LI,1,'first');
        
        if Event_index == 1 % choose birth of state 1 (R1)
            if h2_vector(1) == 0
                x_parent = (1-h2_vector).*nanmean(x_dist_init(find(x_dist_init(:,1)==1),2:no_columns)) + h2_vector.*x_dist(whosnext(1),2:no_columns);
            elseif h2_vector(1) > 0
                x_parent = (1-h2_vector).*nanmean(x_dist(find(x_dist(:,1)==1),2:no_columns)) + h2_vector.*x_dist(whosnext(1),2:no_columns);
            end
            off_std = sqrt(1-h2_vector.^2).*((1-h2_vector).*nanstd(x_dist_init(find(x_dist_init(:,1)==1),2:no_columns))+h2_vector.*nanstd(x_dist(find(x_dist(:,1)==1),2:no_columns))); % offspring trait distribution std
            
            % x_dist(size(x_dist,1)+1,2:no_columns) = V4_pick_individuals(x_parent,off_std,1); % return trait
            x_dist(size(x_dist,1)+1,2:no_columns) = V4_pick_individuals(x_parent,off_std,1); % return trait -- this one increases variance

            x_dist(end,1) = 1;
            % tradeoff! -- overwrite d_min1 value (param 5) according to b_max1 (param 1)
            x_dist(end,5+1) = 2 * x_dist(end,1+1)^2;

        elseif Event_index == 2 % choose death of state 1 (R1)
            x_dist(whosnext(1),:) = []; % reduce dist by lost individual
            
        elseif Event_index == 3 % choose birth of state 2 (R2)
            if h2_vector(3) == 0
                x_parent = (1-h2_vector).*nanmean(x_dist_init(find(x_dist_init(:,1)==2),2:no_columns)) + h2_vector.*x_dist(whosnext(2),2:no_columns);
            elseif h2_vector(3) > 0
                x_parent = (1-h2_vector).*nanmean(x_dist(find(x_dist(:,1)==2),2:no_columns)) + h2_vector.*x_dist(whosnext(2),2:no_columns);
            end
            off_std = sqrt(1-h2_vector.^2).*((1-h2_vector).*nanstd(x_dist_init(find(x_dist_init(:,1)==2),2:no_columns))+h2_vector.*nanstd(x_dist(find(x_dist(:,1)==2),2:no_columns))); % offspring trait distribution std
            %x_dist(size(x_dist,1)+1,2:no_columns) = V4_pick_individuals(x_parent,off_std,1); % return trait
            x_dist(size(x_dist,1)+1,2:no_columns) = V4_pick_individuals(x_parent,off_std,1); % return trait -- this one increases variance
            x_dist(end,1) = 2;
            % tradeoff! -- overwrite d_min2 value (param 7) according to b_max1 (param 3)
            x_dist(end,7+1) = 2 * x_dist(end,3+1)^2;
            
        elseif Event_index == 4 % choose death of state 2 (R2)
            x_dist(whosnext(2),:) = []; % reduce dist by lost individual

        end
        
        R1 = sum(x_dist(:,1)==1); % count rows with State ID == 1
        R2 = sum(x_dist(:,1)==2); % count rows with State ID == 2
        
        if t > time_step
            %time_step
            %[R1 R2]
            pop_slice(1:no_species,time_step_index) = [R1 R2]; % assign current values to sliced standard times
            x_slice(1:no_params,time_step_index) = nanmean(x_dist(:,2:no_columns),1);
            x_var_slice(1:no_params,time_step_index) = nanvar(x_dist(:,2:no_columns),1);
            
            time_step_index = time_step_index + 1; % advance to next standardized time
            time_step = stand_times(time_step_index);
        end
        
        time_advance = exp(-1/CS_vector(end))/(CS_vector(end));
        if isnan(time_advance) == 0
            t = t + time_advance;
        else
            break
        end
    end
    
    % pass in current values to the final standardized time
        pop_slice(1:no_species,time_step_index) = [R1 R2];
        x_slice(1:no_params,end) = nanmean(x_dist(:,2:no_columns),1);
        x_var_slice(1:no_params,end) = nanvar(x_dist(:,2:no_columns),1);

    % pass sliced variable to standard time matrix
        pop_stand(:,:,i) = pop_slice;
        x_stand(:,:,i) = x_slice;
        x_var_stand(:,:,i) = x_var_slice;

end

%% calculate ci's for time series

pop_stand2 = pop_stand; pop_stand(pop_stand ==0) = NaN;
x_stand2 = x_stand; x_stand(x_stand ==0) = NaN;
x_var_stand2 = x_var_stand; x_var_stand(x_var_stand ==0) = NaN;

    upper_ci_level = 75; % choose ci levels
    lower_ci_level = 25; % choose ci levels    
    pop_data_out = V4_medians_and_cis(upper_ci_level,lower_ci_level,pop_stand2); % abundance      
    x_data_out = V4_medians_and_cis(upper_ci_level,lower_ci_level,x_stand2); % trait
    x_var_data_out = V4_medians_and_cis(upper_ci_level,lower_ci_level,x_var_stand2); % variance in trait 
    