clear; clc; clf;
figure(3); clf(3);
figure(4); clf(4);
figure(1); clf(1);
figure(12); clf(12);
figure(5); clf(5);
tic

% NEEDS THESE FUNCTIONS: 
% LTs_eqs_and_stability
% bdCompModel
% GEMv4_bdCompModel_bdTradeoff
% V4_pick_individuals
% V4_medians_and_cis
% jbfill

colors(1:6,1:3) = [[1 0.47 0]; [0.4 0 1]; [0.07 1 0]; [0.9 0 1]; [0.96 1 0]; [1 0 0.2]];
fill_colors = colors.*0.8;

num_states_to_plot = 2; %how many competiting species%
num_GEM_variations = 3; %how many variants of the simulation (evol, no evol, etc.)

% parameters for 2-species birth-death competition model
% values of b_max1 and b_max2 that are always equal initially
% but can lie below, above, or at the ESS,
% which is 1/2s (for s = 2 (default), ESS = 1/2*2 = 0.25)
    b_max1 = 0.25;           params(1)  = b_max1; % resource 1 max birth
    b11 = 0.0001;            params(2)  = b11; % resource 1 density-dependence of birth
    b_max2 = 0.25;           params(3)  = b_max2; % resource 1 max birth
    b22 = 0.0001;            params(4)  = b22; % resource 2 density-dependence of birth  
    d_min1 = 2*b_max1^2;          params(5)  = d_min1; % resource 1 min death
    d11 = 0.0001;            params(6)  = d11; % resource 1 density-dependence of death
    d_min2 = 2*b_max1^2;          params(7)  = d_min2; % resource 2 min death
    d22 = 0.0001;            params(8)  = d22; % resource 2 density-dependence of death  

% parameters for interspecific interactions
    b12 = 0.00012;                 params(9) = b12; % effect on resource 1 of resource 2 (births)
    b21 = 0.00012;                 params(10) = b21; % effect on resource 2 of resource 1 (births)
    d12 = 0.00012;                 params(11) = d12; % effect on resource 1 of resource 2 (deaths)
    d21 = 0.00012;                 params(12) = d21; % effect on resource 2 of resource 1 (deaths)
    
    
% Calculate b-d equivalent of K for each state, helpful for plot limits later
    K1 = (b_max1 - d_min1) / (b11 + d11);
    K2 = (b_max2 - d_min2) / (b22 + d22);
    Ks = [K1 K2]; % create a vector of them to use later
    
%% Find isoclines and equilibria for standard ODE
    % this calls on the LTs_eqs_and_stability function
    figure(12); clf(12);
    hold on; box on;
    [R1_hat_temp, R2_hat_temp, local_stability_temp, R1_iso, R2_vals, R2_iso, R1_vals] = ...
        LTs_eqs_and_stability(b_max1,b11, b_max2, b22, d_min1, d11, d_min2, d22, b12, b21, d12, d21);
    plot(R1_iso, R2_vals, '-k','LineWidth',2);
    plot(R1_vals, R2_iso, '--k','LineWidth',2);
 % add a vector field
    [R1_test,R2_test] = meshgrid(1:Ks(1)/30:Ks(1), 1:Ks(2)/30:Ks(2));
     U = (b_max1 - b11.*R1_test - b12.*R2_test).*R1_test - (d_min1 + d11.*R1_test + d12.*R2_test).*R1_test;
     V = (b_max2 - b22.*R2_test - b21.*R1_test).*R2_test - (d_min2 + d22.*R2_test + d21.*R1_test).*R2_test;
    quiver(R1_test,R2_test,U,V,5,'r')

    shg;
    
%% GEM setup and starting values
    % CVs and h2 for trait variation.
    % Rows are states (competitors), columns are GEM variants (no var, no
    % evol, evol).
    % In variant 1 (no var), CV = 0 and h2 = 0 for both competitors
    % In variant 2 (var, but no evol), CV = 0.2 and h2 = 0 for both
    % competitors
    % In variant 3 (evol), CV = 0.2 and h2 = 0.5 for both competitors
    cv = [0 0.2 0.2; 
          0 0.2 0.2];
    h_2 = [0  0  0.5; 
           0  0  0.5];

    % Specify initial ecological state for both competitors
    % NOTE: change these manually to generate simulation output for
    % different starting conditions. Plots were exported and combined in
    % MS Powerpoint to generate figures 3 and 4 in the main text.

    % y0 = [25, 25]; % below the saddle point equilibrium
    % y0 = [284, 284]; % saddle point for b_max = 0.25 and d_min = 2*b_max^2 
     y0 = [545, 545]; % above the saddle point equilibrium

    num_replicates = 250; % number of GEM simulations
    t_max = 300; % time span to run simulations
    tspan = [0 t_max]; % start and end times
    
%% pre-visualize standard ODE solution model dynamics
    % requires function bdCompModel 
    ode = @(t,y) bdCompModel(t,y,b_max1,b11, b_max2, b22, d_min1, d11, d_min2, d22, b12, b21, d12, d21); % compile function and call
    [t1,y1] = ode45(ode, tspan, y0); % return time and population density vectors
    figure(1);clf(1);
    hold on;
    plot(t1,y1,'-');
    ylim([0 (max(Ks) + 0.15*max(Ks))])
    legend({'R1','R2'});
    
    shg;

%% specify which parameters go with which states
    state_parameter_match = zeros(length(y0),length(params)); % start with all 0s for no
    state_parameter_match(1,[1,2,5,6,9,11]) = 1; % state R1 (1) parameters
    state_parameter_match(2,[3,4,7,8,10,12]) = 1; % state R2 (2) parameters
    
%% set up function for finding NEEAs from abundances
    syms bmax_R1 bmax_R2 R1_S R2_S
        
        LRS_R1 = (bmax_R1 - b11.*R1_S - b12.*R2_S) / (2*bmax_R1.^2 + d11.*R1_S + d12.*R2_S);
        LRS_R2 = (bmax_R2 - b22.*R2_S - b21.*R1_S) / (2*bmax_R2.^2 + d22.*R2_S + d21.*R1_S);
        % take derivative
        d_LRS_R1 = diff(LRS_R1, bmax_R1);
        d_LRS_R2 = diff(LRS_R2, bmax_R2);
        % find where derivative = 0
        bmax_R1_hat = solve(d_LRS_R1 == 0); 
        bmax_R2_hat = solve(d_LRS_R2 == 0);
%% run GEM algorithm for each GEM variation, plot results
% For each GEM variant (no var, no evol, evol)...
for j = 1:num_GEM_variations
    % specify which parameters should evolve and have variability
    % create empty vector of CVs associated with each trait
        cv_vector = zeros(1,length(params)); % start with all 0s for no variance
    % Set the evolving parameter(s)' variation according to actual cv, 
    % which is stored in cv(state, GEM variation)
        cv_vector(1) = cv(1,j); % R1's b_max
        cv_vector(3) = cv(2,j); % R2's b_max
        cv_vector(5) = cv(1,j); % R1's d_min
        cv_vector(7) = cv(2,j); % R2's d_min
    
    % Set the evolving parameter(s)' heritability according to actual h^2, 
    % which is stored in h_2(state, GEM variation)
        h2_vector = zeros(1,length(params)); % start with all 0s for not heritable
    % Now set the evolving parameter(s)' h^2
        h2_vector(1) = h_2(1,j); % R1's b_max
        h2_vector(3) = h_2(2,j); % R2's b_max
        h2_vector(5) = h_2(1,j); % R1's d_min
        h2_vector(7) = h_2(2,j); % R2's d_min
    
    % run GEM %
    [x_stand, pop_stand, stand_times, pop_data_out, x_data_out, x_var_data_out] = GEMv4_bdCompModel_bdTradeoff(params, state_parameter_match, cv_vector, h2_vector, num_replicates, y0, t_max);

    % calculate realized initial values of evolving traits
    % param medians and CIs are in x_data_out(statistic, t, param)
    % where statistic 1 = lower CI, 2 = median, 3 = upper CI
    realized_initial_b_max1 = x_data_out(2,1,1); % R1's b_max
    realized_initial_b_max2 = x_data_out(2,1,3); % R2's b_max
    realized_initial_d_min1 = x_data_out(2,1,5); % R1's d_min
    realized_initial_d_min2 = x_data_out(2,1,7); % R2's d_min
    
    ode = @(t,y) bdCompModel(t,y,realized_initial_b_max1,b11, realized_initial_b_max2,...
        b22, realized_initial_d_min1, d11, realized_initial_d_min2, d22,b12, b21, d12, d21); % compile function and call
    [t1,y1] = ode45(ode, tspan, y0); % return time and population density vectors
      
    % Calculate the cumulative proportion of replicates
    % in which an extinction occured as a function of time,
    % and store in cumExtTraj(time, GEM_variation).
    % While in this loop, calculate the proportion of replicates
    % with an internal equil at each t for each gem variation,
    % store this in PropEquil_t(time, GEM_variation).
    % Note: calculation of PropEquil_t is such that extinct
    % populations get counted- they get "frozen" on their last value
    % (0 or 1). So PropEquil_t is more precisely a measure of the
    % proportion of replicates that either still have an unstable
    % equilibrium or went extinct with one. 
    % Also calculate and plot the local stability of all replicates
    % with an internal equil at each t for each gem variation,
    % store this in Median_LocStab(time, GEM_variation
    % Note: This necessarily excludes communities that went extinct,
    % but still had an unstable equil. when they did.
    if j == 1 % if on first GEM iteration, initialize those arrays
        cumExtTraj = nan(length(stand_times), num_GEM_variations);
        PropEquil_t = nan(length(stand_times), num_GEM_variations);
        PropStEquil_t = nan(length(stand_times), num_GEM_variations); % one for stable equil, too, JIC
        Median_LocStab = nan(length(stand_times), num_GEM_variations);
        LoCI_LocStab = nan(length(stand_times), num_GEM_variations);
        UpCI_LocStab = nan(length(stand_times), num_GEM_variations);
    end
    % initialize similar arrays to hold summaries for each replicate
    Ext = nan(num_replicates, length(stand_times));
    UnstEquil = nan(num_replicates, length(stand_times));
    StEquil = nan(num_replicates, length(stand_times));
    LocStab = nan(num_replicates, length(stand_times));
    Abund_i = nan(num_replicates, length(stand_times));
    Abund_j = nan(num_replicates, length(stand_times));
    Trait_i = nan(num_replicates, length(stand_times));
    Trait_j = nan(num_replicates, length(stand_times));
    % if j = 3 (evolutionary scenario), then need to calculate
    % summary statistics for each replicate along with NEEAs
    if j == 3
        % add arrays to hold the NEEAs at each time for each competitor pop.
        % for each replicate
        NEEA_i = nan(num_replicates, length(stand_times));
        NEEA_j = nan(num_replicates, length(stand_times));
        
        % parallel loop over replicates within this GEM variant j
        parfor ii = 1:num_replicates
            % first store the abundances for the replicate, to use in 
            % calculating NEEAs below
            % need pop_stand(state, t, replicate) for this
            this_rep = pop_stand(:, :, ii);
            % need missing values as NAs for NEEA calc.
            this_rep(isnan(this_rep)) = 0;
            %
            Ext_thisrep = Ext(ii,:);
            UnstEquil_thisrep = UnstEquil(ii,:)
            StEquil_thisrep = StEquil(ii,:)
            LocStab_thisrep = LocStab(ii,:)
            NEEA_thisi = NEEA_i(ii,:);
            NEEA_thisj = NEEA_j(ii,:);
            for t = 1:length(stand_times)
                % find NEEA for each competitor population
                bmax_R1_hat_t = double(subs(bmax_R1_hat,...
                                            [R1_S, R2_S],...
                                            [this_rep(1,t), this_rep(2,t)]));
                bmax_R2_hat_t = double(subs(bmax_R2_hat,...
                                            [R1_S, R2_S],...
                                            [this_rep(1,t), this_rep(2,t)]));
                % end up with a negative and a positive root, store positive one
                NEEA_thisi(t) = bmax_R1_hat_t(bmax_R1_hat_t > 0);
                NEEA_thisj(t) = bmax_R2_hat_t(bmax_R2_hat_t > 0);
                % if both states have non-zero abundances...
                if ~isnan(x_stand(1,t,ii)) && ~isnan(x_stand(3,t,ii))
                    Ext_thisrep(t)= 0; % set rep. extinction = 0
                    %find traits at t for replicate ii and check for
                    % equilibrium (do this using inequality condition,
                    % K1/a12 < K2 && K2/a21 < K1).
                    % Note that K = (bmax - dmin) / (b11 - d11),
                    % and a12 = (b12 + d12) / (b11 + d11),
                    % so check that (b_max1 - d_min1) / (b12 + d12) < 
                    % (b_max2 - dmin2) / (b22 + d22), and vice versa
                    K1_div_a12 = (x_stand(1,t,ii) - x_stand(5,t,ii)) / (x_stand(9,t,ii) + x_stand(11,t,ii));
                    K2 = (x_stand(3,t,ii) - x_stand(7,t,ii)) / (x_stand(4,t,ii) + x_stand(8,t,ii));
                    K2_div_a21 = (x_stand(3,t,ii) - x_stand(7,t,ii)) / (x_stand(10,t,ii) + x_stand(12,t,ii));
                    K1 = (x_stand(1,t,ii) - x_stand(5,t,ii)) / (x_stand(2,t,ii) + x_stand(6,t,ii));
                    if K1_div_a12 < K2 && K2_div_a21 < K1 % if there's an unstable coexistence equilibrium
                        UnstEquil_thisrep(t) = 1;
                        StEquil_thisrep(t) = 0;
                        % calculate local stability for the replicate
                        [R1_hat, R2_hat, local_stability, R1_iso, R2_vals, R2_iso, R1_vals] = ...
                           LTs_eqs_and_stability(x_stand(1,t,ii), x_stand(2,t,ii), x_stand(3,t,ii), x_stand(4,t,ii), ...
                                                 x_stand(5,t,ii), x_stand(6,t,ii), x_stand(7,t,ii), x_stand(8,t,ii),...
                                                 x_stand(9,t,ii), x_stand(10,t,ii), x_stand(11,t,ii), x_stand(12,t,ii));
                        % store the value for the replicate
                        LocStab_thisrep(t) = local_stability;
                    elseif K1_div_a12 > K2 && K2_div_a21 > K1 % if there's a stable coexistence equilibrium
                        UnstEquil_thisrep(t) = 0;
                        StEquil_thisrep(t) = 1;
                    else % if neither of the above are true, there is no equilibrium or either kind
                        UnstEquil_thisrep(t) = 0;
                        StEquil_thisrep(t) = 0;
                    end
                else % otherwise, at least one of the states has zero abundance...
                    Ext_thisrep(t)= 1; % set rep. extinction = 1
                    UnstEquil_thisrep(t) = UnstEquil_thisrep(t-1); % set unst_equil_rep equal to value in previous time step
                    StEquil_thisrep(t) = StEquil_thisrep(t-1); % set st_equil_rep equal to value in previous time step
                    % no need to update eigenvalue-- leave NA
                end 
            end
            Ext(ii, :) = Ext_thisrep;
            UnstEquil(ii,:) =  UnstEquil_thisrep;
            StEquil(ii,:) =  StEquil_thisrep;
            LocStab(ii,:) = LocStab_thisrep;
            NEEA_i(ii,:) = NEEA_thisi;
            NEEA_j(ii,:) = NEEA_thisj;
            Abund_i(ii,:) = pop_stand(1, :, ii);
            Abund_j(ii,:) = pop_stand(2, :, ii);
            Trait_i(ii,:) = x_stand(1,:,ii);
            Trait_j(ii,:) = x_stand(3,:,ii);
        end % concludes loop over replicates for evolutionary GEM variant
    else % if j = 1 or 2 (no variation and no evolution GEM variants)
        parfor ii = 1:num_replicates
            % first store the abundances for the replicate, to use in 
            % calculating NEEAs below
            % need pop_stand(state, t, replicate) for this
            this_rep = pop_stand(:, :, ii);
            % need missing values as NAs for NEEA calc.
            this_rep(isnan(this_rep)) = 0;
            %
            Ext_thisrep = Ext(ii,:);
            UnstEquil_thisrep = UnstEquil(ii,:)
            StEquil_thisrep = StEquil(ii,:)
            LocStab_thisrep = LocStab(ii,:)
            for t = 1:length(stand_times)
                % if both states have non-zero abundances...
                if ~isnan(x_stand(1,t,ii)) && ~isnan(x_stand(3,t,ii))
                    Ext_thisrep(t)= 0; % set rep. extinction = 0
                    %find traits at t for replicate ii and check for
                    % equilibrium (do this using inequality condition,
                    % K1/a12 < K2 && K2/a21 < K1).
                    % Note that K = (bmax - dmin) / (b11 - d11),
                    % and a12 = (b12 + d12) / (b11 + d11),
                    % so check that (b_max1 - d_min1) / (b12 + d12) < 
                    % (b_max2 - dmin2) / (b22 + d22), and vice versa
                    K1_div_a12 = (x_stand(1,t,ii) - x_stand(5,t,ii)) / (x_stand(9,t,ii) + x_stand(11,t,ii));
                    K2 = (x_stand(3,t,ii) - x_stand(7,t,ii)) / (x_stand(4,t,ii) + x_stand(8,t,ii));
                    K2_div_a21 = (x_stand(3,t,ii) - x_stand(7,t,ii)) / (x_stand(10,t,ii) + x_stand(12,t,ii));
                    K1 = (x_stand(1,t,ii) - x_stand(5,t,ii)) / (x_stand(2,t,ii) + x_stand(6,t,ii));
                    if K1_div_a12 < K2 && K2_div_a21 < K1 % if there's an unstable coexistence equilibrium
                        UnstEquil_thisrep(t) = 1;
                        StEquil_thisrep(t) = 0;
                        % calculate local stability for the replicate
                        [R1_hat, R2_hat, local_stability, R1_iso, R2_vals, R2_iso, R1_vals] = ...
                           LTs_eqs_and_stability(x_stand(1,t,ii), x_stand(2,t,ii), x_stand(3,t,ii), x_stand(4,t,ii), ...
                                                 x_stand(5,t,ii), x_stand(6,t,ii), x_stand(7,t,ii), x_stand(8,t,ii),...
                                                 x_stand(9,t,ii), x_stand(10,t,ii), x_stand(11,t,ii), x_stand(12,t,ii));
                        % store the value for the replicate
                        LocStab_thisrep(t) = local_stability;
                    elseif K1_div_a12 > K2 && K2_div_a21 > K1 % if there's a stable coexistence equilibrium
                        UnstEquil_thisrep(t) = 0;
                        StEquil_thisrep(t) = 1;
                    else % if neither of the above are true, there is no equilibrium or either kind
                        UnstEquil_thisrep(t) = 0;
                        StEquil_thisrep(t) = 0;
                    end
                else % otherwise, at least one of the states has zero abundance...
                    Ext_thisrep(t)= 1; % set rep. extinction = 1
                    UnstEquil_thisrep(t) = UnstEquil_thisrep(t-1); % set unst_equil_rep equal to value in previous time step
                    StEquil_thisrep(t) = StEquil_thisrep(t-1); % set st_equil_rep equal to value in previous time step
                    % no need to update eigenvalue-- leave NA
                end 
            end
            Ext(ii, :) = Ext_thisrep;
            UnstEquil(ii,:) =  UnstEquil_thisrep;
            StEquil(ii,:) =  StEquil_thisrep;
            LocStab(ii,:) = LocStab_thisrep;
            Abund_i(ii,:) = pop_stand(1, :, ii);
            Abund_j(ii,:) = pop_stand(2, :, ii);
            Trait_i(ii,:) = x_stand(1,:,ii);
            Trait_j(ii,:) = x_stand(3,:,ii);
        end
    end

    % export Ext as a .csv for survival analysis
    my_field = strcat('j',num2str(j));
    myfilename = append('Extinctions_', my_field, ...
        ".csv");
    writematrix(Ext, myfilename);
    % also export traits and abundances for all replicates
    % trait for pop i
    myfilename2 = append('traits_i_', my_field, ...
        ".csv");
    writematrix(Trait_i, myfilename2);
    % trait for pop j
    myfilename3 = append('traits_j_', my_field, ...
        ".csv");
    writematrix(Trait_j, myfilename3);
    % abund for pop i
    myfilename4 = append('abundances_i_', my_field, ...
        ".csv");
    writematrix(Abund_i, myfilename4);
    % abund for pop j
    myfilename5 = append('abundances_j_', my_field, ...
        ".csv");
    writematrix(Abund_j, myfilename5);

    % Figure 3 %
    % calculate summaries of those values across replicates for each t
    for t = 1:length(stand_times)
        cumExtTraj(t, j) = sum(Ext(:, t)) / num_replicates;
        PropEquil_t(t, j) = sum(UnstEquil(:, t)) / num_replicates;
        PropStEquil_t(t, j) = sum(StEquil(:, t)) / num_replicates; % one for stable equil, too, JIC
        Median_LocStab(t, j) = prctile(LocStab(:, t), 50);
        LoCI_LocStab(t, j) = prctile(LocStab(:, t), 25);
        UpCI_LocStab(t, j) = prctile(LocStab(:, t), 75);
    end
    % make a plot of (A) cumulative proportion of replicates with extinctions,
    % (B) proportion of replicates with an unstable equilibrium, and
    % (C) the real part of the dominant eigenvalue at the unstable equilibrium,
    % in cases where it exists.
    figure(3);
        subplot(3, 1, 1);
            box on; hold on;
            plot(stand_times, cumExtTraj(:, j), '-', 'Color', colors(j,:), 'Linewidth', 3);
            ylim([0 1]);
            if j == num_GEM_variations
                legend({'No variation','No evolution','Evolution'}, 'Location', 'best');
                xlabel('time')
                ylabel('cumulative extinctions')
            end
        subplot(3, 1, 2);
            box on; hold on;
            plot(stand_times, PropEquil_t(:, j), '-', 'Color', colors(j,:), 'Linewidth', 3);
            ylim([0 1]);
            if j == num_GEM_variations
                legend({'No variation','No evolution','Evolution'},'Location', 'best');
                xlabel('time')
                ylabel('prop. replicates with equilibrium')
            end
        subplot(3, 1, 3);
            box on; hold on;
            %jbfill(stand_times, UpCI_LocStab(:, j)', LoCI_LocStab(:, j)', colors(j,:),'w',1,0.2); hold on;
            %patch([stand_times, fliplr(stand_times)], [UpCI_LocStab(:, j)', fliplr(LoCI_LocStab(:, j)')], colors(j,:), 'FaceAlpha', 0.3);
            UpCI_LocStab_trim = UpCI_LocStab(~isnan(UpCI_LocStab(:,j)), j);
            LoCI_LocStab_trim = LoCI_LocStab(~isnan(LoCI_LocStab(:,j)), j);
            stand_times_trim = stand_times(~isnan(UpCI_LocStab(:,j)));
            patch([stand_times_trim, fliplr(stand_times_trim)],...
            [UpCI_LocStab_trim', fliplr(LoCI_LocStab_trim')],...
            colors(j,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            % plot median, store it as an object to later use to specify
            % the legend values
            locstabmed(j) = plot(stand_times, Median_LocStab(:,j),'Color',colors(j,:),'LineWidth',3);
            if j == num_GEM_variations
                legend(locstabmed, {'No variation','No evolution','Evolution'});
                xlabel('time')
                ylabel('local stability')
            end

    figure(4); % plot traits
    % create an identifier to tack onto the name for each plot,
    % so we can call them back up after the GEM variants all run to add
    % the NEEAs and ESS values
    % Can also go ahead and add panels showing the median trait values for
    % cases in which the focal pop (R1) is alone at tmax, 
    % extinct at tmax, or where both pops are still present at tmax
    my_field = strcat('j',num2str(j));
    box on; hold on;
        subplot(4,1,1); % plot for all replicates
            box on; hold on;
            xlim([0 t_max]);
            ylim([0 0.4]);
            % plot CIs using jbfill function
            % recall x_data_out(statistic, time, parameter)
            % want to plot lower (stat = 1) and upper (stat = 3) CIs
            % for focal population i bmax (parameter = 1)
            jbfill(stand_times, x_data_out(1,:,1), x_data_out(3,:,1),colors(j,:),'w',1,0.2);
            hold on;
            % plot median (stat = 2)
            plot4.a.(my_field) = plot(stand_times,x_data_out(2,:,1),'-k','Color',colors(j,:),'LineWidth',2);
            ylabel(['maximum birth rate, b_{' 'max_' '1' '}'],'FontSize',12);
    if j == 3
        % Add the NEEA(t) for bmaxs to plots of trait trajectories.
        % want to subdivide this by abundances at tmax 
        NEEAs_i_noR2 = NEEA_i(isnan(pop_stand(2,end,:)),:);
        NEEAs_j_noR2 = NEEA_j(isnan(pop_stand(2,end,:)),:);
        NEEAs_i_noR1 = NEEA_i(isnan(pop_stand(1,end,:)),:);
        NEEAs_j_noR1 = NEEA_j(isnan(pop_stand(1,end,:)),:);
        NEEAs_i_R1ANDR2 = NEEA_i(~isnan(pop_stand(1,end,:).* pop_stand(2,end,:)),:);
        NEEAs_j_R1ANDR2 = NEEA_j(~isnan(pop_stand(1,end,:).* pop_stand(2,end,:)),:);

        % need to summarize these into medians and CIs
        % Note, also need to do this for the NEEAs for all replicates
        %summary_NEEAs_i = V4_medians_and_cis(75, 25, NEEA_i);
        summary_NEEAs_i = [prctile(NEEA_i, 25, 1);
                           prctile(NEEA_i, 50, 1);
                           prctile(NEEA_i, 75, 1)];
        summary_NEEAs_j = [prctile(NEEA_j, 25, 1);
                           prctile(NEEA_j, 50, 1);
                           prctile(NEEA_j, 75, 1)];
        summary_NEEAs_i_noR2 = [prctile(NEEAs_i_noR2, 25, 1);
                                prctile(NEEAs_i_noR2, 50, 1);
                                prctile(NEEAs_i_noR2, 75, 1)];
        summary_NEEAs_j_noR2 = [prctile(NEEAs_j_noR2, 25, 1);
                                prctile(NEEAs_j_noR2, 50, 1);
                                prctile(NEEAs_j_noR2, 75, 1)];
        summary_NEEAs_i_noR1 = [prctile(NEEAs_i_noR1, 25, 1);
                                prctile(NEEAs_i_noR1, 50, 1);
                                prctile(NEEAs_i_noR1, 75, 1)];
        summary_NEEAs_j_noR1 = [prctile(NEEAs_j_noR1, 25, 1);
                                prctile(NEEAs_j_noR1, 50, 1);
                                prctile(NEEAs_j_noR1, 75, 1)];
        summary_NEEAs_i_R1ANDR2 = [prctile(NEEAs_i_R1ANDR2, 25, 1);
                                   prctile(NEEAs_i_R1ANDR2, 50, 1);
                                   prctile(NEEAs_i_R1ANDR2, 75, 1)];
        summary_NEEAs_j_R1ANDR2 = [prctile(NEEAs_j_R1ANDR2, 25, 1);
                                   prctile(NEEAs_j_R1ANDR2, 50, 1);
                                   prctile(NEEAs_j_R1ANDR2, 75, 1)];

        % also need the median +/-CIs of trait trajectories for each subset
        % need x_stand(trait, t, replicate) for this (bmax1 and 2 are
        % traits 1 and 3)
        % store results as med_bmax_XXXX = (statistic, time),
        % where statistic = 1 = 25th percentile, 2 = median, and 3 = 75th 
        % first for no R2
        med_bmax1_noR2 = V4_medians_and_cis(75, 25, x_stand(1,:,isnan(pop_stand(2,end,:))));
        med_bmax2_noR2 = V4_medians_and_cis(75, 25, x_stand(3,:,isnan(pop_stand(2,end,:))));
        % now for no R1
        med_bmax1_noR1 = V4_medians_and_cis(75, 25, x_stand(1,:,isnan(pop_stand(1,end,:))));
        med_bmax2_noR1 = V4_medians_and_cis(75, 25, x_stand(3,:,isnan(pop_stand(1,end,:))));
        % now for both present
        med_bmax1_R1ANDR2 = V4_medians_and_cis(75, 25, x_stand(1,:,...
                                            ~isnan(pop_stand(1,end,:) .* pop_stand(2,end,:))));
        med_bmax2_R1ANDR2 = V4_medians_and_cis(75, 25, x_stand(3,:,...
                                            ~isnan(pop_stand(1,end,:) .* pop_stand(2,end,:))));

        %recall figure(4) to add trait medians, NEEAs, and ESSs
        figure(4);
        hold on;
            subplot(4,1,1);
                xlim([0 t_max]);
                ylim([0 0.4]);
                xlabel('time', 'FontSize', 12);
                % plot NEEA CIs using jbfill function
                jbfill(stand_times, summary_NEEAs_i(1,:), summary_NEEAs_i(3,:),colors(4,:),'w',1,0.2);
                hold on;
                % save as objects to refer to in the legend
                plot4.a.NEEA = plot(stand_times, summary_NEEAs_i(2,:),'-k','Color',colors(4,:),'LineWidth',2);
                plot4.a.ESS = plot(stand_times, repelem(0.25,length(stand_times)),'-k','Color','black','LineWidth',2);
                legend([plot4.a.j1, plot4.a.j2, plot4.a.j3, plot4.a.NEEA, plot4.a.ESS],...
                    {'No variation', 'No evolution', 'Evolution','NEEA','ESS'},'Location', 'best');
                title('All Replicates');
            subplot(4,1,2);
                box on; hold on;
                xlim([0 t_max]);
                ylim([0 0.4]);
                xlabel('time', 'FontSize', 12);
                % plot CIs for trait using jbfill function
                jbfill(stand_times, med_bmax1_noR2(1,:), med_bmax1_noR2(3,:),colors(j,:),'w',1,0.2);
                hold on;
                % plot median trait
                plot4.b.j3 = plot(stand_times, med_bmax1_noR2(2,:),'-k','Color',colors(j,:),'LineWidth',2);
                xlabel('time', 'FontSize', 12);
                hold on;
                % plot NEEa CIs
                jbfill(stand_times, summary_NEEAs_i_noR2(1,:), summary_NEEAs_i_noR2(3,:),colors(4,:),'w',1,0.2);
                hold on;
                % add median NEEA
                plot4.b.NEEA = plot(stand_times, summary_NEEAs_i_noR2(2,:),'-k','Color',colors(4,:),'LineWidth',2);
                % add ESS
                plot4.b.ESS = plot(stand_times, repelem(0.25,length(stand_times)),'-k','Color','black','LineWidth',2);
                % make a legend
%                 legend([plot4.b.j3, plot4.b.NEEA, plot4.b.ESS],...
%                     {'Evolution','NEEA','ESS'},'Location', 'best');
                title('R2 Extinct');
            subplot(4,1,3);
                box on; hold on;
                xlim([0 t_max]);
                ylim([0 0.4]);
                xlabel('time', 'FontSize', 12);
                % need to use something other than jbfill here b/c median
                % eventually becomes NaN (all these replicates have r1 extinct at
                % max_t)
                UpCI_med_bmax1_noR1 = med_bmax1_noR1(3, ~isnan(med_bmax1_noR1(3,:)));
                LoCI_med_bmax1_noR1 = med_bmax1_noR1(1, ~isnan(med_bmax1_noR1(3,:)));
                stand_times_trim = stand_times(~isnan(med_bmax1_noR1(3,:)));
                patch([stand_times_trim, fliplr(stand_times_trim)],...
                    [LoCI_med_bmax1_noR1, fliplr(UpCI_med_bmax1_noR1)],...
                    colors(j,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                % add median trait
                plot4.c.j3 = plot(stand_times, med_bmax1_noR1(2,:),'-k','Color',colors(j,:),'LineWidth',2);
                hold on;
                % plot NEEa CIs
                jbfill(stand_times, summary_NEEAs_i_noR1(1,:), summary_NEEAs_i_noR1(3,:),colors(4,:),'w',1,0.2);
                hold on;
                % add median NEEA
                plot4.c.NEEA = plot(stand_times, summary_NEEAs_i_noR1(2,:),'-k','Color',colors(4,:),'LineWidth',2);
                plot4.c.ESS = plot(stand_times, repelem(0.25,length(stand_times)),'-k','Color','black','LineWidth',2);
%                 legend([plot4.c.j3, plot4.c.NEEA, plot4.c.ESS],...
%                     {'Evolution','NEEA','ESS'},'Location', 'best');
                title('R1 Extinct');
            subplot(4,1,4);
                box on; hold on;
                xlim([0 t_max])
                ylim([0 0.4]);
                xlabel('time', 'FontSize', 12);
                % plot CIs using jbfill function
                jbfill(stand_times, med_bmax1_R1ANDR2(1,:), med_bmax1_R1ANDR2(3,:),colors(j,:),'w',1,0.2);
                hold on;
                % add median trait
                plot4.d.j3 = plot(stand_times, med_bmax1_R1ANDR2(2,:),'-k','Color',colors(j,:),'LineWidth',2);
                % plot NEEa CIs
                jbfill(stand_times, summary_NEEAs_i_R1ANDR2(1,:), summary_NEEAs_i_R1ANDR2(3,:),colors(4,:),'w',1,0.2);
                hold on;
                % add median NEEA
                plot4.d.NEEA = plot(stand_times, summary_NEEAs_i_R1ANDR2(2,:),'-k','Color',colors(4,:),'LineWidth',2);
                plot4.d.ESS = plot(stand_times, repelem(0.25,length(stand_times)),'-k','Color','black','LineWidth',2);
%                 legend([plot4.d.j3, plot4.d.NEEA, plot4.d.ESS],...
%                     {'Evolution','NEEA','ESS'},'Location', 'best');
                title('Neither Extinct');
               hold off;       
        
        % For the evolutionary scenario, can also
        % plot competitor trait values against each other, colored by 
        % coexistence duration
        figure(5);
        box on; hold on;
        xlim([0 0.4]);
        ylim([0 0.4]);
        scatter(Trait_i(:,length(Trait_i)), Trait_j(:,length(Trait_j)));
    end 
end
toc