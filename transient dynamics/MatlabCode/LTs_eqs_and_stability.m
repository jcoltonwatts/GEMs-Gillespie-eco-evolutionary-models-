function [R1_hat1, R2_hat1, eigen_v_max, R1_at_R2, R2_range, R2_at_R1, R1_range] = ...
    LTs_eqs_and_stability(b_max1, b11, b_max2, b22, d_min1, d11, d_min2, d22, b12, b21, d12, d21)

%% set up Jacobian

% Declare variables
syms R1s R2s b_max1s b11s b_max2s b22s d_min1s d11s d_min2s d22s b12s b21s d12s d21s

% Find Jacobian from equations
    % Needs to be of the form:
    % jac = jacobian([resource1_equation,resource2_equation],[R1,R2])
jacobian_matrix = jacobian([R1s*(b_max1s - b11s*R1s - b12s*R2s) - R1s*(d_min1s + d11s*R1s + d12s*R2s),... 
    R2s*(b_max2s - b22s*R2s - b21s*R1s) - R2s*(d_min2s + d22s*R2s + d21s*R1s)], [R1s, R2s]);

%% find equilibria via isoclines
% For each resource, only need to consider values of the other resource
% from 0 to K/alpha
R2_range = 0:1:((b_max1 - d_min1)/(b12 + d12)); % pick a R2 range 
R1_at_R2 = ((b_max1 - d_min1)/(b11 + d11)) - ((b12 + d12)/(b11 + d11)).*R2_range; % calculate R1 from isocline
R1_range = 0:1:((b_max2 - d_min2)/(b21 + d21)); % pick a R1 range  
R2_at_R1 = (b_max2 - d_min2)/(b22 + d22) - ((b21 + d21)/(b22 + d22)).*R1_range; % calculate R2 from isocline

% Can uncomment these to check the first offending b_max value
% R2_range = 0:1:((b_maxs(1) - d_min1)/(b12 + d12)); % pick a R2 range 
% R1_at_R2 = ((b_maxs(1) - d_min1)/(b11 + d11)) - ((b12 + d12)/(b11 + d11)).*R2_range; % calculate R1 from isocline
% R1_range = 0:1:((b_maxs(1) - d_min2)/(b21 + d21)); % pick a R1 range  
% R2_at_R1 = (b_maxs(1) - d_min2)/(b22 + d22) - ((b21 + d21)/(b22 + d22)).*R1_range; % calculate R2 from isocline
% 
% figure(99);
% plot(R1_range, R2_at_R1, R1_at_R2, R2_range)

% find the intersections of the two isoclines
[R1_hat1_maybe, R2_hat1_maybe] = intersections(R1_range, R2_at_R1, R1_at_R2, R2_range);
if isempty([R1_hat1_maybe, R2_hat1_maybe])
    %disp('No internal equilibrium')
    eigen_v_max = NaN;
    eigen_v_min = NaN;
    R1_hat1 = NaN;
    R2_hat1 = NaN;
elseif length(R1_hat1_maybe) > 1
    % are there different equilibria (isoclines are the same/overlap)?
    % check if there are as many equilibria as points checked
    if length(R1_hat1_maybe) == length(R2_range)
        %disp('isoclines overlap');
        eigen_v_max = NaN;
        eigen_v_min = NaN;
    % if not, "intersections.m" has returned duplicates. Just keep the
    % first
    else
        R1_hat1 = R1_hat1_maybe(1);
        R2_hat1 = R2_hat1_maybe(1);
        
    %% use isoclines and Jacobian to calculate local stability

    % substitute rep i's parameters into jacobian matrix
    num_jac = subs(jacobian_matrix,...
        [R1s, R2s, b_max1s, b11s, b_max2s, b22s, d_min1s, d11s, d_min2s, d22s, b12s, b21s, d12s, d21s],...
        [R1_hat1,R2_hat1,b_max1, b11, b_max2, b22, d_min1, d11, d_min2, d22, b12, b21, d12, d21]);

   % num_jac = subs(jacobian_matrix,...
   %     [R1s, R2s, b_max1s, b11s, b_max2s, b22s, d_min1s, d11s, d_min2s, d22s, b12s, b21s, d12s, d21s],...
   %     [R1_hat1,R2_hat1,x_data_out(2,end,1), x_data_out(2,end,2), x_data_out(2,end,3), x_data_out(2,end,4), ...
   %                                  x_data_out(2,end,5), x_data_out(2,end,6), x_data_out(2,end,7), x_data_out(2,end,8),...
   %                                  x_data_out(2,end,9), x_data_out(2,end,10), x_data_out(2,end,12), x_data_out(2,end,12)]);

    % Find eigenvalues
    eig(num_jac);
        eigen_v_max = double(max(real(eig(num_jac))));
        eigen_v_min = double(min(real(eig(num_jac))));
    end
elseif length(R1_hat1_maybe) == 1
    [R1_hat1, R2_hat1] = intersections(R1_range, R2_at_R1, R1_at_R2, R2_range);

    %% use isoclines and Jacobian to calculate local stability

    % substitute rep i's parameters into jacobian matrix
    num_jac = subs(jacobian_matrix,...
        [R1s, R2s, b_max1s, b11s, b_max2s, b22s, d_min1s, d11s, d_min2s, d22s, b12s, b21s, d12s, d21s],...
        [R1_hat1,R2_hat1,b_max1, b11, b_max2, b22, d_min1, d11, d_min2, d22, b12, b21, d12, d21]);

   % num_jac = subs(jacobian_matrix,...
   %     [R1s, R2s, b_max1s, b11s, b_max2s, b22s, d_min1s, d11s, d_min2s, d22s, b12s, b21s, d12s, d21s],...
   %     [R1_hat1,R2_hat1,x_data_out(2,end,1), x_data_out(2,end,2), x_data_out(2,end,3), x_data_out(2,end,4), ...
   %                                  x_data_out(2,end,5), x_data_out(2,end,6), x_data_out(2,end,7), x_data_out(2,end,8),...
   %                                  x_data_out(2,end,9), x_data_out(2,end,10), x_data_out(2,end,12), x_data_out(2,end,12)]);

    % Find eigenvalues
    eig(num_jac);
        eigen_v_max = double(max(real(eig(num_jac))));
        eigen_v_min = double(min(real(eig(num_jac))));
end

