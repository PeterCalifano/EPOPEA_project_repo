close all
clear
clc


%% Earth --> Mars transfer optimization
Is = 3000; % Specific Impulse of the thruster

n_sol = 100; % Number of points for the trajectory 
n_integrator = 100;

dep_dates.initial = [2025, 1, 1, 12, 0, 0]; 
dep_dates.final = [2030, 1, 1, 12, 0, 0];

time_unit = 365.5;

dep_JD.initial = date2mjd2000(dep_dates.initial);
dep_JD.final = date2mjd2000(dep_dates.final);

arr_JD.initial = date2mjd2000(dep_dates.initial) + time_unit;
arr_JD.final = date2mjd2000(dep_dates.final) + 3*time_unit;

% Number of points of the grid
n_points = 100;
% Final mass
M_end = 1000;

grid_launch_win = linspace(dep_JD.initial, dep_JD.final, n_points);
grid_arrival_win = linspace(arr_JD.initial, arr_JD.final, n_points);

m_cost = nan(length(grid_launch_win), length(grid_arrival_win));
T_max = nan(length(grid_launch_win), length(grid_arrival_win));


planet1 = 3; % Earth
planet2 = 4; % Mars


N = 1; % Number of revs


for i = 1:length(grid_launch_win)
    tic
    JD_launch = grid_launch_win(i);
    parfor j = 1:length(grid_arrival_win)

        JD_arrival = grid_arrival_win(j);

        if (JD_arrival - JD_launch) > 0

            % Find RI, VI, RF, VF
            [X_1, V_1, ~] = SSephem_1800_2050(planet1, JD_launch);
            [X_2, V_2, ~] = SSephem_1800_2050(planet2, JD_arrival);

            TOF = 3600*24*(JD_arrival - JD_launch); % [s]
            % Objective function
            [l, time, T, m, r, flag] = lambLTj_plan(X_1, X_2, V_1, V_2, TOF, N, M_end, n_sol, n_integrator, Is);

            m_cost(i, j) = m(1) - m(end); % Propellant mass required
            T_max(i, j) = max(T); % Max of the thrust history

            if abs(m_cost(i, j)) > 2000 || abs(imag(m_cost(i,j))) > 1e-8 || abs(imag(T_max(i,j))) > 1e-8
                m_cost(i, j) = nan;
                T_max(i, j) = nan;
            end
        else
            m_cost(i, j) = nan;
            T_max(i, j) = nan;

        end
    end

    fprintf('Number of launch dates analyzed: %3.0f out of %3.0f \n', i, length(grid_launch_win));
    toc
end

%%
% m_prop = m(1) - m(end); % First cost

min_mass_prop = min(m_cost, [], 'all');
[rowid , colid] = find(m_cost == min_mass_prop, 1);

T_max_mincost = T_max(rowid, colid);
min_LaunchJD = grid_launch_win(rowid);
min_ArrivalJD = grid_arrival_win(colid);

fprintf('\nMinimum propellant mass found: %3.3f [kg]', min_mass_prop)
fprintf('\nMax required Thrust for optimal cost transfer: %3.3f [N]', T_max_mincost)
fprintf('\nDeparture date of minimum cost transfer: %3.3f', min_LaunchJD)
fprintf('\nArrival date of minimum cost transfer: %3.3f\n', min_ArrivalJD)

surf(grid_launch_win, grid_arrival_win, m_cost)