% Prova plot interplanetray for relative positions
clearvars; close all; clc

[R_PROP,TIMES,FLAGS, r_DSM, r_planets] = plot_trajectory(2);

% find WHEN DSM and FLYBYS occur
save_DSM = zeros(1, length(r_DSM));
save_pl = zeros(1, length(r_planets));

c = 0;
for j = 1:length(R_PROP)
    for i = 1:length(r_DSM)
        check = r_DSM(:,i) - R_PROP(:,j);
        if abs(check) < 1e-12*ones(3,1)
            c = c + 1;
            save_DSM(c) = j; 
        end
    end
end

save_pl(1) = 1;
c = 1;
for j = 2:length(R_PROP)
    for i = 2:length(r_planets)
        check = r_planets(:,i) - R_PROP(:,j);
        if abs(check) < 1e-12*ones(3,1)
            c = c + 1;
            save_pl(c) = j; 
        end
    end
end

ind = zeros(1,9);
for i = 1:8
    if rem(i,2) ~= 0
        ind(i) = save_pl(i);
    else
        ind(i) = save_DSM(i);
    end
end
ind(9) = save_pl(end);

% total duration in years
tof_years = TIMES(end)-TIMES(1);
fprintf('\nOverall Transfer: %.2f years\n', tof_years);

% number of weeks
weeks = round((tof_years*365.25)/7);

% Elapsed time in hour
TIMES = TIMES*365.25*24;
tof_h = zeros(1,8);
tof_d = zeros(1,8);
tof_y = zeros(1,8);
week_each = zeros(1,8);
for i = 2:length(ind)
    t2 = TIMES(ind(i));
    t1 = TIMES(ind(i-1));

    tof_h(i-1) = t2 - t1;          % tof in hour
    tof_d(i-1) = tof_h(i-1)/24;      % tof in days
    tof_y(i-1) = tof_d(i-1)/365.25;  % tof in years

    week_each(i-1) = round(tof_d(i-1)/7);
end
week_all = sum(week_each);
