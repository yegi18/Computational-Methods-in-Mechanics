% Constants and coefficients
alp = 2;
beta = 0.001;
gam = 20;
delt = 0.002;

% Time variables
dt = 0.01;
time = 0;

% Initial conditions
u = [5000, 5000.249865174713, 5000.499746719312, 5000.749649394678];
v = [100, 99.96587343229123, 99.93175987764373, 99.89765868126679];
rel_ = [];  % Store relative errors
time_steps = [];  % Store time points for plotting

% Define the ODEs
uf = @(u, v, alp, beta) u * (alp - beta * v);  % du/dt
vf = @(u, v, gam, delt) -v * (gam - delt * u);  % dv/dt

% Coefficients
b = {[1], [0, 1], [-3/2, 3, -1/2], [-10/3, 6, -2, 1/3]};
a = [1, 2, 3, 4];
a_ = [1, 2/3, 6/11, 12/25];
b_ = {[1], [4/3, -1/3], [18/11, -9/11, 2/11], [48/25, -36/25, 16/25, -3/25]};

max_value = 1e10;  % Set a maximum threshold for u and v

% Main time loop
while time <= 100
    k = floor(time / dt);
    
    if k >= 4  % Apply multistep method after enough time steps
        % Get the last 4 values for predictor and corrector steps
        uvec = flip(u(end-3:end));  % Use the last 4 values for u
        vvec = flip(v(end-3:end));  % Use the last 4 values for v

        % Predictor step
        unew_pred = dot(b{4}, uvec) + a(4) * dt * uf(u(end), v(end), alp, beta);
        vnew_pred = dot(b{4}, vvec) + a(4) * dt * vf(u(end), v(end), gam, delt);

        % Corrector step
        unew = dot(b_{4}, uvec) + a_(4) * dt * uf(unew_pred, vnew_pred, alp, beta);
        vnew = dot(b_{4}, vvec) + a_(4) * dt * vf(unew_pred, vnew_pred, gam, delt);

        % Calculate relative errors
        relu = abs((unew - unew_pred) / unew);
        relv = abs((vnew - vnew_pred) / vnew);
        rel_ = [rel_, max(relu, relv)];
        time_steps = [time_steps, time];

        % Append new values of u and v
        u = [u, unew];
        v = [v, vnew];
    end

    time = time + dt;
end

% Plot u and v over time
figure;
plot(time_steps, u(5:end), 'b-', 'LineWidth', 1.5);  % Plotting u starting from 5th value
hold on;
plot(time_steps, v(5:end), 'r-', 'LineWidth', 1.5);  % Plotting v starting from 5th value

xlabel('Time');
ylabel('Values of u and v');
title('u and v over Time with Time Step = 0.01');
legend('u', 'v');
hold off;

% Plot the relative error over time
figure;
plot(time_steps, rel_, 'g-', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Relative Error');
title('Relative Error over Time with Time Step = 0.01');
