alp = 2;
beta = 0.001;
gam = 20;
delt = 0.002;

uf = @(u, v, alp, beta) u * (alp - beta * v);
vf = @(u, v, gam, delt) -1 * v * (gam - delt * u);

time = 0;
fc = [1/4, 3/8, 12/13, 1, 1/2];
fal = {[1/5], [3/40, 9/40], [44/45, -56/15, 32/9], ...
       [19372/6561, -25360/2187, 64448/6561, -212/729], ...
       [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656], ...
       [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]};

fb = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
fb_ = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];

uold = 5000;
vold = 100;
dt = 0.01;
u = [];
v = [];

t_end = 100;  % Simulation end time
tol = 1e-4;   % Tolerance for error control

tic;  % Start timing

while time <= t_end
    % Dormand-Prince step calculation
    k1 = uf(uold, vold, alp, beta);
    k1_ = vf(uold, vold, gam, delt);
    
    k2 = uf(uold + k1 * fal{1}(1) * dt, vold + k1_ * fal{1}(1) * dt, alp, beta);
    k2_ = vf(uold + k1 * fal{1}(1) * dt, vold + k1_ * fal{1}(1) * dt, gam, delt);
    
    k3 = uf(uold + k1 * fal{2}(1) * dt + k2 * fal{2}(2) * dt, vold + k1_ * fal{2}(1) * dt + k2_ * fal{2}(2) * dt, alp, beta);
    k3_ = vf(uold + k1 * fal{2}(1) * dt + k2 * fal{2}(2) * dt, vold + k1_ * fal{2}(1) * dt + k2_ * fal{2}(2) * dt, gam, delt);
    
    k4 = uf(uold + k1 * fal{3}(1) * dt + k2 * fal{3}(2) * dt + k3 * fal{3}(3) * dt, vold + k1_ * fal{3}(1) * dt + k2_ * fal{3}(2) * dt + k3_ * fal{3}(3) * dt, alp, beta);
    k4_ = vf(uold + k1 * fal{3}(1) * dt + k2 * fal{3}(2) * dt + k3 * fal{3}(3) * dt, vold + k1_ * fal{3}(1) * dt + k2_ * fal{3}(2) * dt + k3_ * fal{3}(3) * dt, gam, delt);
    
    k5 = uf(uold + k1 * fal{4}(1) * dt + k2 * fal{4}(2) * dt + k3 * fal{4}(3) * dt + k4 * fal{4}(4) * dt, vold + k1_ * fal{4}(1) * dt + k2_ * fal{4}(2) * dt + k3_ * fal{4}(3) * dt + k4_ * fal{4}(4) * dt, alp, beta);
    k5_ = vf(uold + k1 * fal{4}(1) * dt + k2 * fal{4}(2) * dt + k3 * fal{4}(3) * dt + k4 * fal{4}(4) * dt, vold + k1_ * fal{4}(1) * dt + k2_ * fal{4}(2) * dt + k3_ * fal{4}(3) * dt + k4_ * fal{4}(4) * dt, gam, delt);
    
    k6 = uf(uold + k1 * fal{5}(1) * dt + k2 * fal{5}(2) * dt + k3 * fal{5}(3) * dt + k4 * fal{5}(4) * dt + k5 * fal{5}(5)*dt, vold + k1_ * fal{5}(1) * dt + k2_ * fal{5}(2) * dt + k3_ * fal{5}(3) * dt + k4_ * fal{5}(4) * dt + k5_ * fal{5}(5)*dt, alp, beta);
    k6_ = vf(uold + k1 * fal{5}(1) * dt + k2 * fal{5}(2) * dt + k3 * fal{5}(3) * dt + k4 * fal{5}(4) * dt + k5 * fal{5}(5)*dt, vold + k1_ * fal{5}(1) * dt + k2_ * fal{5}(2) * dt + k3_ * fal{5}(3) * dt + k4_ * fal{5}(4) * dt + k5_ * fal{5}(5)*dt, gam, delt);
    
    k7 = uf(uold + k1 * fal{6}(1) * dt + k2 * fal{6}(2) * dt + k3 * fal{6}(3) * dt + k4 * fal{6}(4) * dt + k5 * fal{6}(5)*dt + k6 * fal{6}(6)*dt, vold + k1_ * fal{6}(1) * dt + k2_ * fal{6}(2) * dt + k3_ * fal{6}(3) * dt + k4_ * fal{6}(4) * dt + k5_ * fal{6}(5)*dt + k6_ * fal{6}(6)*dt, alp, beta);
    k7_ = vf(uold + k1 * fal{6}(1) * dt + k2 * fal{6}(2) * dt + k3 * fal{6}(3) * dt + k4 * fal{6}(4) * dt + k5 * fal{6}(5)*dt + k6 * fal{6}(6)*dt, vold + k1_ * fal{6}(1) * dt + k2_ * fal{6}(2) * dt + k3_ * fal{6}(3) * dt + k4_ * fal{6}(4) * dt + k5_ * fal{6}(5)*dt + k6_ * fal{6}(6)*dt, gam, delt);

    % Combine the k-values for the next step
    k = [k1, k2, k3, k4, k5, k6, k7];
    k_ = [k1_, k2_, k3_, k4_, k5_, k6_, k7_];

    % Calculate the next values of u and v using fb and fb_
    unew = uold + dt * sum(fb .* k);
    vnew = vold + dt * sum(fb .* k_);
    unew_ = uold + dt * sum(fb_ .* k);
    vnew_ = vold + dt * sum(fb_ .* k_);

    % Relative error estimates for adaptive time step
    relu = abs((unew - unew_) / (unew + 1e-10));
    relv = abs((vnew - vnew_) / (vnew + 1e-10));
    
    % Calculate the local error and adjust time step
    error_estimate = max(relu, relv);
    if error_estimate > tol
        dt = dt * (tol / error_estimate) ^ 0.2;
    else
        dt = dt * (tol / error_estimate) ^ 0.2;
        time = time + dt;
        u = [u, unew]; %#ok<AGROW>
        v = [v, vnew]; %#ok<AGROW>
        uold = unew;
        vold = vnew;
    end
    
    fprintf('Time: %.6f, u: %.6f, v: %.6f, dt: %.6e, error: %.2e\n', time, unew, vnew, dt, error_estimate);
end

toc;  % End timing
figure;
plot(u)