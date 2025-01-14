alp = 2;
beta = 0.001;
gam = 20;
delt = 0.002;

uf = @(u, v, alp, beta) u * (alp - beta * v);
vf = @(u, v, gam, delt) -1 * v * (gam - delt * u);

time = 0;
fc = [1/4, 3/8, 12/13, 1, 1/2];
fal = {[1/5], [3/40, 9/40], [3/10, -9/10 , 6/5], ...
       [-11/54    ,       5/2      ,    -70/27   ,       35/27], ...
       [1631/55296 ,   175/512    ,   575/13824   ,  44275/110592   ,  253/4096] ...
       };

fb = [37/378          ,  0     ,     250/621     ,    125/594     ,          0 ,512/1771 ];
fb_ = [2825/27648    ,    0     ,    18575/48384 ,    13525/55296   ,   277/14336  ,1/4];

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
    
   

    % Combine the k-values for the next step
    k = [k1, k2, k3, k4, k5, k6];
    k_ = [k1_, k2_, k3_, k4_, k5_, k6_];

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