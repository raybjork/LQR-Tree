function [x_d, u_d, dt] = collocate_trajectory(x_0, x_f, N, system)
%% COLLOCATE_TRAJECTORY(dynamics, x_0, x_f, u_max, N)
%   executes a direct collocation optimization program to find an input
%   sequence to drive the cartpole system from x_0 to x_f.
%
%   @param x_0: the state at the start of the trajectory; n_x by 1 vector
%   @param x_f: the state at the emd of the trajectory; n_x by 1 vector
%
%   @output z: decision variable vector containing the x_i and u_i
%   @output Aeq: matrix from linear constrant Aeq z = beq
%   @output beq: (column) vector from linear constrant Aeq z = beq
%   @output lb: lower bound of constraint lb <= z <= ub; n_z by 1 vector
%   @output ub: upper bound of constraint lb <= z <= ub; n_z by 1 vector
%   @output z0: initial guess for z; n_z by 1 vector
    
    % get information about the system
    u_max = system.u_max;
    f = system.dynamics();

    % constants defining time parameters of the program
    dt_guess = .1;
    dt_max = 1;
    dt_cost = 1;

    % define dimensions of the problem
    nx = length(x_0);
    nu = 1;

    % Add constraints to Aeq, beq to enforce starting at x_0 and
    % ending at x_f
    Aeq = zeros(2*nx, N * (nx + nu) + 1);
    Aeq(1 : nx, 1 : nx) = diag(ones(1, nx));
    Aeq(nx + 1 : end, (N - 1) * (nx + nu) + 1 : (N - 1) * (nx + nu) + nx) = ...
        diag(ones(1, nx));
    beq = [x_0 ; x_f];

    % bounding box constraints u \in [-M,M]^nu
    lb = -inf(N * (nx + nu) + 1, 1);
    ub = inf(N * (nx + nu) + 1, 1);
    for i=1:N
        % bounding box constraints for u_i
        [~,u_i_inds] = sample_indices(i, nx, nu);
        lb(u_i_inds) = -u_max;
        ub(u_i_inds) = u_max;
    end
    lb(end) = 0; % lower bound for time step
    ub(end) = dt_max; % upper bound for time step


    % make initial guess for z
    z0 = reshape([interp1([1 2],[x_0'; x_f'], linspace(1, 2, N))'; ...
        zeros(1, N)], [N * (nx + nu), 1]);
    z0 = [z0; dt_guess];
    
    options = optimoptions('fmincon','SpecifyObjectiveGradient',false,...
        'SpecifyConstraintGradient', false,'Display','iter');
    problem.objective = @(z) trajectory_cost(z, N, nx, nu, dt_cost);

    problem.x0 = z0;
    problem.options = options;
    problem.nonlcon = @(z) all_constraints(f, z, N, nx, nu);
    problem.solver = 'fmincon';
    problem.Aeq = Aeq;
    problem.beq = beq;
    problem.lb = lb;
    problem.ub = ub;

    z = fmincon(problem);
    
    dt = z(end);
    u = zeros(nu,0);
    x_d = zeros(nx,0);
    for i=1:N
       [x_i_inds,u_i_inds] = sample_indices(i, nx, nu);
       x_d(:,i) = z(x_i_inds);
       u(:,i) = z(u_i_inds);
    end
    u_d = interp1(0:dt:(dt*(N-1)), u', 'linear', 'pp');
    x_d = spline(0:dt:(dt*(N-1)), x_d);
    u_d = @(t) ppval(u_d, t);
    x_d = @(t) ppval(x_d, t);
end

function [c, ceq, dC, dCeq] = all_constraints(dynamics, z, N, nx, nu)
    [ceq, dCeq] = dynamics_constraints(dynamics, z, N, nx, nu);
    c = zeros(0,1);
    dC = zeros(0,numel(z));
    dC = sparse(dC)';
    dCeq = sparse(dCeq)';
end

function [x_i_inds, u_i_inds] = sample_indices(i, nx, nu)
%% SAMPLE_INDICES 
%   calculates indices of z such that z(x_i_inds)=x_i and z(u_i_inds)=u_i.
%   @param i: sample number; scalar
%   @param nx: dimension of state vector, x; scalar
%   @param nu: dimension of input vector, u; scalar
%
%   @output x_i_inds: indices such that z(x_i_inds) = x_i; 1 by n_x vector
%   @output u_i_inds: indices such that z(u_i_inds) = u_i; 1 by n_u vector
    x_i_inds = (i-1)*(nx+nu) + 1: (i-1)*(nx+nu) + nx;
    u_i_inds = (i-1)*(nx+nu) + 1 + nx;
end

function [g,dG] = trajectory_cost(z, N, nx, nu, dt_cost)
%% TRAJECTORY_COST(z) 
%   computes the cost and cost jacobian.
%   @param z: decision variable (column) vector containing the x_i and u_i
%   @param N: number of sample points; scalar
%   @param nx: dimension of state vector, x; scalar
%   @param nu: dimension of input vector, u; scalar
%   @param dt: \Delta t, the inter-sample interval duration; scalar
%
%   @output g: total accrued cost; scalar
%   @output dG: gradient of total accrued cost; nz by 1 vector

    g = 0;
    dG = zeros(N*(nx + nu) + 1, 1);
    for i=1:(N-1)
        [~, u_i] = sample_indices(i, nx, nu);
        [~, u_ip1] = sample_indices(i+1, nx, nu);
        g = g + ((z(u_i)^2) + (z(u_ip1)^2))*(z(end) / 2);
        dG(u_i) = 2 * z(u_i) * z(end);
    end
    g = g + dt_cost * z(end);
    [~, u_1] = sample_indices(1, nx, nu);
    [~, u_n] = sample_indices(N, nx, nu);
    dG(u_1) = z(u_1) * z(end);
    dG(u_n) = z(u_n) * z(end);
    dG(end) = g/z(end) + dt_cost;
end

function [h,dH] = dynamics_constraints(dynamics, z, N, nx, nu)
%% DYNAMICS_CONSTRAINTS(z) 
%   compiles the dynamics constraints generated by dynamics_constraint_with_derivative.
%   @param z: decision variable (column) vector containing the x_i and u_i
%   @param N: number of sample points; scalar
%   @param nx: dimension of state vector, x; scalar
%   @param nu: dimension of input vector, u; scalar
%   @param dt: \Delta t, the inter-sample interval duration; scalar

%   @output h: compiled h_i from dynamics_constraint_with_derivative;
%   (N-1)*nx by 1 vector
%   @output dH_i: compiled dH_i from dynamics_constraint_with_derivative;
%   (N-1)*nx by nz matrix

    h = zeros((N-1)*nx, 1);
    dH = zeros((N-1)*nx, N*(nx + nu) + 1);
    for i=1:(N-1)
        start = (i-1)*(nx*nu) + 1;
        [x_i, u_i] = sample_indices(i, nx, nu);
        [x_ip1, u_ip1] = sample_indices(i + 1, nx, nu);
        [h_i, dH_i] = dynamics_constraint_with_derivative(dynamics, ...
            z(x_i), z(u_i), z(x_ip1), z(u_ip1), z(end));
        h(start: start + nx - 1) = h_i;
        dH(start: start + nx - 1, (i-1)*(nx+nu) + 1: (i+1)*(nx+nu)) = dH_i(:,1:end-1);
        dH(start: start + nx - 1, end) = dH_i(:,end);
    end
end

function [h_i, dH_i] = dynamics_constraint_with_derivative(dynamics, x_i, u_i, x_ip1, u_ip1, dt)
%% DYNAMICS_CONSTRAINT_WITH_DERIVATIVE(x_i, u_i, x_ip1, u_ip1, dt) 
%   returns and computes the gradient of the vector constraint asssociated
%   with dynamics_constraint(x_i, u_i, x_ip1, u_ip1, dt).
%
%   @param x_i: the state at the start of the interval; nx by 1 vector
%   @param u_i: the input at the start of the interval; nu by 1 vector
%   @param x_ip1: the state at the end of the interval; nx by 1 vector
%   @param u_ip1: the input at the end of the interval; nu by 1 vector
%   @param dt: \Delta t, the duration of the interval; scalar
%
%   @output h_i: constraint value from dynamics_constraint; nx by 1 vector
%   @output dH_i: jacobian of h_i w.r.t. [x_i; u_i; x_ip1; u_ip1]; nx by
%   (2nx + 2nu) matrix

    h_i = evaluate_dynamics_constraint(dynamics, x_i, u_i, x_ip1, u_ip1, dt);
    if nargout > 1
      % use numerical derivatives to compute dH
      % dH = [dh/dx0 dh/du0 dh/dx1 dh/du1]
      % where the partial derivatives are written (dh/dx0)_ij = dh_i/dx0_j
      delta = 1e-8;
      dH_i = zeros(numel(x_i), 2*(numel(x_i)+numel(u_i))+1);
      for j=1:numel(x_i)
          dx = zeros(numel(x_i),1);
          dx(j) = delta;
          dHx_i_j = evaluate_dynamics_constraint(dynamics, x_i + dx, u_i, x_ip1, u_ip1, dt) - h_i;
          dHx_ip1_j = evaluate_dynamics_constraint(dynamics, x_i, u_i, x_ip1 + dx, u_ip1, dt) - h_i;
          dH_i(:,j) = dHx_i_j/delta;
          dH_i(:,j + numel(x_i) + numel(u_i)) = dHx_ip1_j/delta;
      end

      for j=1:numel(u_i)
          du = zeros(numel(u_i),1);
          du(j) = delta;
          dHu_i_j = evaluate_dynamics_constraint(dynamics, x_i, u_i + du, x_ip1, u_ip1, dt) - h_i;
          dHu_ip1_j = evaluate_dynamics_constraint(dynamics, x_i, u_i, x_ip1, u_ip1 + du, dt) - h_i;
          dH_i(:,j + numel(x_i)) = dHu_i_j/delta;
          dH_i(:,j + numel(x_i) + numel(u_i) + numel(x_ip1)) = dHu_ip1_j/delta;
      end
      
      dHu_i_j = evaluate_dynamics_constraint(dynamics, x_i, u_i, x_ip1, u_ip1, dt + delta) - h_i;
      dH_i(:,end) = dHu_i_j/delta;      
    end
end

function h_i = evaluate_dynamics_constraint(f, x_i, u_i, x_ip1, u_ip1, dt)
%% EVALUATE_DYNAMICS_CONSTRAINT(x_i, u_i, x_ip1, u_ip1, dt)
%   computes the vector contstraint h_i(x_i, u_i, x_ip1, u_ip1) = 0
%
%   @param x_i: the state at the start of the interval; nx by 1 vector
%   @param u_i: the input at the start of the interval; nu by 1 vector
%   @param x_ip1: the state at the end of the interval; nx by 1 vector
%   @param u_ip1: the input at the end of the interval; nu by 1 vector
%   @param dt: \Delta t, the duration of the interval; scalar
%
%   @output h_i: quantity derived in Problem 2(c); nx by 1 vector

    s_midway = .5*(x_i + x_ip1)-(dt/8)*(f(x_ip1, u_ip1) - f(x_i, u_i));
    s_dot_midway = (1.5/dt)*(x_ip1 - x_i)-.25*(f(x_ip1, u_ip1) + f(x_i, u_i));
    h_i = s_dot_midway - f(s_midway, .5*(u_i + u_ip1));
end