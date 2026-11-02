function [y,history,iterations,stop] = solve_reconstruction(y,bhat,s,lambda,beta,gamma,s_aux,symbol,valid,step_tolerance)
% Modified Newton method with Armijo backtracking, for this fixed-parameter study.
% Evaluate f_{beta,gamma} directly, in the manuscript residual units.
% Gradient tolerance and diagonal damping are scaled by 1/N.
y=validate_solver_inputs(y,bhat,s,lambda,beta,gamma,s_aux,symbol,valid);
N = numel(bhat);
gradient_tolerance = 1e-3/N;
initial_damping = 1e-6/N;
minimum_damping = 1e-16/N;
maximum_damping = 1e8/N;
if nargin < 10, step_tolerance = 1e-12; end
validateattributes(step_tolerance,{'numeric'},{'real','scalar','finite','positive'},mfilename,'step_tolerance');
history = zeros(numel(y),500);
stop = 0;
[f,g,B] = reconstruction_objective(y,bhat,s,lambda,beta,gamma,s_aux,symbol,valid);
assert(isreal(f) && isreal(g) && isreal(B) && isfinite(f) && ...
    all(isfinite(g(:))) && all(isfinite(B(:))), ...
    'MargRegVarPro:InvalidInitialObjective', ...
    'The initial objective, gradient, and curvature matrix must be real and finite.');
for iterations = 1:500
    history(:,iterations) = y;
    if norm(g) < gradient_tolerance
        stop = 1;
        break
    end
    matrix = B;
    matrix = (matrix+matrix.')/2;
    damping = initial_damping;
    [R,flag] = chol(matrix);
    while flag ~= 0
        matrix = matrix + damping*eye(numel(y));
        damping = max(10*damping,minimum_damping);
        [R,flag] = chol(matrix);
        if damping > maximum_damping, break; end
    end
    if flag == 0
        direction = -R \ (R' \ g);
    else
        % The fallback direction uses the N-scaled gradient.
        direction = -N*g;
    end

    t = .25;
    while true
        y_trial = y + t*direction;
        if norm(y_trial-y) < 1e-16
            t = .5*t;
            if t < 1e-14
                y_trial=y; f_trial=f; g_trial=g; B_trial=B;
                break
            end
            continue
        end
        [f_trial,g_trial,B_trial] = reconstruction_objective( ...
            y_trial,bhat,s,lambda,beta,gamma,s_aux,symbol,valid);
        if isfinite(f_trial) && f_trial <= f + 1e-4*t*(g.'*direction)
            break
        end
        t = .5*t;
        if t < 1e-14
            y_trial=y; f_trial=f; g_trial=g; B_trial=B;
            break
        end
    end
    step = norm(y_trial-y);
    y=y_trial; f=f_trial; g=g_trial; B=B_trial;
    if step < step_tolerance
        stop = 2;
        break
    end
end
history = history(:,1:iterations);
end
