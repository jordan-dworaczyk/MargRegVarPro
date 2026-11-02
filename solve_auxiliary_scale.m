function [y,history,iterations,stop] = solve_auxiliary_scale(y,bhat,s,lambda,beta,gamma,s_aux,symbol,valid,step_tolerance)
% Modified Newton method for the auxiliary-scale sensitivity experiment.
% Evaluate the objective, gradient, and curvature as N*f, N*g, and N*B.
% This gives gradient tolerance 1e-3 and initial diagonal damping 1e-6.
y=validate_solver_inputs(y,bhat,s,lambda,beta,gamma,s_aux,symbol,valid);
N = numel(bhat);
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
    gradient = N*g;
    if norm(gradient) < 1e-3
        stop = 1;
        break
    end
    matrix = N*B;
    matrix = (matrix+matrix.')/2;
    damping = 1e-6;
    [R,flag] = chol(matrix);
    while flag ~= 0
        matrix = matrix + damping*eye(numel(y));
        damping = max(10*damping,1e-16);
        [R,flag] = chol(matrix);
        if damping > 1e8, break; end
    end
    if flag == 0
        direction = -R \ (R' \ gradient);
    else
        direction = -gradient;
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
        if isfinite(f_trial) && N*f_trial <= N*f + 1e-4*t*(gradient.'*direction)
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
