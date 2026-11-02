function [y,history,iterations,stop] = solve_gcv_blur(y,bhat,s,lambda,beta,gamma,s_aux,symbol,valid,eta)
% Modified Newton step at fixed lambda during conditional-GCV alternation.
% Use precision units eta*f for the fixed gradient tolerance and damping.
% The inner solve allows at most 50 modified Newton steps.
y=validate_solver_inputs(y,bhat,s,lambda,beta,gamma,s_aux,symbol,valid);
validateattributes(eta,{'numeric'},{'real','scalar','finite','positive'},mfilename,'eta');
objective_scale = eta;
step_tolerance = 1e-12;
history = zeros(numel(y),50);
stop = 0;
[f,g,B] = reconstruction_objective(y,bhat,s,lambda,beta,gamma,s_aux,symbol,valid);
assert(isreal(f) && isreal(g) && isreal(B) && isfinite(f) && ...
    all(isfinite(g(:))) && all(isfinite(B(:))), ...
    'MargRegVarPro:InvalidInitialObjective', ...
    'The initial objective, gradient, and curvature matrix must be real and finite.');
for iterations = 1:50
    history(:,iterations) = y;
    gradient = objective_scale*g;
    if norm(gradient) < 1e-3
        stop = 1;
        break
    end
    matrix = objective_scale*B;
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
        if isfinite(f_trial) && objective_scale*f_trial <= objective_scale*f + 1e-4*t*(gradient.'*direction)
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

