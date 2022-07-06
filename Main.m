function max_dist = Main(gamma, delta, info)
    if nargin < 3
        info = true;
    end

    e = exp(1);
    Gamma = 1 - 1/e - gamma;
    assert(0.5 <= Gamma && Gamma <= 1 - 1/e);
    assert(0 < delta);

    % R is the root of (1 - e^R * (1 - Gamma)).
    % R is denoted as lowercase r in the paper.
    R = -log(1 - Gamma);
    assert(0 <= R && R <= 1);

    % Compute BetaStar.
    beta_star = GetBetaStar(Gamma, delta, info);
    if beta_star < 0
        % Cannot compute beta_star. Exit.
        max_dist = -inf;
        if info
            fprintf("[INFO] gamma = %.10f, delta = %.10f is not valid.\n", gamma, delta);
        end
        return;
    end

    if info
        fprintf("[INFO] start with delta = %f, gamma = %f\n", delta, gamma);
    end

    function v = fLowerBound(alpha)
        % Obtain the better lower bound for f(alpha) among the two bounds.
        bound_simple = LowerBoundSimple(alpha, Gamma);
        bound_h = LowerBoundH(alpha, Gamma, delta, beta_star);
        v = max(bound_simple, bound_h);
    end

    if info
        fprintf("[INFO] show some bounds:\n");
    end

    % max_dist is the (approximate) maximum value of bound_h - bound_simple.
    % It is not necessary for verifying the result.
    % It can be used for debugging / finding a good delta for fixed Gamma.
    % See Opt.m for more details about optimization.
    max_dist = -inf;
    for x = 0 : 0.05 : 1
        bound_simple = LowerBoundSimple(x, Gamma);
        bound_h = LowerBoundH(x, Gamma, delta, beta_star, info);
        if info
            fprintf("[DEBUG] x = %.10f, bound_simple = %.10f, bound_h = %.10f, " ...
                    + "difference = %.10f\n", x, bound_simple, bound_h, bound_h - bound_simple);
        end
        max_dist = max(max_dist, bound_h - bound_simple);
    end

    % Obtain the lower and upper bounds for integral of f on [0, R],
    % and try to make a contradiction.
    int_f_upper_bound = IntegralFUpperBound(R, Gamma);
    fLowerBoundVec = @(x) arrayfun(@fLowerBound, x);
    int_f_lower_bound = integral(fLowerBoundVec, 0, R);

    if info
        fprintf("[DEBUG] max_dist: %.10f\n", max_dist);
        fprintf("\n");
        fprintf("Program finished.\n  delta = %.10f\n  R = %.10f\n", delta, R);
        fprintf("  beta_star = %.10f\n", beta_star);
        fprintf("  Lower bound of integral of f on [0, R]: %.10f\n", int_f_lower_bound);
        fprintf("  Upper bound of integral of f on [0, R]: %.10f\n", int_f_upper_bound);
        if int_f_lower_bound > int_f_upper_bound
            fprintf("[CONCLUSION] Gamma = (1 - 1/e) - (%.10f) is not achievable.\n", gamma);
            fprintf("  This is because for the specified Gamma, the lower bound is greater than\n");
            fprintf("  the upper bound, which leads to a contradiction.\n");
        else
            fprintf("[CONCLUSION] Failed to make a contradiction.\n");
        end
    end
end


function v = LowerBoundH(alpha, Gamma, delta, beta_star, info)
    % Recall that in the paper we obtained the inequality f(alpha) >= h(alpha).
    % This function computes h(alpha) which is a lower bound of f.
    if nargin < 5
        info = false;
    end

    e = exp(1);
    gamma = 1 - (1 / e) - Gamma;
    R = -log(1 - Gamma);
    if (alpha > R)
        v = -inf;
        return;
    end

    term1 = (Gamma - alpha) / (1 - Gamma) * ...
            (alpha - (1 - Gamma) * (exp(alpha) - 1));
    term2_coeff = (Gamma - beta_star) / (1 - Gamma);

    function v = InnerFunction(x)
        % InnerFunction is the integrand in the second term of h(alpha).
        branch1 = 1 - exp(x) * (1 - Gamma);
        branch2 = (1 - exp(alpha) * (1 - Gamma)) / ((1 + delta) * Gamma) * ...
                  (1 - exp(x - alpha) * (1 - Gamma));
        v = min(branch1, branch2);
        assert(all(branch2 >= 0));
    end

    InnerFunctionVec = @(x) arrayfun(@InnerFunction, x);
    term2_int = integral(InnerFunctionVec, alpha, R);
    v = term1 + term2_coeff * term2_int;
    if info
        % fprintf("[DEBUG] LowerBoundH: term1: %f, term2_coeff: %f, term2_int: %f\n", ...
        %         term1, term2_coeff, term2_int);
    end
end


function v = LowerBoundSimple(alpha, Gamma)
    % This function is the simpler lower bound for f(alpha).
    v = 1 - exp(alpha) * (1 - Gamma);
end


function v = IntegralFUpperBound(beta, Gamma)
    % The upper bound for integral of f(x) on [0, beta].
    v = (beta + 1 - exp(beta - 1) - Gamma);
end


function beta_star = GetBetaStar(Gamma, delta, info)
    % This function computes beta_star, which is the unique root of
    % a function defined in the paper.

    if nargin < 3
        info = true;
    end

    e = exp(1);
    gamma = (1 - 1/e) - Gamma;

    function val = SubFunction(beta)
        % SubFunction is the function of beta whose unique root defines
        % beta_star.
        val = (delta + e * gamma + delta * e * gamma) * ...
              (exp(beta - 1) - 1/e) + gamma - beta * delta;
    end

    % We have SubFunction(0) == gamma.
    assert(abs(SubFunction(0) - gamma) <= 1e-8);

    if SubFunction(1) >= 0
        % SubFunction does not have a unique root in [0, 1].
        % Our method cannot obtain a valid bound. Exit.
        beta_star = -1;
        return;
    end

    % Find the unique root in [0, 1].
    root = fzero(@SubFunction, [0, 1]);
    if root >= Gamma
        % The unique root of SubFunction does not satisfy [root < Gamma].
        % Our method cannot obtain a valid bound. Exit.
        beta_star = -1;
        return;
    end

    beta_star = root;

    if info
        % fprintf("[DEBUG] beta_star: %.10f\n", beta_star);
        % fprintf("[DEBUG] subfunction(1): %.10f\n", SubFunction(1));
    end
end
