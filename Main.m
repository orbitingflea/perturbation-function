function Main()
    gamma = 0.1;
    e = exp(1);
    Gamma = 1 - 1/e - gamma;
    delta = 1;
    R = -log(1 - Gamma);
    assert(0 <= R && R <= 1);

    function v = fLowerBound(x)
        bound1 = Bound1(x, Gamma);
        bound3 = Bound3(x, Gamma, delta);
        v = max(bound1, bound3);
    end

    fprintf("[INFO] show some bounds:\n");
    for x = 0 : 0.05 : 1
        bound1 = Bound1(x, Gamma);
        bound3 = Bound3(x, Gamma, delta);
        fprintf("x = %.10f,  bound1 = %.10f,  bound3 = %.10f,  difference = %.10f\n", ...
                x, bound1, bound3, bound3 - bound1);
    end
    F_upper = UpperBound2(R, Gamma);
    fprintf("Upper bound: %.10f\n", F_upper);

    flb = @(x) arrayfun(@fLowerBound, x);
    int_lower_bound = integral(flb, 0, R);
    fprintf("Lower bound (integral): %.10f\n", int_lower_bound);

    if int_lower_bound > F_upper
        fprintf("[CONCLUSION] Good.\n");
    else
        fprintf("[CONCLUSION] Failed.\n");
    end
end


function v = FunInner(theta, y, Gamma, delta)
    branch1 = 1 - exp(y) .* (1 - Gamma);
    branch2 = (1 + delta) ./ Gamma .* (1 - exp(theta) .* (1 - Gamma)) .* ...
              (1 - exp(y - theta) .* (1 - Gamma));
    v = max(0, min(branch1, branch2));
end


function v = Bound3(theta, Gamma, delta)
    e = exp(1);
    gamma = 1 - (1 / e) - Gamma;
    R = -log(1 - Gamma);
    if (theta > R)
        v = -inf;
        return;
    end

    term1 = (Gamma - theta) ./ (1 - Gamma) .* ...
            (theta - (1 - Gamma) .* (exp(theta) - 1));
    term2_coeff = 1 / (1 - Gamma) .* ...
        (Gamma - (2 * e) / (e - 1) .* (2 ./ delta + 1) .* gamma);

    function v = InnerFunOfY(y)
        v = FunInner(theta, y, Gamma, delta);
    end

    term2_int = integral(@InnerFunOfY, theta, R);
    v = term1 + term2_coeff .* term2_int;
    % fprintf("term1: %f, term2_coeff: %f, term2_int: %f\n", term1, term2_coeff, term2_int);
end


function v = Bound1(theta, Gamma)
    v = (1 - exp(theta) + Gamma * exp(theta));
end


function v = UpperBound2(alph, Gamma)
    v = (alph + 1 - exp(alph - 1) - Gamma);
end
