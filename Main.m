function max_dist = Main(gamma, delta, info)
    if nargin < 3
        info = true;
    end

    e = exp(1);
    Gamma = 1 - 1/e - gamma;
    R = -log(1 - Gamma);
    assert(0.5 <= Gamma && Gamma <= 1 - 1/e);
    assert(0 < delta);
    assert(0 <= R && R <= 1);

    if info
        fprintf("[INFO] start with delta = %f, gamma = %f\n", delta, gamma);
    end

    function v = fLowerBound(x)
        bound1 = Bound1(x, Gamma);
        bound3 = Bound3(x, Gamma, delta);
        v = max(bound1, bound3);
    end

    if info
        fprintf("[INFO] show some bounds:\n");
    end
    max_dist = -inf;
    for x = 0 : 0.05 : 1
        bound1 = Bound1(x, Gamma);
        bound3 = Bound3(x, Gamma, delta, info);
        if info
            fprintf("x = %.10f,  bound1 = %.10f,  bound3 = %.10f,  difference = %.10f\n", ...
                    x, bound1, bound3, bound3 - bound1);
        end
        max_dist = max(max_dist, bound3 - bound1);
    end
    F_upper = UpperBound2(R, Gamma);
    flb = @(x) arrayfun(@fLowerBound, x);
    int_lower_bound = integral(flb, 0, R);

    if info
        fprintf("\n");
        fprintf("Program finished.\ndelta = %.10f\nR = %.10f\n", delta, R);
        fprintf("Lower bound integral to R: %.10f\n", int_lower_bound);
        fprintf("Upper bound of F(R): %.10f\n", F_upper);
        fprintf("Max distance between bound3 and bound1: %.10f\n", max_dist);
        if int_lower_bound > F_upper
            fprintf("[CONCLUSION] Good.\n");
        else
            fprintf("[CONCLUSION] Failed.\n");
        end
    end
end


function v = FunInner(theta, y, Gamma, delta)
    branch1 = 1 - exp(y) .* (1 - Gamma);
    branch2 = 1 ./ (1 + delta) ./ Gamma .* (1 - exp(theta) .* (1 - Gamma)) .* ...
              (1 - exp(y - theta) .* (1 - Gamma));
    v = max(0, min(branch1, branch2));
end


function v = Bound3(theta, Gamma, delta, info)
    if nargin < 4
        info = false;
    end

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
    if info
        fprintf("term1: %f, term2_coeff: %f, term2_int: %f\n", term1, term2_coeff, term2_int);
    end
end


function v = Bound1(theta, Gamma)
    v = (1 - exp(theta) + Gamma * exp(theta));
end


function v = UpperBound2(alph, Gamma)
    v = (alph + 1 - exp(alph - 1) - Gamma);
end
