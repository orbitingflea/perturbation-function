function Main()
    gamma = 1e-10;
    e = exp(1);
    Gamma = 1 - 1/e - gamma;
    delta = 1e-5;

    function v = fLowerBound(x)
        bound1 = Bound1(x, Gamma);
        bound3 = Bound3(x, Gamma, delta);
        v = max(bound1, bound3);
    end

    fprintf("[INFO] show some bounds:\n");
    for x = 0 : 0.01 : 1
        bound1 = Bound1(x, Gamma);
        bound3 = Bound3(x, Gamma, delta);
        fprintf("x=%.10f  %.10f  %.10f  delta = %.10f\n", ...
                x, bound1, bound3, bound3 - bound1);
    end
    F_upper = UpperBound2(1, Gamma);
    fprintf("Upper bound: %.10f\n", F_upper);

    flb = @(x) arrayfun(@fLowerBound, x);
    int_lower_bound = integral(flb, 0, 1);
    fprintf("Lower bound (integral): %.10f\n", int_lower_bound);
end


function v = FunInner(theta, y, Gamma, delta)
    branch1 = 1 - exp(y) .* (1 - Gamma);
    branch2 = (1 + delta) ./ Gamma .* (1 - exp(theta) .* (1 - Gamma)) .* ...
              (1 - exp(y - theta) .* (1 - Gamma));
    v = min(branch1, branch2);
end


function v = Bound3(theta, Gamma, delta)
    e = exp(1);
    gamma = 1 - (1 / e) - Gamma;
    term1 = (Gamma - theta) ./ (1 - Gamma) .* ...
            (theta - (1 - Gamma) .* (exp(theta) - 1));
    term2_coeff = 1 / (1 - Gamma) .* ...
        (Gamma - (2 * e) / (e - 1) .* (2 ./ delta + 1) .* gamma);
    function v = InnerFunOfY(y)
        v = FunInner(theta, y, Gamma, delta);
    end
    term2_int = integral(@InnerFunOfY, theta, 1);
    v = term1 + term2_coeff .* term2_int;
end


function v = Bound1(theta, Gamma)
    v = (1 - exp(theta) + Gamma * exp(theta));
end


function v = UpperBound2(alph, Gamma)
    v = (alph + 1 - exp(alph - 1) - Gamma);
end
