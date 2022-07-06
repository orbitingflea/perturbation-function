% This file is to find the optimal delta for a specified Gamma.
% This is not needed if you only want to verify the result.

function Opt(gamma)
    % Find the optimal delta in (0, 1) for gamma.
    global results;
    logdelta = optimizableVariable('logdelta', [-10, 4]);
    loss_fun = @(x)(-Main(gamma, exp(x.logdelta), false));
    results = bayesopt(loss_fun, [logdelta]);

    delta = exp(results.XAtMinObjective.logdelta);
    Main(gamma, delta, true);
end
