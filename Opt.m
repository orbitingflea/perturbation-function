function Opt(gamma)
% find the optimal delta in (0, 1).
    global results;
    logdelta = optimizableVariable('logdelta', [-10, 4]);
    loss_fun = @(x)(-Main(gamma, exp(x.logdelta), false));
    results = bayesopt(loss_fun, [logdelta]);

    % disp(results);
    delta = exp(results.XAtMinObjective.logdelta);
    Main(gamma, delta, true);
end
