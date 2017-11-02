function params_scaled = scalePeak(params, mp, ages)
% params_scaled = scalePeak(params, mp, ages)

params_scaled = params;
for ss = 1:size(params,1)
    % Calculate how much higher the peak is than the mean over this age
    % range, and then scale gaussian peak parameter appropriately
    params_scaled(ss,3) =  params(ss,3).*mp/mean(evalgaussian1d(params(ss,:),ages));
    % Test to make sure the new mean is equal to what it should be within
    % numerical precision
    assert(mean(evalgaussian1d(params_scaled(ss,:),ages)) - mp < 10^-16);
end