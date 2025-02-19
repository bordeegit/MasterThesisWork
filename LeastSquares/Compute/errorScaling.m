function Err_scaled = errorScaling(Err, type, printFlag)

if strcmp(type, 'log')
    Err_scaled = sign(Err) .* log1p(abs(Err));
elseif strcmp(type, 'tanh') % Alternatevely we can do hyperbolic tangent scaling
    lambda = 0.1; % Adjust this parameter for different levels of compression
    Err_scaled = tanh(lambda * Err) / tanh(lambda);
else
    error("Scaling type not found, select between log and tanh");
end

if printFlag
    figure;
    subplot(1,2,1);
    scatter(Err(:,1), Err(:,2), 20, linspace(1,10,length(Err)), 'filled');
    title('Original Error');
    
    subplot(1,2,2);
    scatter(Err_scaled(:,1), Err_scaled(:,2), 20, linspace(1,10,length(Err)), 'filled');
    title(strcat(type, ' Scaled Error'));
end 