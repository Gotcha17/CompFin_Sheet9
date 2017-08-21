function [X_exact, X_Euler, X_Milshtein] = Sim_Paths_GeoBM(X0, mu, sigma, T, N)
    delta_t = T/N;
    // Simulating N std. nor. values
    delta_W = grand(N, 1, "nor", 0, sqrt(delta_t));
    // Setting starting parameters
    a = mu - sigma^2/2;
    b = sigma;
    X_exact(1) = X0;
    X_Euler(1) = X0;
    X_Milshtein(1) = X0; 
    for i=2:N+1
        Y = sum(a*delta_t + b*delta_W(1:i-1));
        // Exact simulation
        X_exact(i) = X0*exp(Y);
        // Euler simulation
        X_Euler(i) = X_Euler(i-1)+X_Euler(i-1)*mu*delta_t+...
                     X_Euler(i-1)*sigma*delta_W(i-1);
        // Milshstein simulation
        X_Milshtein(i) = X_Milshtein(i-1)+X_Milshtein(i-1)*mu*delta_t+...
                         X_Milshtein(i-1)*sigma*delta_W(i-1)+...
                         X_Milshtein(i-1)*0.5*b^2*(delta_W(i-1)^2-delta_t);
    end
endfunction
X0=100; mu=0.05; sigma=0.2; T=2;
N = [10, 100, 1000, 10000];
// Prepare plotting setups
clf
scf(0)
a = gca()
a.tight_limits = 'on'
a.auto_scale = 'on'

// Call function with test data and Plot results in the matrix view.
// We are walking throught the for loop, calling function for each N.
for i = 1 : length(N)
    
    [X_exact, X_Euler, X_Milshtein] = Sim_Paths_GeoBM(X0, mu, sigma, T, N(i))
    subplot(2,2,i)
    plot2d([0 : N(i)],[X_exact, X_Euler, X_Milshtein],[1, 2, 3])
    xtitle('Simulated paths for ' + string(N(i)) + ' steps', 'N', 'X')
    legends(['Exact';'Euler';'Milshtein'],[1, 2, 3], opt="lr")
    
end



