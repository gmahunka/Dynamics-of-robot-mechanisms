clc; clear; close all;
holvan=0;
% Define parameter space
p_values = 0:0.005:0.16;
d_values = 0:0.005:0.9;
[P, D] = meshgrid(p_values, d_values);  % Create grid
stability_map = zeros(size(P));  % Store stability (1 = stable, 0 = unstable)

% Symbolic variables
syms p d mu eta
A = [mu, -1,  0,  0;
      0, mu-1, -1, -1/2;
      0, 0, mu-1, -1;
     -d, (p+d), 0, mu];

% Solve characteristic equation for eigenvalues
eq = det(A) == 0;
eigenvalues = solve(eq, mu);  % Solve for mu

% Iterate over grid points
for i = 1:length(p_values)
    for j = 1:length(d_values)
        p1 = p_values(i);
        d1 = d_values(j);

        % Substitute (p,d) values into eigenvalues
        lambda_vals = double(subs(eigenvalues, [p, d], [p1, d1]));

         if any(abs(lambda_vals) >= 1)
            stability_map(j, i) = 0; % Unstable
            continue;
        end
        

        eq_eta = lambda_vals == (eta + 1)/(eta - 1);
        ctr=0;
        for m=1:length(lambda_vals)
            eta_sol = solve(eq_eta(m), eta);
            eta_real_part = double(real(eta_sol));
            if eta_real_part < 0
                ctr=ctr+1;
            end
        end
       
        if ctr==length(lambda_vals)
            stability_map(j, i) = 1; % Stable
        else
            stability_map(j, i) = 0; % Unstable
        end
        holvan=holvan+1
    end
end

% Plot stability map
figure;
imagesc(p_values, d_values, stability_map);
set(gca,'YDir','normal'); % Fix the axis orientation
colormap([1 0 0; 0 0 1]); % Red for unstable, Blue for stable
xlabel('p [-]');
ylabel('d [-]');
title('Stability Map (Blue: Stable, Red: Unstable)');
hold on;
h3 =@(p,d) - 16*d^3 - 8*d^2*p - 64*d^2 + 64*d - 16*p^2 - 96*p;
fimplicit(h3, [0, 1.5, 0, 1.1], 'LineWidth', 3, 'Color', 'k');
plot([0, 0], [0,0.8284 ], 'LineWidth', 3, 'Color', 'k'); 

