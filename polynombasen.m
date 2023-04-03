% Plots der Monom-, Lagrange- und Newton-Basis auf dem Intervall [-1,1]

% Definitionsbereich
x = linspace(-1, 1, 1000);

% Ordnung der Polynome
n = 5;

% Basisfunktionen speichern
B = zeros(length(x), n);

% Generiere Monom-Basis
for i = 1:n
    B(:,i) = x.^(i-1);
end

% Generiere Lagrange-Basis
L = zeros(length(x), n);
for i = 1:n
    % Berechne i-tes Lagrange-Polynom
    L_i = 1;
    for j = 1:n
        if j ~= i
            L_i = L_i .* (x - (-1 + 2*(j-1)/(n-1))) ./ (2/(n-1)*(i-1) - 2/(n-1)*(j-1));
        end
    end
    L(:,i) = L_i;
end

% Generiere Newton-Basis
N = zeros(length(x), n);
for i = 1:n
    % Berechne i-tes Newton-Polynom
    N_i = ones(size(x));
    for j = 1:i-1
        N_i = N_i .* (x - (-1 + 2*(j-1)/(n-1)));
    end
    N(:,i) = N_i;
end

% Plot Basis-Funktionen
figure;
% Monom-Basis
ax1 = subplot(1, 3, 1);
for i = 1:n
    plot(x, B(:,i), 'LineWidth', 2);
    hold on;
    grid on;
    title('Monom-Basis');
end
% Lagrange-Basis
ax2 = subplot(1, 3, 2);
for i = 1:n
    plot(x, L(:,i), 'LineWidth', 2);
    hold on;
    grid on;
    title('Lagrange-Basis');
end
% Newton-Basis 
ax3 = subplot(1, 3, 3);
for i = 1:n
    plot(x, N(:,i), 'LineWidth', 2);
    hold on;
    grid on;
    title('Newton-Basis');
    ylim([-1 2])
end
linkaxes([ax1 ax2 ax3],'y');