% Romberg-Verfahren zur numerischen Berechnung des Integrals einer Funktion
% f über das Intervall [a,b].
% Extrapolation der summierten Trapezregel mit alpha = 2.
% Auswertung mit Hilfe des Neville-Verfahrens.

a = 1;
b = 2;
f = @(x) exp(x)./x + sin(x*pi);
n = 8; % Anzahl der berechneten Folgegliedern
k = 10.^(1:n); % Anzahl der Stützstellen für die Trapezregel
h = (b-a) ./ k; % Intervallbreite
phi = zeros(1,n);
res = integral(f,a,b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Berechnung Phi(h_k) mithilfe der Trapezregel
for i = 1:n
    phi(i) = 0.5 * h(i) * (f(a) + 2*sum(f(a+h(i):h(i):b-h(i))) + f(b));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Neville-Verfahren
x = h.^2; % Stützstellen für Interpolationspolynom p_n
y = phi; % Funktionswerte für Interpolationspolynom p_n
for m = 1:n-1
    for j = 1:n-m
        y(j) = ((0 - x(j))*y(j+1) - (0 - x(j+m))*y(j)) / (x(j+m)-x(j));
    end
end
% Output: Approximation y(1) = p_n(0)
fprintf(['\nApproximation für Integral: %f\n' ...
    'Relativer Fehler mit Neville: %f\n' ...
    'Relativer Fehler ohne Neville: %f\n'] ...
    ,y(1),abs(1-y(1)/res),abs(1-phi(n)/res));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot

loglog(h,abs(phi-res),'o-');
hold on
plot(h,h.^2);
hold off

title('Romberg-Verfahren (Extrapolation der summierten Trapezregel)', ...
    'FontSize',20);
legend('Fehler $|\Phi (h)-\int_a^b f(x)dx|$', ...
    'Konvergenzrate $\alpha = 2$','FontSize',12,'Location','northwest')
xlabel('Intervallbreite $h$','FontSize',15);
ylabel('Absoluter Fehler','FontSize',15);