% experimentelle Konvergenzordnungen verschiedener Verfahren
% zur Berechnung der Quadratwurzel einer positiven reellen Zahl x.

x = 5; % Quadratwurzel von x soll berechnet werden
n = 10; % Anzahl der berechneten Folgenglieder


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Heron-Verfahren: x_{n+1} = (x_n + x/x_n) / 2
x0 = x; % Anfangswert x0

xn = zeros(1,n); % Heron-Folge
xn(1) = x0;
for i = 2:n
    xn(i) = (xn(i-1) + x/xn(i-1)) / 2;
end

en = abs(xn - sqrt(x)); % Fehler
alpha = zeros(1,n-2); % Konvergenzordnung
for i = 1:(n-2)
    alpha(i) = log(en(i+1)/en(i+2)) / log(en(i)/en(i+1));
end

plot(1:(n-2),alpha,'ro-','LineWidth',2,'MarkerSize',9)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixpunktverfahren x_{n+1} = phi(x_n) := 1/2 - (x_n)^2/(2*x)  + x_n
% phi'(sqrt(x)) = 1 - 1/sqrt(x) < 1
xn(1) = 2;
for i = 2:n
    xn(i) = 1/2 - xn(i-1)^2/x/2  + xn(i-1);
end

en = abs(xn - sqrt(x)); % Fehler
alpha = zeros(1,n-2); % Konvergenzordnung
for i = 1:(n-2)
    alpha(i) = log(en(i+1)/en(i+2)) / log(en(i)/en(i+1));
end

plot(1:(n-2),alpha,'bx-','LineWidth',2,'MarkerSize',9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Bisektionsverfahren
ab = [0 x]; % Intervall
for i = 1:n
    mid = mean(ab);
    if mid^2 < x^2
        ab(1) = mid; % untere Grenze wird angepasst
    else
        ab(2) = mid; % obere Grenze wird angepasst 
    end
    xn(i) = mean(ab); % Mitte des Intervalls als Approximation der Wurzel
end

en = abs(xn - sqrt(x)); % Fehler
alpha = zeros(1,n-2); % Konvergenzordnung
for i = 1:(n-2)
    alpha(i) = log(en(i+1)/en(i+2)) / log(en(i)/en(i+1));
end

plot(1:(n-2),alpha,'go-','LineWidth',2,'MarkerSize',9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sekantenverfahren
xn(1) = 0;
xn(2) = x;
for i = 3:n
    xn(i) = xn(i-1) - (xn(i-1)-xn(i-2))/(xn(i-1)^2-xn(i-2)^2)*(xn(i-1)^2-x)
end

en = abs(xn - sqrt(x)); % Fehler
alpha = zeros(1,n-2); % Konvergenzordnung
for i = 1:(n-2)
    alpha(i) = log(en(i+1)/en(i+2)) / log(en(i)/en(i+1));
end

plot(1:(n-2),alpha,'cx-','LineWidth',2,'MarkerSize',9)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



plot(1:n,2*ones(1,n),'r:','LineWidth',2) % alpha = 2
plot(1:n,ones(1,n),'b:','LineWidth',2) % alpha = 1
plot(1:n,1/2*ones(1,n),'g:','LineWidth',2) % alpha = 1/2
plot(1:n,(1+sqrt(5))/2*ones(1,n),'c:','LineWidth',2) % alpha = (1+sqrt(5))/2
hold off

title(sprintf(['Numerische Verfahren zur Bestimmung der ' ...
    'Quadratwurzel von x=%d'],x),"FontSize",20)
legend('Heron-Verfahren','$\Phi(z)=z+1/2-\frac{z^2}{2x}$', ...
    'Bisektionsverfahren','Sekantenverfahren','$\alpha=2$', ...
    '$\alpha=1$','$\alpha=1/2$','$\alpha=\frac{1+\sqrt{5}}{2}$', ...
    'Location','northeast','FontSize',12)
xlabel('Anzahl der Iterationen n','FontSize',15)
ylabel('experimentelle Konvergenzordnung $\alpha$','FontSize',15)
ylim([0 3])
xlim([1 n-2])
xticks(1:n)