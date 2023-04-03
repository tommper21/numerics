% Dividierte Differenzen zur Berechnung der Koeffizienten des
% Interpolationspolynoms p bzgl. der Newton-Basis und
% anschliessende Auswertung mit dem Horner-Schema


y = [1 2 3 4 5 6 7]; % zu interpolierende Funktionswerte
a = -2; % linke Intervallgrenze
b = 5; % rechte Intervallgrenze
%x = linspace(a,b,n); % äquidistante Stützstellen
x = [-1 0 1 2 2.5 3 5]; % beliebige Stützstellen
t = 2.2; % Stelle, an der Interpolationspolynom p ausgewertet werden soll


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(length(x)==length(y),'x und y haben nicht gleiche Dimension')
for i = 1:length(x)
    assert((a <= x(i)) & (x(i) <= b),'Stützstellen nicht im Intervall')
end
assert((a <= t) & (t <= b),'Auswertungsstelle nicht im Intervall')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dividierte Differenzen
n = length(y);
for m = 1:n-1
    for j = n:-1:m+1
        y(j) = (y(j) - y(j-1)) / (x(j) - x(j-m));
    end
end

%% Horner-Schema
w = y(n);
for k = n-1:-1:1
    w = (t - x(k))*w + y(k);
end
% Ergebnis: w = p(t)

