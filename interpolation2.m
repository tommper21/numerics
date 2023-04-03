% Erstellung einer Vandermonde-Matrix X mit beliebigen Stützstellen
% x = [x1, ..., xn] bezüglich Monombasis und Plot des Interpolationspoylnom


y = [1 2 3 4 5 6 7]; % zu interpolierende Funktionswerte
a = -2; % linke Intervallgrenze
b = 5; % rechte Intervallgrenze
%x = linspace(a,b,n); % äquidistante Stützstellen
x = [-1 0 1 2 2.5 3 5]; % beliebige Stützstellen


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(length(x)==length(y),'x und y haben nicht gleiche Dimension')
for i = 1:length(x)
    assert((a <= x(i)) & (x(i) <= b),'Stützstellen nicht im Intervall')
end

n = length(y);
X = x' * ones(1,n);
X = X .^ ( ones(n,1) * (0:(n-1)) ); % Vandermonde-Matrix
c = X\y'; % Koeffizientenvektor bzgl. Monom-Basis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolationspolynom mit Interpolationspunkten plotten
xx = linspace(a,b);
plot(xx,polyval(c(end:-1:1),xx),x,y,'o');
legend('Interpolationspolynom','zu interpolierende Funktionswerte $y_i$')
title(['Polynominterpolation mit ',num2str(n), ...
    ' St\"utzstellen auf dem Intervall [', ...
    num2str(a),',',num2str(b),']'])