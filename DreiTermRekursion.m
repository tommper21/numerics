% 3-Term Rekursion für Orthogonalpolynome für \omega(x) = 1 liefert
% Nullstellen vom n+1 Orthogonalpolynom (<f,g> := int(f*g,[a b])).
clear
a = -1; % linke Intervallgrenze
b = 1; % rechte Intervallgrenze
n = 2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3-Term-Rekursion
p = {};
p{1} = [1];
beta(1) = (a+b)/2;
p{2} = [1 -beta(1)];

for k=2:n+1
tmp1 = polyint(conv(conv([1 0],p{k}),p{k})); % x * p_k * p_k
tmp2 = polyint(conv(p{k},p{k})); % p_k * p_k
beta(k) = diff(polyval(tmp1,[a b])) / diff(polyval(tmp2, [a b]));

tmp3 = polyint(conv(p{k-1},p{k-1})); % p_k-1 * p_k-1
% ||p_k|| / ||p_k-1||
gamma(k) = sqrt(diff(polyval(tmp2, [a b])) / diff(polyval(tmp3, [a b])));

p{k+1} = add_poly(conv([1 -beta(k)],p{k}),-gamma(k)^2*p{k-1});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tridiagonalmatrix und Nullstellen
T = zeros(n+1);
% Setzen der Diagonal- und Nebendiagonalelemente
for i = 1:n+1
    T(i,i) = beta(i);
    if i > 1
        T(i,i-1) = -gamma(i);
        T(i-1,i) = -gamma(i);
    end
end

% Eigenwerte von T sind die Nullstellen des n+1-Orthogonalpolynoms, die
% innerhalb des Intervalls [a,b] liegen.
nullstellen = eig(T)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Gewichte der Quadraturformel
tau = ones(1,n+1);
for j = 2:n+1
    tau(j) = tau(j-1) * (-1) / gamma(j);
end
for j=1:n+1
    tmp = 0;
    for k=1:n+1
        tmp = tmp + (tau(k) * polyval(p{k},nullstellen(j)))^2;
    end

    omega(j) = (b-a) / tmp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Output: Nullstellen x_j und Gewichte \omega_j für die Gauss-Quadratur
nullstellen
omega

%% Anwendung der Gauss-Quadratur-Formel
f = @(x) exp(x)./(x.^2+1); % Beispielfunktion, die integriert werden soll
res = sum(f(nullstellen).*omega); % Approximation des Integrals int(f,a,b)
%mittels Gauss-Quadratur


