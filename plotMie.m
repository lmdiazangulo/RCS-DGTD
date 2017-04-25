thP = linspace(0,pi,250); % Stores theta components.
phP = pi;%linspace(-pi,pi,150);% Stores phi components.

for p=1:numel(thP)
    for q=1:numel(phP)
        scaledFreq = intFreq * sphereRadius;
        [es_theta, es_phi] = ...
         mie_pec(1, scaledFreq, thP(p), phP(q), 50);
        RCSMie(p,q) = abs(es_theta^2 + es_phi^2) *  sphereRadius^2;
    end
end

plot(thP, 10*log10(RCSMie(:,1)),'.-r');
grid on;