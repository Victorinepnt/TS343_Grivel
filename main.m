clc;
close all;
clear;

%% Parametres 
N = 100;      % length of signal
sigma = 0.1;
u = sqrt(sigma) * (randn(1,N)); % + 1i * randn(1,N));

Nfft = 2^10;
fe = 8 * 10^3;
f_axis = (-1/2 : 1/Nfft : 1/2 - 1/Nfft) * fe;
P = 3; % ordre du processus

%% processus AR - filtrage d'un BBG
b = 1;
a = [1, 0.9, 0.9, 0.9];
x = filter(b,a,u);

% Plus le nombre de pole augmente et plus la richesse frequentielle
% diminue. Le bruit devient colore. Toutes les frequences n'ont pas la meme
% puissance.

% Spectre de puissance d'une realisation du bruit
TF_S = fftshift(fft(x, Nfft));
Spectre_P = abs(TF_S).^2 /Nfft;  % --> observation des resonances, comparer avec th

%% Estimation des paramètres AR et de la variance du processus générateur.

% Yules-Walker
rx = xcorr(x, P, 'unbiased');
rx = rx(length(rx)+1-P:end).';

R = xcorr(x, P-1, 'unbiased');
Rx = zeros(P,P);
for k=1: length(R)
    vect = ones(1,P-abs(k-P)) * R (k);
    diago = diag(vect,P-k);
    Rx = Rx + diago;
end
 
a_est_yw = -inv(Rx) * rx;
a_est_yw = round(abs(a_est_yw),2);

% estimation du bruit
rx = xcorr(x, P, 'unbiased');
rx = rx(1:P).';
sigma_est_yw = a_est_yw.' * rx;


%% LMS normalise

% for i=1:Nit
%     for k=1:N
%         X = x(k:k+P-1);
%         H_k = H(k:k+P-1);
% 
%         h_d(i) = H_k.' * X ;
% 
%         e(i) = d - h_d(i);
%         X = h_d(i) + e(i);
%         H_k = H_k + alpha * (X * e)/norm(X,"fro")^2;
% 
%     end
% end
Nit=5000;
e = zeros(1,Nit);
x_norm  = x/(norm(x)^2);
d = x - u;
H = ones(1,P);
alpha = 0.1;

for k = 1:Nit
    for i = P+1:N

        X = x_norm(i-P:i-1);
        hat_d = H * X.';

        % Calcul de l'erreur
        if k<N
            e(k) = d(k) - hat_d;
        end
%         else
%             e(k) = d()

        % Mise à jour des coefficients AR avec la méthode LMS
        H = H + alpha * e(k) * X;
    end
    % Calcul de la variance du bruit à chaque itération
    %H_mat(k,:) = H;
    %e_lms = d(k) - H * X.';
    sigma_lms_norm(k) = var(e,1);%mean(abs(e).^2);
end

sigma_lms = (norm(x)^2) * sigma_lms_norm;


% evolution des parametres au cours du temps,
% tracer la valeur attendue sur le meme graphe
% observer la fluctuation autour de la vrai valeur (vitesse de convergence et convergence)
%% Figures
% figure,
% plot(f_axis, Spectre_P);
% grid on
% 
% % Affichage des résultats
% disp('Paramètres AR estimés (Yules-Walker) :');
% disp(a_est_yw);
% disp(['Variance du bruit estimée (Yules-Walker) : ', num2str(sigma_est_yw)]);


% Tracez la convergence de la variance du bruit avec la méthode LMS
figure;
plot(1:Nit, sigma_lms, 'o-', 'LineWidth', 1.5);
xlabel('Itération');
ylabel('Variance du bruit estimée');
title('Convergence de la variance du bruit avec la méthode LMS');
grid on;

