close all;
%% Visualisation
figure;
imagesc(scanning_phi, scanning_theta, max(abs(beams), [], 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Пространство угол места-азимут (оторбражение max дальность)');
grid on;

figure;
mesh(scanning_phi, scanning_theta, max(abs(beams), [], 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Пространство угол места-азимут (отображение rms дальность)');
grid on;

figure;
mesh(scanning_phi, scanning_theta, rms(abs(beams), 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Пространство угол места-азимут (отображение rms дальность)');
grid on;


figure;
mesh(scanning_phi, scanning_theta, max(abs(filtered_beams), [], 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Отфильтрованные угловые координаты цели (отображение max дальность)');
grid on;

figure;
mesh(scanning_phi, scanning_theta, rms(abs(filtered_beams), 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Отфильтрованные угловые координаты цели (отображение rms дальность)');
grid on;

figure;
plot(scanning_phi, max(abs(filtered_beams), [], 3)')
title('Семейство графиков зависимости отлика от азимута при разных углах места')
xlabel('Азимут, градусы')
ylabel('Амплитуда')
xlim("tight")

figure;
hold on; grid on;
plot(scanning_phi, max(abs(filtered_beams(elev_idx, :, :)), [], 3), 'DisplayName', 'max')
plot(scanning_phi, rms(filtered_beams(elev_idx, :, :), 3), 'DisplayName', 'rms')
xlabel('Азимут (градусы)')
xlim("tight")
ylabel('Амплитуда')
title('Зависимость амплитуды от азимута (отображение среза угол места 15 градусов)')
legend('Location', 'northwest');

figure;
mesh(1 : size(filtered_beams, 3), ...
    scanning_phi, ...
    reshape(max(abs(beams), [], 1), size(beams, 2), size(beams, 3)));
xlabel('Дальность (отсчеты)');
ylabel('Азимут (градусы)');
title('Пространство азимут-дальность(отображение max угол места)');
grid on;

figure;
azimuth_range_filtered = reshape(max(abs(filtered_beams), [], 1), size(filtered_beams, 2), size(filtered_beams, 3));
mesh(1 : size(filtered_beams, 3), ...
    scanning_phi, ...
    azimuth_range_filtered);
xlabel('Дальность (отсчеты)');
ylabel('Азимут (градусы)');
title('Отфильтрованное пространство азимут-дальность(отображение max угол места)');
grid on;

figure;
imagesc(r_idx :  750, ...
    scanning_phi, ...
    azimuth_range_filtered(:, r_idx : 750));
xlabel('Дальность (отсчеты)');
ylabel('Азимут (градусы)');
title('Отфильтрованное пространство азимут-дальность (отображение max угол места)');
grid on;


