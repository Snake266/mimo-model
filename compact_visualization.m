close all;
%% Visualisation
figure;
subplot(2, 1, 1)
imagesc(scanning_phi, scanning_theta, max(abs(beams), [], 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Пространство угол места-азимут (оторбражение max дальность)');
grid on;
subplot(2, 1, 2)
mesh(scanning_phi, scanning_theta, rms(abs(beams), 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Пространство угол места-азимут (отображение rms дальность)');
grid on;



figure;
subplot(2, 1, 1)
mesh(scanning_phi, scanning_theta, max(abs(filtered_beams), [], 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Отфильтрованные угловые координаты цели (отображение max дальность)');
grid on;

subplot(2, 1, 2)
mesh(scanning_phi, scanning_theta, rms(abs(filtered_beams), 3));
xlabel('Азимут (градусы)');
ylabel('Угол места (градусы)');
title('Отфильтрованные угловые координаты цели (отображение rms дальность)');
grid on;

figure;
subplot(2, 1, 1)
plot(scanning_phi, max(abs(filtered_beams), [], 3)')
title('Семейство графиков зависимости отлика от азимута при разных углах места')
xlabel('Азимут, градусы')
ylabel('Амплитуда')
xlim("tight")

subplot(2, 1, 2)
hold on; grid on;
plot(scanning_phi, max(abs(filtered_beams(elev_idx, :, :)), [], 3), 'DisplayName', 'max')
plot(scanning_phi, rms(filtered_beams(elev_idx, :, :), 3), 'DisplayName', 'rms')
xlabel('Азимут (градусы)')
xlim("tight")
ylabel('Амплитуда')
title('Зависимость амплитуды от азимута (отображение среза угол места 15 градусов)')
legend('Location', 'northwest');

figure;
azimuth_range_filtered = reshape(max(abs(filtered_beams), [], 1), size(filtered_beams, 2), size(filtered_beams, 3));
subplot(2, 1, 1)
mesh(1 : size(filtered_beams, 3), ...
    scanning_phi, ...
    azimuth_range_filtered);
xlabel('Дальность (отсчеты)');
ylabel('Азимут (градусы)');
title('Отфильтрованное пространство азимут-дальность (отображение max угол места)');
grid on;

subplot(2, 1, 2)
imagesc(r_idx :  750, ...
    scanning_phi, ...
    azimuth_range_filtered(:, r_idx : 750));
xlabel('Дальность (отсчеты)');
ylabel('Азимут (градусы)');
title('Отфильтрованное пространство азимут-дальность (отображение max угол места)');
grid on;


