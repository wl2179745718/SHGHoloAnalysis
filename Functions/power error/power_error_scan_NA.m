function [] = power_error_scan_NA(fx, FX, FY, dfx, dfy, Re_fz, E_MLB_image_Oli, E_MLB_image, Er)
index_0 = find(fx>0,1);
index_2 = find(fx>1.2,1);
NA = fx(index_0:index_2);

parfor ii = 1:size(NA,2)
filter_out = abs( sqrt(FX.^2+FY.^2) ) < NA(ii);
Er_power = dfx * dfy * sum(sum( filter_out.*Re_fz.*abs(Er).^2 ));
error_power_Oli = dfx * dfy * sum(sum( filter_out.*Re_fz.*abs(Er - E_MLB_image_Oli).^2 ));
error_power_MLB = dfx * dfy * sum(sum( filter_out.*Re_fz.*abs(Er - E_MLB_image).^2 ));
err_pwr_ratio_Oli(ii) = error_power_Oli/Er_power
err_pwr_ratio_MLB(ii) = error_power_MLB/Er_power
end

figure
hold on
plot(NA, err_pwr_ratio_Oli)
plot(NA, err_pwr_ratio_MLB)
xlabel('NA')
legend('Olivier','MLB')
set(gcf, 'Position', get(0, 'Screensize'));
end