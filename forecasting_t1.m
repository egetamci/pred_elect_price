function forec_s = forecasting_t1(p_rf, periodf_arr, so_euro_appf, N22, N11, YY)
st11 = N11 - N22;
for i = 1:p_rf
for j = 1:N11
forec_s(i + (j - 1) * p_rf, 1) = YY(i + (j - 1) * p_rf, 1) + periodf_arr(i, 1);
end
end
forec_s(1:p_rf * N22, 1) = forec_s(1:p_rf * N22, 1) + so_euro_appf;
so_euro_app_ort = zeros(p_rf, 1);
for i = 1:N22
so_euro_app_ort = so_euro_app_ort + so_euro_appf(1 + (p_rf * (i -1)):p_rf * i, 1);
end
so_euro_app_ort22 = so_euro_app_ort ./ N22;
for kk = 1:st11
so_euro_appf1(1 + (kk - 1) * p_rf:kk * p_rf, 1) = so_euro_app_ort22;
end
forec_s(p_rf * N22 + 1:p_rf * N11, 1) = forec_s(p_rf * N22 + 1:p_rf *N11, 1) + so_euro_appf1;
end