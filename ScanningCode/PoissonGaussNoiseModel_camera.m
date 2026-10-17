function noisy_out = PoissonGaussNoiseModel_camera(clean_in, noise_params)

% noise model for camera

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Parsing

% noise_params = check_noise_params(noise_params);                           % Check the noise parameter struct for missing elements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw noisy outputs
% 修复1：替换 poissrnd。利用正态分布近似泊松分布，应对新版MATLAB
noisy_out = noise_params.gamma * max(0, round(clean_in + sqrt(clean_in) .* randn(size(clean_in))));

noisy_out = noisy_out  + noise_params.offset;

% 修复2：替换 normrnd。用基础 randn 生成相同尺寸的正态分布噪声
noisy_out = noisy_out + noise_params.sigma_p * randn(size(noisy_out, 1), size(noisy_out, 2));

end