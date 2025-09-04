function result = current_and_lineloss(r1, theta1_deg, r2, theta2_deg, real_part3, imaginary_part3)
    % Convert phase angles from degrees to radians
    theta1_rad = deg2rad(theta1_deg);
    theta2_rad = deg2rad(theta2_deg);

    % Compute real and imaginary parts of the first complex number
    z1_real = r1 * cos(theta1_rad);
    z1_imaginary = r1 * sin(theta1_rad);

    % Compute real and imaginary parts of the second complex number
    z2_real = r2 * cos(theta2_rad);
    z2_imaginary = r2 * sin(theta2_rad);

    % Subtract the real and imaginary parts of the second complex number from the first one
    sub_real = z1_real - z2_real;
    sub_imaginary = z1_imaginary - z2_imaginary;

    % Create the third complex number in rectangular form
    z3 = real_part3 + 1i * imaginary_part3;

    % Divide the result by the third complex number
    result_rectangular = (sub_real + 1i * sub_imaginary) / z3;

    % Convert the result from rectangular to polar form
    result_magnitude = abs(result_rectangular);
    result_angle_rad = angle(result_rectangular);
    result_angle_deg = rad2deg(result_angle_rad);

    % Construct the result in polar form
    result = [result_magnitude, result_angle_deg];

    %real loss
    real_loss= result_magnitude.^2*real_part3*100;
    reactive_loss= result_magnitude.^2* imaginary_part3*100;
    result=[result real_loss reactive_loss];
end
