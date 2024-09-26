f = @(y,x) -(1/6 + (pi.*sin(pi.*x))/(1.6 - cos(pi.*x))).*y;




function ret = do_euler(f, start_x, start_y, step_size, end_x, tolerance)

    prev_end_y = 0;

    for i = 1:100
        x = start_x:step_size:end_x;
        y = zeros(size(x));
        y(1) = start_y;

        disp(["Number of steps: ", length(x)]);

        for n = 1:length(y)
            y(n+1) = y(n) + step_size * f(y(n), x(n));
        end

        ret = y;
        
        
        error = abs(prev_end_y - y(end));

        %Exit condition
        if error < tolerance
            if(i ~= 1)
                break;
            end
        end
        
        prev_end_y = y(end);
        step_size = step_size * 0.5; %0.5 = 1/2. Datorn räknar multiplikation snabbare än division
    end
end



do_euler(f, 0, 0, 0.5, 6, 10^(-4))