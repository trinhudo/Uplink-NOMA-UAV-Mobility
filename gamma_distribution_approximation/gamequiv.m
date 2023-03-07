function [w, t] = gamequiv(m, s, tol)
    M = 2;
    warning off;
    
    epsilon = m - floor(m);
    y_true = (1+s/m).^(-epsilon);
    
    syms x;
    while(1)
        weight= sym('w',[1, M]); weight(M) = 1 - sum(weight(1:M-1));
        scale = sym('t',[1, M]);
        
        y_pred = 0;
        for k = 1:M
            y_pred = y_pred + weight(k) ./ (1+scale(k)*x);
        end
        
        lower_bounds = [-Inf*ones(1,M) zeros(1,M)];
        
        options = fitoptions('Method', 'NonlinearLeastSquares',...
                    'Lower', lower_bounds);
        [f, gof, ~] = fit(s, y_true, char(y_pred), options);
        sprintf('RMSE: %f.3', gof.rmse)
        
        if gof.rmse < tol
            break;
        end
        M = M + 1;
    end
    
    coeff = coeffvalues(f);
    t = coeff(1:M);
    w = coeff(M+1:2*M-1); w = [w 1-sum(w)];
    f
end