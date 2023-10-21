function x = NewtonRaphson(f, x0)

    x = x0;
    MinRelErr = 0.0001;
    fx = f(x);
    
    while true
        
        f1x = Derivative(f,x);

        xnew = x - fx/f1x;
        fxnew = f(xnew);
        
        EstRelErr = abs((xnew-x)/x);
        
        x = xnew;
        fx = fxnew;
        
        if EstRelErr<MinRelErr
            
            break;
            
        end 
    end
end