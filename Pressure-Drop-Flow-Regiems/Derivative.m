function df = Derivative(f, x)
    
    dx = 1e-6;
    
    df = (f(x+dx)-f(x-dx))/(2*dx);

end