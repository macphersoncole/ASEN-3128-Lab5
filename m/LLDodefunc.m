function dydt = LLDodefunc(t,state_vec,A)
    
    dydt = A*state_vec;

end