function table = nonDimensionalize(SM,m,g,V,I_y)

    table = [SM(1,1)/m SM(2,1)/m 0 -g;...
        SM(1,2)/(m-SM(5,2)) SM(2,2)/(m-SM(5,2)) (SM(3,2)+m*V)/(m-SM(5,2)) 0;...
        (SM(1,3)+(SM(5,3)*SM(1,2)/(m-SM(5,2))))/I_y (SM(2,3)+(SM(5,3)*SM(2,2)/(m-SM(5,2))))/I_y (SM(3,3)+(SM(5,3)*(SM(3,2)+m*V)/(m-SM(5,2))))/I_y 0;...
        0 0 1 0];

end