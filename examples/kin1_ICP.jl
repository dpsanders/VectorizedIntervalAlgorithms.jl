using Symbolics



function f( (x1, x2, x3, x4, x5, x6 ) )

    SA[-0.4077 + cos(x2)*cos(x6) + cos(x3)*cos(x6) + cos(x4)*cos(x6) + cos(x5)*sin(x2)*sin(x6) - cos(x5)*sin(x3)*sin(x6) - cos(x5)*sin(x4)*sin(x6),
           -1.9115 + cos(x5)*sin(x1) + cos(x1)*cos(x2)*sin(x5) + cos(x1)*cos(x3)*sin(x5) + cos(x1)*cos(x4)*sin(x5),
           -1.9791 + sin(x2)*sin(x5) + sin(x3)*sin(x5) + sin(x4)*sin(x5),
           -4.0616 + 3*cos(x1)*cos(x2) + 2*cos(x1)*cos(x3) + cos(x1)*cos(x4),
           -1.7172 + 3*cos(x2)*sin(x1) + 2*cos(x3)*sin(x1) + cos(x4)*sin(x1),
           -3.9701 + 3*sin(x2) + 2*sin(x3) + sin(x4)
    ]
end



Cs = Contractor.(Ref(vars), exs)

@variables x2, x3, x4

Contractor(-3.9701 + 3*sin(x2) + 2*sin(x3) + sin(x4))


