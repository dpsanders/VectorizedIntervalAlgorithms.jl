
function V_ij(r²)
    r⁻⁶ = (r²)^(-3)

    4*((r⁻⁶ - 0.5)^2) - 1.0  # - 0.25
end


function V3(x2, x3, y3)

    r12² = (x2)^2
    r13² = (x3)^2 + (y3)^2
    r23² = (x3 - x2)^2 + (y3)^2   # z3 = z2 = 0

    return V_ij(r12²) + V_ij(r13²) + V_ij(r23²)
end

V3(X) = V3(X[1], X[2], X[3])


function V4(x2, x3, x4, y3, y4, z4)

    r14² = (x4)^2 + (y4)^2 + (z4)^2
    r24² = (x4 - x2)^2 + (y4)^2 + (z4)^2
    r34² = (x4 - x3)^2 + (y4 - y3)^2 + (z4)^2

    return V3(x2, x3, y3) + V_ij(r14²) + V_ij(r24²) + V_ij(r34²)

end

V4(X) = V4(X[1], X[2], X[3], X[4], X[5], X[6])


function V5(x2, x3, x4, x5, y3, y4, y5, z4, z5)

    V4_value = V4(x2, x3, x4, y3, y4, z4)

    r15² = (x5)^2 + (y5)^2 + (z5)^2
    r25² = (x5 - x2)^2 + (y5)^2 + (z5)^2
    r35² = (x5 - x3)^2 + (y5 - y3)^2 + (z5)^2
    r45² = (x5 - x4)^2 + (y5 - y4)^2 + (z5 - z4)^2

    return V4_value + V_ij(r15²) + V_ij(r25²) + V_ij(r35²) + V_ij(r45²)
end

V5(X) = V5(X[1], X[2], X[3], X[4], X[5], X[6], X[7], X[8], X[9])


function V6(x2, x3, x4, x5, x6, y3, y4, y5, y6, z4, z5, z6)

    V5_value = V5(x2, x3, x4, x5, y3, y4, y5, z4, z5)

    r16² = (x6)^2 + (y6)^2 + (z6)^2
    r26² = (x6 - x2)^2 + (y6)^2 + (z6)^2
    r36² = (x6 - x3)^2 + (y6 - y3)^2 + (z5)^2
    r46² = (x6 - x4)^2 + (y6 - y4)^2 + (z5 - z4)^2
    r56² = (x6 - x5)^2 + (y6 - y5)^2 + (z6 - z5)^2

    return V5_value + V_ij(r16²) + V_ij(r26²) + V_ij(r36²) + V_ij(r46²) + V_ij(r56²)
end

V6(X) = V6(X[1], X[2], X[3], X[4], X[5], X[6], X[7], X[8], X[9], X[10], X[11], X[12])
