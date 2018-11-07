#=
\[\vec{F}\equiv-\nabla V(\vec{r})\]

For \(V(\vec{r}) = \frac{Z}{r}\), \(\vec{F}=-\partial_r\frac{Z}{r}\vec{r}=\frac{Z}{r^2}\vec{r}\).

For \(V(\vec{r}) = \frac{Z}{\sqrt{r^2+a^2}}\), \(\vec{F}=\frac{Zr}{(r^2+a^2)^{3/2}}\vec{r}\).
=#

function soft_coulomb(Q::Unitful.Charge, q::Unitful.Charge, a::Unitful.Length)
    Z = Q*q
    (U = r::Unitful.Length -> kₑ*Z/√(r^2 + a^2),
     F = r::Unitful.Length -> (kₑ*Z*r/(r^2+a^2)))
end

export soft_coulomb
