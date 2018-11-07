# 𝐅 = q𝐄 (+ q𝐯 × 𝐁)
function lorentz(E::ElectricFields.LinearField, q::Unitful.Charge=-1.60217662e-19u"C", m::Unitful.Mass=1u"me")
    (du,u,p,t) -> begin
        du[1] = u[2]
        du[2] = q*E(t)/m
    end
end

function lorentz(E::ElectricFields.LinearField, U::Function, ndims::Integer,
                 q::Unitful.Charge=-ec, m::Unitful.Mass=1u"me")
    Fᵤ = U(q).F
    ξ = vcat(1,zeros(ndims-1))
    (du,u,p,t) -> begin
        𝐫 = u.x[1]
        r = norm(𝐫)
        𝐫̂ = 𝐫 * (ustrip(r) == 0 ? 0/u"m" : 1/r)
        du.x[1] .= u.x[2]
        du.x[2] .= (q*E(t)*ξ + Fᵤ(r)*𝐫̂./u"m")/m
    end
end

function integrate(E::ElectricFields.LinearField,
                   interval::NTuple{2, Unitful.Time},
                   x₀::Unitful.Length, v₀::Velocity,
                   args...; solver=Tsit5, kwargs...)
    prob = ODEProblem(lorentz(E), ArrayPartition([x₀],[v₀]), interval)
    solve(prob, solver(), args...; kwargs...)
end

function integrate(E::ElectricFields.LinearField, U::Function,
                   interval::NTuple{2, Unitful.Time},
                   𝐫₀::Vector{<:Unitful.Length}, 𝐯₀::Vector{<:Velocity},
                   args...; solver=Tsit5, kwargs...) 
    prob = ODEProblem(lorentz(E,U,length(𝐫₀)), ArrayPartition(𝐫₀,𝐯₀), interval)
    solve(prob, solver(), args...; kwargs...)
end

integrate(E::ElectricFields.LinearField, U::Function,
          interval::NTuple{2, Unitful.Time},
          𝐫₀::Unitful.Length, 𝐯₀::Velocity,
          args...; kwargs...) =
              integrate(E, U, interval, [𝐫₀], [𝐯₀], args...; kwargs...)

export integrate
