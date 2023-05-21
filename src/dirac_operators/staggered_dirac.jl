struct StaggeredDiracOperator <: AbstractDiracOperator
    mass::Float64
    BC::Vector{Int64}
end

function (D::StaggeredDiracOperator)(U::Gaugefield)
    NX, NT, _ = size(U)
    mass = D.mass
    BC_x, BC_t = D.BC

    len = NX * NT
    expiU = cis(U)

    Dstagg = zeros(ComplexF64, len, len)
    return Dstagg
end