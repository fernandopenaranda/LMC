"""
McCann-Koshino n layer model
k- adimensional momentum (module)
ξ is the valley
ϕ -> k_vect = k * (cos(θ), sin(θ))
θ in radians.

The model is written in meV, fs, and K. 
The module Optics_in_the_length_gauge is written in eV, s, K, so there is a unit convention fixer
"""
const s0 = SA[1 0; 0 1]
const σ0 = SA[1 0; 0 1]
const σx = SA[0 1; 1 0]
const σy = SA[0 -1im; 1im 0]
const σz = SA[1 0; 0 -1]
const a0 = 2.46 # lattice constant in Å
const d = 3.35


@with_kw struct Params_rhombohedral
    ξ::Int #+1 or -1
    μ::Float64
    γ0::Float64
    γ1::Float64
    γ2::Float64
    γ3::Float64
    γ4::Float64
    Delta_Ez::Float64
    Valley_asym::Float64
end
params_rhombohedral() = Params_rhombohedral(1, 0, 3160, 390, -20, 315, 44, 0, 0)
#____________________________________________________________________________________________
#
#             N - Layer model 
#
#____________________________________________________________________________________________
"""
returns all the n1,n2,n3 integers that satisfy N = n1 + 2n2 + 3n3 with N the number of layers
"""

abc_Nlayer(N, k::Array, p = params_rhombohedral()) = abc_Nlayer(N, norm(k), atan(k[2], k[1]), p)

function abc_Nlayer(N, k::Number, θ, p = params_rhombohedral()) 
    Ezterm = p.Delta_Ez *[1 0;0 -1] # oppossite signs top and bottom layer. Preserves PHS
    Valleypolarization = p.Valley_asym * [1 0; 0 1] * p.ξ/2
    triads = find_triads(N)
    [-p.μ X_NLG(triads, k, θ, p); conj(X_NLG(triads, k, θ, p)) -p.μ] +  y(p)*(k)^2 .* 1I + Ezterm + Valleypolarization 
end
y(p) = 3/2*(p.γ0*p.γ4/p.γ1)

function X_NLG(triads, k, θ, p)
    θ += ifelse(p.ξ==1, 0, π) # to account for the two valleys
    val = 0+0im
    for (n1, n2, n3) in triads
        val += fact_triad(n1,n2,n3) * 1/(-p.γ1)^(n1+n2+n3-1) *
            (√3/2 * p.γ0*k*cis(θ))^n1 * (√3/2 * p.γ3*k*cis(-θ))^n2 * (p.γ2/2)^n3
    end
    return val
end

function find_triads(N::Int)
    triads = []
    for n3 in 0:div(N, 3)
        for n2 in 0:div(N - 3n3, 2)
            n1 = N - 2n2 - 3n3
            if n1 >= 0
                push!(triads, (n1, n2, n3))
            end
        end
    end
    return triads
end

#____________________________________________________________________________________________
#       N derivative
#____________________________________________________________________________________________
dhxNlg(N, k::Array, p) = a0 * dhxNlg(N, norm(k), atan(k[2], k[1]), p)
dhyNlg(N, k::Array, p) = a0 * dhyNlg(N, norm(k), atan(k[2], k[1]), p)

function dhxNlg(N, k::Number, θ, p)
    θ += ifelse(p.ξ==1, 0, π)
    triads = find_triads(N)
    m11 = 2*k*y(p)*cos(θ)
    m22 = m11
    m12 = 0im
    for (n1,n2,n3) in triads
        m12 += mN12x(n1,n2,n3, k, θ, p)
    end
    return [m11 m12; conj(m12) m22]
end

function dhyNlg(N, k::Number, θ, p) 
    θ += ifelse(p.ξ==1, 0, π)
    triads = find_triads(N)
    m11 = 2*k*y(p)*sin(θ)
    m22 = m11
    m12 = 0im
    for (n1,n2,n3) in triads
        m12 += mN12y(n1,n2,n3, k, θ, p)
    end
    return [m11 m12; conj(m12) m22]
end

# lmc and qah are ok.
# mN12y(n1,n2,n3, k, θ, p) = dkcommon(n1,n2,n3, k, θ, p) * (cos(θ)*(n1+n2) - sin(θ)*p.ξ*1im*(n1-n2))
# mN12x(n1,n2,n3, k, θ, p) =  dkcommon(n1,n2,n3, k, θ, p) * (sin(θ)*(n1+n2) +cos(θ)*p.ξ*1im*(n1-n2))
# dkcommon(n1,n2,n3, k, θ, p) = prefactor(n1,n2,n3, p) * k^(n1+n2-1) * cis(p.ξ * θ*(n1-n2))

# lmc and qah are ok. Wierd way to introduce the xi dependency but robust
mN12y(n1,n2,n3, k, θ, p) =  dkcommon(n1,n2,n3, k, θ, p)  * (cos(θ-ifelse(p.ξ==1, 0, π)) *  (n1+n2) * cis(θ*(n1-n2)) - 
    sin(θ-ifelse(p.ξ==1, 0, π)) * 1im *(n1-n2) * cis(θ*(n1-n2))) #parece que los ejes est'an girados por eso x->y
mN12x(n1,n2,n3, k, θ, p) =  dkcommon(n1,n2,n3, k, θ, p) * p.ξ * (sin(θ-ifelse(p.ξ==1, 0, π)) * (n1+n2) * cis(θ*(n1-n2)) + 
   cos(θ-ifelse(p.ξ==1, 0, π)) * 1im *(n1-n2) * cis(θ*(n1-n2)))
dkcommon(n1,n2,n3, k, θ, p) = (prefactor(n1,n2,n3, p) * k^(n1+n2-1)) 


# both ok rest is ook [1]
# mN12y(n1,n2,n3, k, θ, p) =  (cos(θ-ifelse(p.ξ==1, 0, π)) * dkXNLG(n1,n2,n3, k, θ, p) - 
#     sin(θ-ifelse(p.ξ==1, 0, π)) * dθXNLG(n1,n2,n3, k, θ, p)) #parece que los ejes est'an girados por eso x->y
# mN12x(n1,n2,n3, k, θ, p) =  p.ξ*(sin(θ-ifelse(p.ξ==1, 0, π)) * dkXNLG(n1,n2,n3, k, θ, p) + 
#     cos(θ-ifelse(p.ξ==1, 0, π)) * dθXNLG(n1,n2,n3, k, θ, p))
# dkXNLG(n1,n2,n3, k, θ, p) = prefactor(n1,n2,n3, p) * k^(n1+n2-1)* (n1+n2) * cis(θ*(n1-n2))
# dθXNLG(n1,n2,n3, k, θ, p) = prefactor(n1,n2,n3, p) * k^(n1+n2-1) * 1im *(n1-n2) * cis(θ*(n1-n2))

prefactor(n1,n2,n3, p) = fact_triad(n1,n2,n3)*(1/(-p.γ1)^(n1+n2+n3-1))*(p.γ2/2)^n3 * (√3/2*p.γ0)^n1 * (√3/2*p.γ3)^n2 

fact_triad(n1,n2,n3) = 
    factorial(n1+n2+n3)/(factorial(n1)*factorial(n2)*factorial(n3))

#____________________________________________________________________________________________
# N second derivative vxx
#____________________________________________________________________________________________
dhxxNlg(N, k::Array, p) = a0 * dhxxNlg(N, norm(k), atan(k[2], k[1]), p)
function dhxxNlg(N, k::Number, θ, p)
    θ += ifelse(p.ξ==1, 0, π)
    triads = find_triads(N)
    m11 = 2*y(p)
    m22 = m11
    m12 = 0im
    for (n1,n2,n3) in triads
        m12 += mN12xx(n1,n2,n3, k, θ, p)
    end
    return [m11 m12; conj(m12) m22]
end

# new: works too
mN12xx(n1,n2,n3, k, θ, p) =  sin(θ-ifelse(p.ξ==1, 0, π)) * mN12xxk(n1,n2,n3, k, θ, p) + cos(θ-ifelse(p.ξ==1, 0, π))/k * mN12xxθ(n1,n2,n3, k, θ, p)
mN12xxk(n1,n2,n3, k, θ, p) = prefactor(n1,n2,n3,p) * cis(θ * (n1-n2)) *
     (sin(θ-ifelse(p.ξ==1, 0, π)) * (n1+n2) + cos(θ-ifelse(p.ξ==1, 0, π)) * 1im * (n1-n2)) * (n1+n2-1) * k^(n1+n2-2)
mN12xxθ(n1,n2,n3, k, θ, p) = prefactor(n1,n2,n3,p) * k^(n1+n2-1) * cis(θ * (n1-n2)) * (
     (n1+n2) * (1im *(n1-n2)*sin(θ-ifelse(p.ξ==1, 0, π)) + cos(θ-ifelse(p.ξ==1, 0, π))) + 1im *(n1-n2)* (1im *(n1-n2)*cos(θ-ifelse(p.ξ==1, 0, π)) -sin(θ-ifelse(p.ξ==1, 0, π)))
)


# old: works ok [1]
# mN12xx(n1,n2,n3, k, θ, p) =  prefactor(n1,n2,n3,p) * k^(n1+n2-2) * cis(θ*(n1-n2)) * (
#     1im * (n1 - n2)^2*cos(θ-ifelse(p.ξ==1, 0, π))^2 + (n1 - n2) * (-3 + n1 + n2) * cos(θ-ifelse(p.ξ==1, 0, π))*sin(θ-ifelse(p.ξ==1, 0, π)) +
#     1/2 * (n1+n2) *(n1+n2-1 - (-3 + n1 + n2) *cos(2*(θ-ifelse(p.ξ==1, 0, π))) +  1im*(n1-n2) * sin(2*(θ-ifelse(p.ξ==1, 0, π)))))


rzNlg(N, ψs) = d * (N-1) * ψs' * σz * ψs # in Angstroms
#############
# Pentalayer
#############
abc_pentalayer(k::Array, p = params_rhombohedral()) = abc_Nlayer(5, norm(k), atan(k[2], k[1]), p)

# abc_pentalayer(k::Array, p = params_rhombohedral()) = abc_pentalayer(norm(k), atan(k[2], k[1]), p)

# function abc_pentalayer(k::Number, θ, p = params_rhombohedral()) 
#     Ezterm = p.Delta_Ez *[1 0;0 -1] # oppossite signs top and bottom layer. Preserves PHS
#     Valleypolarization = p.Valley_asym * [1 0; 0 1] * p.ξ/2
#     [-p.μ X_5LG(k, θ, p); conj(X_5LG(k, θ, p)) -p.μ] +  y(p)*(k)^2 .* 1I + Ezterm + Valleypolarization + 0*[0 1+1im;1-1im 0]
# end

# function X_5LG(k, θ, p)
#     θ += ifelse(p.ξ==1, 0, π) # to account for the two valleys
#     a(p) * k^5 * cis(5θ) - 
#     b(p) * k^4 * cis(2θ) +
#     c(p) * k^3 * cis(-θ) +  # this two are the weaker
#     d(p) * k^2 * cis(2θ) - 
#     ϵ(p) * k * cis(-θ) # this two are the weaker
# end

# a(p) = p.γ0^5 / p.γ1^4 
# b(p) = 4p.γ0^3*p.γ3/p.γ1^3
# c(p) = 3p.γ0*p.γ3^2/p.γ1^2
# d(p) = 3/2 * p.γ2*p.γ0^2/p.γ1^2
# ϵ(p) = p.γ2*p.γ3/p.γ1
# y(p) = 2*(p.γ0*p.γ4/p.γ1)

# dhx5lg(k::Array, p) = dhx5lg(norm(k), atan(k[2], k[1]), p)
# dhy5lg(k::Array, p) = dhy5lg(norm(k), atan(k[2], k[1]), p)

# function dhx5lg(k::Number, θ, p)
#     θ += ifelse(p.ξ==1, 0, π)
#     m11 = 2*k*y(p)*cos(θ)
#     m22 = m11
#     m12 = m12x(k, θ, p)
#     return [m11 m12; conj(m12) m22]
# end

# function dhy5lg(k::Number, θ, p) 
#     θ += ifelse(p.ξ==1, 0, π)
#     m11 = 2*k*y(p)*sin(θ)
#     m22 = m11
#     m12 = m12y(k, θ, p)
#     return [m11 m12; conj(m12) m22]
# end

# m12y(k, θ, p) =  cos(θ-ifelse(p.ξ==1, 0, π)) * dkX5LG(k, θ, p) - 
#     sin(θ-ifelse(p.ξ==1, 0, π))/k * dθX5LG(k, θ, p) #parece que los ejes est'an girados por eso x->y
# m12x(k, θ, p) =  sin(θ-ifelse(p.ξ==1, 0, π)) * dkX5LG(k, θ, p) + 
#     cos(θ-ifelse(p.ξ==1, 0, π))/k * dθX5LG(k, θ, p)

# dkX5LG(k, θ, p) = 5a(p) * k^4 * cis(5θ) - 
#     4b(p) * k^3 * cis(2θ) +
#     3c(p) * k^2 * cis(-θ) +
#     2d(p) * k * cis(2θ) - 
#     ϵ(p)  * cis(-θ) 

# dθX5LG(k, θ, p) = 1im*(5a(p) * k^5 * cis(5θ) - 
#     2b(p) * k^4 * cis(2θ) +
#     -c(p) * k^3 * cis(-θ) +
#     2d(p) * k^2 * cis(2θ) - 
#     -ϵ(p) * k * cis(-θ))

# rz5lg(ψs) = rzNlg(5, ψs)

 
# TEST sigma yyy is 0. Success
#test to compute yyy ,mN12xx is really mN12yy i have changed x by y in teh velocities too # success :) sigma yyy is 0
# mN12x(n1,n2,n3, k, θ, p) =  dkcommon(n1,n2,n3, k, θ, p)  * (cos(θ-ifelse(p.ξ==1, 0, π)) *  (n1+n2) * cis(θ*(n1-n2)) - 
#     sin(θ-ifelse(p.ξ==1, 0, π)) * 1im *(n1-n2) * cis(θ*(n1-n2))) #parece que los ejes est'an girados por eso x->y
# mN12y(n1,n2,n3, k, θ, p) =  dkcommon(n1,n2,n3, k, θ, p) * p.ξ * (sin(θ-ifelse(p.ξ==1, 0, π)) * (n1+n2) * cis(θ*(n1-n2)) + 
#    cos(θ-ifelse(p.ξ==1, 0, π)) * 1im *(n1-n2) * cis(θ*(n1-n2)))
# dkcommon(n1,n2,n3, k, θ, p) = (prefactor(n1,n2,n3, p) * k^(n1+n2-1)) 
# mN12xx(n1,n2,n3, k, θ, p) =  cos(θ-ifelse(p.ξ==1, 0, π)) * mN12yyk(n1,n2,n3, k, θ, p) - sin(θ-ifelse(p.ξ==1, 0, π))/k * mN12yyθ(n1,n2,n3, k, θ, p)
# mN12yyk(n1,n2,n3, k, θ, p) = prefactor(n1,n2,n3,p) * cis(θ * (n1-n2)) *
#      (cos(θ-ifelse(p.ξ==1, 0, π)) * (n1+n2) - sin(θ-ifelse(p.ξ==1, 0, π)) * 1im * (n1-n2)) * (n1+n2-1) * k^(n1+n2-2)
# mN12yyθ(n1,n2,n3, k, θ, p) = prefactor(n1,n2,n3,p) * k^(n1+n2-1) * cis(θ * (n1-n2)) * (
#      (n1+n2) * (1im *(n1-n2)*cos(θ-ifelse(p.ξ==1, 0, π)) - sin(θ-ifelse(p.ξ==1, 0, π))) + 1im *(n1-n2)* (1im *(n1-n2)*sin(θ-ifelse(p.ξ==1, 0, π)) +cos(θ-ifelse(p.ξ==1, 0, π)))
# )
#test