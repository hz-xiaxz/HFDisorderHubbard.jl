# Consider using BloqadeLattice


# Maybe use BloqadeLattice as backend for it's more easy to represent Rydberg atoms in a lattice
# does it have convient support for pbc lattice? leave it as TODO

struct CubicLattice
    L::Int
    N::Int
    neigh::AbstractVector
    function CubicLattice(L::Int)
        N = L^3
        neigh = Vector{NTuple{6,Int64}}(undef, N)
        cubic_lat = reshape(1:N, L, L, L)
        left = circshift(cubic_lat, (-1, 0, 0))
        right = circshift(cubic_lat, (1, 0, 0))
        front = circshift(cubic_lat, (0, 1, 0))
        back = circshift(cubic_lat, (0, -1, 0))
        up = circshift(cubic_lat, (0, 0, 1))
        down = circshift(cubic_lat, (0, 0, -1))
        for i = 1:N
            neigh[i] = (left[i], right[i], front[i], back[i], up[i], down[i])
        end
        return new(L, N, neigh)
    end
end
