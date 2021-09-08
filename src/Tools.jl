"""
    function addpos(data::Array, pos::AbstractPoint)
    function addpos(data::StructArray, pos::AbstractPoint)

Add `pos` to `:Pos` of all particles
"""
function addpos(data::Array, pos::AbstractPoint)
    for i in 1:length(data)
        data[i] = setproperties!!(data[i], Pos = data[i].Pos + pos)
    end
end

function addpos(data::StructArray, pos::AbstractPoint)
    Pos = data.Pos
    for i in 1:length(data)
        Pos[i] += pos
    end
end

"""
    function addvel(data::Array, vel::AbstractPoint)
    function addvel(data::StructArray, vel::AbstractPoint)

Add `vel` to `:Vel` of all particles
"""
function addvel(data::Array, vel::AbstractPoint)
    for i in 1:length(data)
        data[i] = setproperties!!(data[i], Vel = data[i].Vel + vel)
    end
end

function addvel(data::StructArray, vel::AbstractPoint)
    Vel = data.Vel
    for i in 1:length(data)
        Vel[i] += vel
    end
end

"""
    function setpos(data::Union{Array, StructArray}, pos::AbstractPoint)

Set system center (middle value) to `pos`
"""
function setpos(data::Union{Array, StructArray}, pos::AbstractPoint)
    p0 = median(data, :Pos)
    addpos(data, pos - p0)
end

"""
    function setvel(data::Union{Array, StructArray}, vel::AbstractPoint)

Set system velocity (mass weighted average) to `vel`
"""
function setvel(data::Union{Array, StructArray}, vel::AbstractPoint)
    v0 = averagebymass(data, :Vel)
    addvel(data, vel - v0)
end