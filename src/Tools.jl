"""
    function addpos(data::Array, pos::AbstractPoint)

Add `pos` to `:Pos` of all particles
"""
function addpos(data::Array, pos::AbstractPoint)
    for i in 1:length(data)
        data[i] = setproperties!!(data[i], Pos = data[i].Pos + pos)
    end
end

"""
    function addvel(data::Array, vel::AbstractPoint)

Add `vel` to `:Vel` of all particles
"""
function addvel(data::Array, vel::AbstractPoint)
    for i in 1:length(data)
        data[i] = setproperties!!(data[i], Vel = data[i].Vel + vel)
    end
end

"""
    function setpos(data, pos::AbstractPoint)

Set system center (middle value) to `pos`
"""
function setpos(data::Array, pos::AbstractPoint)
    p0 = median(data, :Pos)
    addpos(data, pos - p0)
end

function setpos(data::Array, pos::AbstractPoint)
    p0 = median(data, :Pos)
    for k in keys(data)
        addpos(data, pos - p0)
    end
end

"""
    function setvel(data, vel::AbstractPoint)

Set system velocity (mass weighted average) to `vel`
"""
function setvel(data::Array, vel::AbstractPoint)
    v0 = averagebymass(data, :Vel)
    addvel(data, vel - v0)
end

function setvel(data::Dict, vel::AbstractPoint)
    v0 = averagebymass(data, :Vel)
    for k in keys(data)
        addvel(data, vel - v0)
    end
end

