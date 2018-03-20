module BasisSpecification

using ..GraphSignal, ..GraphPartition

export BasisSpec #, bsfull, bs_haar, bs_level

# BasisSpec data structure and constructors
"""
    BS = BasisSpec(levlist, levlengths, c2f, description)

is a data structure for a BasisSpec object contaning the following fields:
* `levlist::Vector{UInt8}`: the integer sequence that specifies a particular basis
* `levlengths::Vector{UInt16}`: the integer sequence that specifies the length of each basis block in `levlist` (optional)
* `c2f::Bool`: if true (default), this indicates that the basis comes from the coarse-to-fine dictionary; if false, this indicates that the basis comes from the fine-to-coarse dictionary
* `description::String`: a description of the specified basis


Copyright 2015 The Regents of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito) |
Translated to julia and revised by Naoki Saito, Feb. 22, 2017
"""
type BasisSpec
    levlist :: Vector{Int}
    levlengths :: Vector{Int}
    c2f::Bool
    description::String

    # An inner constructor here.
    function BasisSpec(levlist :: Vector{Int};
                       levlengths :: Vector{Int} = Vector{Int}(0),
                       c2f::Bool = true,
                       description::String = "")

        # Sanity check here.
        if Base.length(levlengths) != Base.length(levlist) && Base.length(levlengths) != 0
            warn("length(levlengths) != length(levlist). Hence, use null levlengths")
            new( levlist, Vector{Integer}(0), c2f, description )
        else
            new( levlist, levlengths, c2f, description )
        end
    end # of BasisSpec constructor

end # of type BasisSpec

end # of module BasisSpecification
