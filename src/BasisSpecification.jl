module BasisSpecification

using ..GraphSignal, ..GraphPartition, LinearAlgebra

export BasisSpec #, bsfull, bs_haar, bs_level

# BasisSpec data structure and constructors
"""
    BS = BasisSpec(dvec, dvec_loc, c2f, description)

is a data structure for a BasisSpec object contaning the following fields:
* `levlist::Vector{Tuple{Int, Int}}`: the integer sequence that specifies a particular basis
* `c2f::Bool`: if true (default), this indicates that the basis comes from the coarse-to-fine dictionary; if false, this indicates that the basis comes from the fine-to-coarse dictionary
* `description::String`: a description of the specified basis


Copyright 2015 The Regents of the University of California

Implemented by Jeff Irion (Adviser: Dr. Naoki Saito) |
Translated to julia and revised by Naoki Saito, Feb. 22, 2017
Modified by Yiqun Shao, May. 20, 2018
"""

mutable struct BasisSpec
    levlist::Vector{Tuple{Int, Int}}
    c2f::Bool
    description::String

    # An inner constructor here.
    function BasisSpec(levlist::Vector{Tuple{Int, Int}};
                       c2f::Bool = true,
                       description::String = "")
        new(levlist, c2f, description )
    end # of BasisSpec constructor

end # of type BasisSpec

end  # of module BasisSpecification
