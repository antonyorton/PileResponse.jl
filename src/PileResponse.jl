module PileResponse

export
    helloworld,
    square,
    smooth_data


include("functions.jl")


function helloworld()
    print("Hi there PileResponse.js world, let's do this again.\n")
    return 123
end

end
