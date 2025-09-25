# Submitted as a patch to StatisticalMeasures.jl:
#    https://github.com/JuliaAI/StatisticalMeasures.jl/pull/53

# Use this file until the patch is merged.

module Willmott

using StatisticalMeasures: API, @trait, @fix_show, register, LearnAPI,
    Infinite, Score, DOC_INFINITE, docstring, aggregate, multimeasure,
    LPSumLoss, Sum

# -------------------------------------------------------------------------
# Willmott index of agreement (d)

# type for measure without argument checks:
struct _WillmottD end

function (::_WillmottD)(yhat, y)
    μ = aggregate(y)  # mean
    # numerator: Σ_i (ŷ_i - y_i)^2
    num = LPSumLoss(p=2)(yhat, y)
    # denominator: Σ_i (|ŷ_i - μ| + |y_i - μ|)^2
    den = multimeasure((yhat, y) -> (abs(yhat - μ) + abs(y - μ))^2; mode=Sum())(yhat, y)
    return den == 0 ? (num == 0 ? 1.0 : 0.0) : 1 - num/den
end

WillmottD() = _WillmottD() |> API.robust_measure |> API.fussy_measure
const WillmottDType = API.FussyMeasure{<:API.RobustMeasure{<:_WillmottD}}

@trait(
    _WillmottD,
    consumes_multiple_observations = true,
    kind_of_proxy = LearnAPI.Point(),
    observation_scitype = Union{Missing,Infinite},
    orientation = Score(),
    human_name = "Willmott index of agreement (d)",
)

@fix_show WillmottD::WillmottDType

register(WillmottD, "willmott_d")

const WillmottDDoc = """
Return Willmott index of agreement (d)

``d = 1 - \\dfrac{\\sum (ŷ_i - y_i)^2}{\\sum (|ŷ_i - \\bar y| + |y_i - \\bar y|)^2}``,

where ``\\bar y`` is the mean of the targets. The value lies in ``[0,1]`` with higher
being better.

References: Willmott (1981) https://doi.org/10.1080/02723646.1981.10642213
"""

"$WillmottDDoc"
WillmottD
"$WillmottDDoc"
const willmott_d = WillmottD()

end # module
