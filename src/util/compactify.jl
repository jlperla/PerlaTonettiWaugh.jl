struct Compactifier
    transformer::ContinuousTransformations.ComposedTransformation
    Compactifier(lb::Float64, ub::Float64) = lb >= ub ? error("ub should be strictly greater than lb.") : new(bridge(ℝ, Segment(lb, ub), RealCircle()))
    Compactifier(lb::Float64, ub::Float64, transformer_type::ContinuousTransformations.UnivariateTransformation) = lb >= ub ? error("ub should be strictly greater than lb.") : new(bridge(ℝ, Segment(lb, ub), transformer_type))
    (f::Compactifier)(x) = f.transformer(x)
end

# define Decompactifier
struct Decompactifier
    transformer::ContinuousTransformations.ComposedTransformation
    Decompactifier(compactifier::Compactifier) = new(ContinuousTransformations.inverse(compactifier.transformer))
end
(f::Decompactifier)(y) = f.transformer(y)