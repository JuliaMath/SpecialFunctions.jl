module SpecialFunctionsSymbolicsExt

using Symbolics, SpecialFunctions

@variables a x
Symbolics.@register_symbolic gamma(a, x)
Symbolics.@register_symbolic gamma_inc(a, x)
Symbolics.@register_symbolic lowergamma(a, x)
Symbolics.@register_symbolic uppergamma(a, x)

end