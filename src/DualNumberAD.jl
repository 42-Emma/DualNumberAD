module DualNumberAD

export Dual, deriv1, deriv2

struct Dual{T}
  value::T
  delta::T
end

# Operator Overloading

#Addition
Base.:+(a::Dual, b::Dual) = Dual(a.value + b.value, a.delta + b.delta) # Dual + Dual
Base.:+(a::Dual, m::Number) = Dual(a.value + m, a.delta) # Dual + Number
Base.:+(m::Number, a::Dual) = Dual(m + a.value, a.delta) # Number + Dual

# Subtraction
Base.:-(a::Dual, b::Dual) = Dual(a.value - b.value, a.delta - b.delta) # Dual - Dual
Base.:-(a::Dual, m::Number) = Dual(a.value - m, a.delta) # Dual - Number
Base.:-(m::Number, a::Dual) = Dual(m - a.value, -a.delta) # Number - Dual
Base.:-(a::Dual) = Dual(-a.value, -a.delta) # unary minus: -Dual

#Multiplication
Base.:*(a::Dual, b::Dual) = Dual(
  a.value * b.value,
  a.delta * b.value + a.value * b.delta
) # Dual * Dual (product rule)
Base.:*(a::Dual, m::Number) = Dual(a.value * m, a.delta * m) # Dual * Number
Base.:*(m::Number, a::Dual) = Dual(m * a.value, m * a.delta) # Number * Dual

# Division
Base.:/(a::Dual, b::Dual) = Dual(
    a.value / b.value,
    (a.delta * b.value - a.value * b.delta) / (b.value * b.value)
) # Dual / Dual (Quotient rule)
Base.:/(a::Dual, m::Number) = Dual(a.value / m, a.delta / m) # Dual / Number
Base.:/(m::Number, a::Dual) = Dual(
    m / a.value,
    (-m * a.delta) / (a.value * a.value)
) # Number / Dual

# Powers
Base.:^(a::Dual, n::Number) = Dual(
    a.value^n,
    n * a.value^(n - 1) * a.delta
) # Dual to a real number
Base.:^(a::Number, b::Dual) = Dual(
    a^b.value,
    a^b.value * log(a) * b.delta
) # real to a dual power

#Elementary Functions
# Sine
Base.sin(a::Dual) = Dual(sin(a.value), cos(a.value) * a.delta) # sin(D) = sin(x) + cos(x)*δ

#Cos
Base.cos(a::Dual) = Dual(cos(a.value), -sin(a.value) * a.delta) # cos(D) = cos(x) - sin(x)*δ ε

#Exponential
Base.exp(a::Dual) = Dual(exp(a.value), exp(a.value) * a.delta) # exp(D) = exp(x) + exp(x)*δ ε

#Square root
Base.sqrt(a::Dual) = Dual(sqrt(a.value), (a.delta) / (2 * sqrt(a.value))) # sqrt(D) = sqrt(x) + (1/(2*sqrt(x)))*δ ε

# First derivative function
deriv1(f, x0::Number) = f(Dual(x0, one(x0))).delta

# Second Derivative Function
function deriv2(f, x0::Number)
    x = Dual(Dual(x0, one(x0)), Dual(one(x0), zero(x0)))
    y = f(x)
    return y.delta.delta
end

end # module DualNumberAD
