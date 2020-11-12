abstract type NestedEnum end

function Base.:(<=)(a::T, b::T) where {T<:NestedEnum}
    return a.value % b.value == 0
end

function Base.:(<=)(a::T, b::U) where {T<:NestedEnum,U<:NestedEnum}
    return false
end

const PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 
    103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 
    367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 
    503, 509, 521, 523, 541]

# Store the item defined in expr at position i
function _store!(expr::Symbol, i, names, parent_pos, depths)
    names[i] = expr
    parent_pos[i] = 0 # No parent
    depths[i] = 0 # No parent
    return
end

# Store the item defined in expr at position i
function _store!(expr::Expr, i, names, parent_pos, depths)
    if i == 1 # parent can be a curly expression e.g. Duty{Variable}
        expr.head == :curly || error("Syntax error : parent can be a Symbol or a curly expression.")
        names[i] = expr
        parent_pos[i] = 0
        depths[i] = 0
        return
    end
    expr.head == :call || error("Syntax error :  Child <= Parent ")
    expr.args[1] == :(<=) || error("Syntax error : Child <= Parent ")
    i > 1 || error("First element cannot have a parent.")

    name = expr.args[2]
    parent_name = expr.args[3]

    r = findall(n -> n == parent_name, names[1:i-1])
    length(r) == 0 && error("Unknow parent $(parent_name).")
    length(r) > 1 && error("$(parent_name) registered more than once.")
    pos = r[1]
    parent_pos[i] = pos
    names[i] = name
    depths[i] = depths[pos] + 1
    return
end

# Compute the value of each item. The value is equal to the multiplication of 
# the prime numbers assigned to the item and its ancestors.
function _compute_values!(values, parent_pos, primes)
    for i in 1:length(parent_pos)
        factor = 1
        j = parent_pos[i]
        if j != 0
            factor = values[j]
        end
        values[i] = primes[i] * factor
    end
    return
end

# Update parent_pos array in function of permutation p
function _update_parent_pos!(parent_pos, p)
    permute!(parent_pos, p) # We still use positions of the old order.
    inv_p = invperm(p)
    for (i, pos) in enumerate(parent_pos)
        if pos != 0
            parent_pos[i] = inv_p[pos]
        end
    end
    return
end

function _build_expression(names, values, export_symb::Bool = false)
    len = length(names)
    root_name = names[1]
    enum_expr = Expr(:block, :())
    # We define a new type iif the root name is a Symbol
    # If the root name is a curly expression, the user must have defined the
    # template type inheriting from NestedEnum in its code.
    if root_name isa Symbol
        push!(enum_expr.args, :(struct $root_name <: Coluna.ColunaBase.NestedEnum value::UInt end))
    end
    for i in 2:len
        push!(enum_expr.args, :(const $(names[i]) = $(root_name)(UInt($(values[i])))))
        if export_symb
            push!(enum_expr.args, :(export $(names[i])))
        end
    end
    return enum_expr
end

function _build_print_expression(names, values)
    root_name = names[1]
    print_expr = Expr(:function)
    push!(print_expr.args, :(Base.print(io::IO, obj::$(root_name)))) #signature

    # build the if list in reverse order
    prev_cond = :(print(io, "UNKNOWN_DUTY"))
    for i in length(names):-1:2
        head = (i == 2) ? :if : :elseif
        msg = string(names[i])
        cond = Expr(head, :(obj == $(names[i])), :(print(io, $msg)), prev_cond)
        prev_cond = cond
    end
    push!(print_expr.args, Expr(:block, prev_cond))
    return print_expr
end

function _assign_values_to_items(expr)
    Base.remove_linenums!(expr)

    expr.head == :block || error("Block expression expected.")

    len = length(expr.args)
    names = Array{Union{Symbol, Expr}}(undef, len)
    parent_pos = zeros(Int, len) # Position of the parent.
    depths = zeros(Int, len) # Depth of each item
    values = zeros(UInt32, len) # The value is the multiplication of primes of the item and its ancestors.
    primes = PRIMES[1:len]

    name_values = Dict{Union{Symbol, Expr}, Int}() 
    for (i, arg) in enumerate(expr.args)
        _store!(arg, i, names, parent_pos, depths)
    end

    p = sortperm(depths)
    permute!(names, p)
    _update_parent_pos!(parent_pos, p)
    _compute_values!(values, parent_pos, primes)
    return names, values
end

"""

    @nestedenum block_expression

Create a `NestedEnum` subtype such as :

# Example

```@meta
DocTestSetup = quote
    using Coluna
end
```

```jldoctest nestedexample
Coluna.ColunaBase.@nestedenum begin 
    TypeOfItem
    ItemA <= TypeOfItem
        ChildA1 <= ItemA
            GrandChildA11 <= ChildA1
            GrandChildA12 <= ChildA1
        ChildA2 <= ItemA
    ItemB <= TypeOfItem
    ItemC <= TypeOfItem
end

# output

```

Create a nested enumeration with items `ItemA`, `ChildA1`, `ChildA2`, `GrandChildA11`, 
`GrandChildA12`, `ItemB`, and `ItemC` of type `TypeOfItem`.
The operator `<=` indicates the parent of the item.

```jldoctest nestedexample
julia> GrandChildA11 <= ItemA
true
```

```jldoctest nestedexample
julia> GrandChildA11 <= ItemC
false
```

```@meta
DocTestSetup = nothing
```
"""
macro nestedenum(expr)
    return _nestedenum(expr, false)
end

"Create a nested enumeration and export all the items."
macro exported_nestedenum(expr)
    return _nestedenum(expr, true)
end

function _nestedenum(expr, export_names)
    names, values = _assign_values_to_items(expr)
    enum_expr = _build_expression(names, values, export_names)
    print_expr = _build_print_expression(names, values)
    final_expr = quote
        $enum_expr
        $print_expr
    end
    return esc(final_expr)
end
