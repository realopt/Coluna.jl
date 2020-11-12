```@meta
CurrentModule = Coluna.MathProg
DocTestSetup = quote
    using Coluna.MathProg
end
```

# Formulation

```@docs
Formulation
```

## Duties of formulations

```@docs
Original
DwMaster
BendersMaster
DwSp
BendersSp
```

## Attributes of formulations

```@docs
haskey
getvar
getconstr
getvars
getconstrs
```

## Variables

```@docs
setvar!
```

## Constraints

```@docs
setconstr!
```

## Duties of variables and constraints

```@docs
Duty
```

## Attributes of variables and constraints

*Performance note* : use a variable or a constraint rather than its id.

```@docs
getperencost
getcurcost
setcurcost!
getperenlb
getcurlb
setcurlb!
getperenub
getcurub
setcurub!
getperenrhs
getcurrhs
setcurrhs!
getperenkind
getcurkind
setcurkind!
getperensense
getcursense
setcursense!
getperenincval
getcurincval
setcurincval!
isperenactive
iscuractive
activate!
deactivate!
isexplicit
getname
```

```@meta
CurrentModule = nothing
DocTestSetup = nothing
```