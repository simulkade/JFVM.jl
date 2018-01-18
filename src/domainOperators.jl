# cell value math functions
for f in [:sin, :cos, :tan, :cot, :abs, :exp, :log, :log10]
  @eval function $f(a::CellValue)
      return CellValue(a.domain, $f.(a.value))
  end
end

# face value math functions
for f in [:sin, :cos, :tan, :cot, :abs, :exp, :log, :log10]
  @eval function $f(a::FaceValue)
      return FaceValue(a.domain, 
                  $f.(a.xvalue),
                  $f.(a.yvalue),
                  $f.(a.zvalue))
  end
end

# cell value operators
for f in [:+, :-, :*, :/, :^, :(==), :>, :(>=), :<, :(<=)]
  @eval function $f(a::CellValue, b::CellValue)
      return CellValue(a.domain, $f.(a.value, b.value))
  end
end

for f in [:+, :-, :*, :/, :^, :(==), :>, :(>=), :<, :(<=)]
  @eval function $f(a::CellValue, b::Real)
      return CellValue(a.domain, $f.(a.value, b))
  end
end

for f in [:+, :-, :*, :/, :^, :(==), :>, :(>=), :<, :(<=)]
  @eval function $f(a::Real, b::CellValue)
      return CellValue(b.domain, $f.(a, b.value))
  end
end

# face value operators
for f in [:+, :-, :*, :/, :^, :(==), :>, :(>=), :<, :(<=)]
  @eval function $f(a::FaceValue, b::FaceValue)
      return FaceValue(a.domain, 
        $f.(a.xvalue, b.xvalue),
        $f.(a.yvalue, b.yvalue),
        $f.(a.zvalue, b.zvalue))
  end
end

for f in [:+, :-, :*, :/, :^, :(==), :>, :(>=), :<, :(<=)]
  @eval function $f(a::FaceValue, b::Real)
      return FaceValue(a.domain, 
        $f.(a.xvalue, b),
        $f.(a.yvalue, b),
        $f.(a.zvalue, b))
  end
end

for f in [:+, :-, :*, :/, :^, :(==), :>, :(>=), :<, :(<=)]
  @eval function $f(a::Real, b::FaceValue)
      return FaceValue(b.domain, 
        $f.(a, b.xvalue),
        $f.(a, b.yvalue),
        $f.(a, b.zvalue))
  end
end

function -(a::FaceValue)
  FaceValue(a.domain,
    -a.xvalue,
    -a.yvalue,
    -a.zvalue)
end

function -(a::CellValue)
  CellValue(a.domain,
    -a.value)
end

# # cell value operators
# for f in [:+, :-, :*, :/, :^, :(==), :>, :(>=), :<, :(<=)]
#   @eval function $f(a::CellValue, b::CellValue)
#       return CellValue(a.domain, $f.(a.value, b.value))
#   end
# end

# function broadcast(::typeof(*), a::FaceValue, b::FaceValue)
#   FaceValue(a.domain,
#     a.xvalue.*b.xvalue,
#     a.yvalue.*b.yvalue,
#     a.zvalue.*b.zvalue)
# end
  
# function broadcast(::typeof(/), a::FaceValue, b::FaceValue)
#   FaceValue(a.domain,
#     a.xvalue./b.xvalue,
#     a.yvalue./b.yvalue,
#     a.zvalue./b.zvalue)
# end

# function broadcast(::typeof(*), a::CellValue, b::CellValue)
# CellValue(a.domain,
#   a.value.*b.value)
# end

# function broadcast(::typeof(/), a::CellValue, b::CellValue)
# CellValue(a.domain,
#   a.value./b.value)
# end

# function broadcast(::typeof(*), a::CellVector, b::CellVector)
# CellVector(a.domain,
#   a.xvalue.*b.xvalue,
#   a.yvalue.*b.yvalue,
#   a.zvalue.*b.zvalue)
# end

# function broadcast(::typeof(/), a::CellVector, b::CellVector)
# CellVector(a.domain,
#   a.xvalue./b.xvalue,
#   a.yvalue./b.yvalue,
#   a.zvalue./b.zvalue)
# end


## To be thrown away

# function <(a::FaceValue, b::FaceValue)
#   FaceValue(a.domain,
#     a.xvalue.<b.xvalue,
#     a.yvalue.<b.yvalue,
#     a.zvalue.<b.zvalue)
# end

# function <=(a::FaceValue, b::FaceValue)
#   FaceValue(a.domain,
#     a.xvalue.<=b.xvalue,
#     a.yvalue.<=b.yvalue,
#     a.zvalue.<=b.zvalue)
# end

# function ==(a::FaceValue, b::FaceValue)
#   FaceValue(a.domain,
#     a.xvalue.==b.xvalue,
#     a.yvalue.==b.yvalue,
#     a.zvalue.==b.zvalue)
# end

# function >(a::FaceValue, b::FaceValue)
#   FaceValue(a.domain,
#     a.xvalue.>b.xvalue,
#     a.yvalue.>b.yvalue,
#     a.zvalue.>b.zvalue)
# end

# function >=(a::FaceValue, b::FaceValue)
#   FaceValue(a.domain,
#     a.xvalue.>=b.xvalue,
#     a.yvalue.>=b.yvalue,
#     a.zvalue.>=b.zvalue)
# end

# function +(a::FaceValue, b::FaceValue)
# FaceValue(a.domain,
#   a.xvalue+b.xvalue,
#   a.yvalue+b.yvalue,
#   a.zvalue+b.zvalue)
# end

# function *(a::FaceValue, b::FaceValue)
# FaceValue(a.domain,
#   a.xvalue.*b.xvalue,
#   a.yvalue.*b.yvalue,
#   a.zvalue.*b.zvalue)
# end

# function /(a::FaceValue, b::FaceValue)
# FaceValue(a.domain,
#   a.xvalue./b.xvalue,
#   a.yvalue./b.yvalue,
#   a.zvalue./b.zvalue)
# end

# function -(a::FaceValue, b::FaceValue)
# FaceValue(a.domain,
#   a.xvalue-b.xvalue,
#   a.yvalue-b.yvalue,
#   a.zvalue-b.zvalue)
# end

# 
# function +(a::Real, b::FaceValue)
# FaceValue(b.domain,
#   a.+b.xvalue,
#   a.+b.yvalue,
#   a.+b.zvalue)
# end

# function +(a::FaceValue, b::Real)
# FaceValue(a.domain,
#   b.+a.xvalue,
#   b.+a.yvalue,
#   b.+a.zvalue)
# end

# function -(a::Real, b::FaceValue)
# FaceValue(b.domain,
#   a.-b.xvalue,
#   a.-b.yvalue,
#   a.-b.zvalue)
# end

# function -(a::FaceValue)
# FaceValue(a.domain,
#   -a.xvalue,
#   -a.yvalue,
#   -a.zvalue)
# end

# function -(a::FaceValue, b::Real)
# FaceValue(a.domain,
#   a.xvalue.-b,
#   a.yvalue.-b,
#   a.zvalue.-b)
# end

# function *(a::Real, b::FaceValue)
# FaceValue(b.domain,
#   a.*b.xvalue,
#   a.*b.yvalue,
#   a.*b.zvalue)
# end

# function *(a::FaceValue, b::Real)
# FaceValue(a.domain,
#   b.*a.xvalue,
#   b.*a.yvalue,
#   b.*a.zvalue)
# end

# function /(a::Real, b::FaceValue)
# FaceValue(b.domain,
#   a./b.xvalue,
#   a./b.yvalue,
#   a./b.zvalue)
# end

# function /(a::FaceValue, b::Real)
# FaceValue(a.domain,
#   a.xvalue./b,
#   a.yvalue./b,
#   a.zvalue./b)
# end

# # Cell Value operators
# function ==(a::CellValue, b::CellValue)
#   CellValue(a.domain,
#     a.value.==b.value)
# end
  
# function <(a::CellValue, b::CellValue)
#   CellValue(a.domain,
#     a.value.<b.value)
# end

# function <=(a::CellValue, b::CellValue)
#   CellValue(a.domain,
#     a.value.<=b.value)
# end

# function >(a::CellValue, b::CellValue)
#   CellValue(a.domain,
#     a.value.>b.value)
# end

# function >=(a::CellValue, b::CellValue)
#   CellValue(a.domain,
#     a.value.>=b.value)
# end

# function +(a::CellValue, b::CellValue)
# CellValue(a.domain,
#   a.value+b.value)
# end

# function *(a::CellValue, b::CellValue)
# CellValue(a.domain,
#   a.value.*b.value)
# end

# function /(a::CellValue, b::CellValue)
# CellValue(a.domain,
#   a.value./b.value)
# end

# function -(a::CellValue, b::CellValue)
# CellValue(a.domain,
#   a.value-b.value)
# end



# function +(a::Real, b::CellValue)
# CellValue(b.domain,
#   a.+b.value)
# end

# function +(a::CellValue, b::Real)
# CellValue(a.domain,
#   b.+a.value)
# end

# function -(a::Real, b::CellValue)
# CellValue(b.domain,
#   a-b.value)
# end

# function -(a::CellValue)
# CellValue(a.domain,
#   -a.value)
# end

# function -(a::CellValue, b::Real)
# CellValue(a.domain,
#   a.value-b)
# end

# function *(a::Real, b::CellValue)
# CellValue(b.domain,
#   a*b.value)
# end

# function *(a::CellValue, b::Real)
# CellValue(a.domain,
#   b*a.value)
# end

# function /(a::Real, b::CellValue)
# CellValue(b.domain,
#   a./b.value)
# end

# function /(a::CellValue, b::Real)
# CellValue(a.domain,
#   a.value/b)
# end

# # Cell Vector operators
# function +(a::CellVector, b::CellVector)
# CellVector(a.domain,
#   a.xvalue+b.xvalue,
#   a.yvalue+b.yvalue,
#   a.zvalue+b.zvalue)
# end

# function -(a::CellVector, b::CellVector)
# CellVector(a.domain,
#   a.xvalue-b.xvalue,
#   a.yvalue-b.yvalue,
#   a.zvalue-b.zvalue)
# end



# function +(a::Real, b::CellVector)
# CellVector(b.domain,
#   a.+b.xvalue,
#   a.+b.yvalue,
#   a.+b.zvalue)
# end

# function +(a::CellVector, b::Real)
# CellVector(a.domain,
#   b.+a.xvalue,
#   b.+a.yvalue,
#   b.+a.zvalue)
# end

# function -(a::Real, b::CellVector)
# CellVector(b.domain,
#   a.-b.xvalue,
#   a.-b.yvalue,
#   a.-b.zvalue)
# end

# function -(a::CellVector)
# CellVector(a.domain,
#   -a.xvalue,
#   -a.yvalue,
#   -a.zvalue)
# end

# function -(a::CellVector, b::Real)
# CellVector(a.domain,
#   a.xvalue-b,
#   a.yvalue-b,
#   a.zvalue-b)
# end

# function *(a::Real, b::CellVector)
# CellVector(b.domain,
#   a*b.xvalue,
#   a*b.yvalue,
#   a*b.zvalue)
# end

# function *(a::CellVector, b::Real)
# CellVector(a.domain,
#   b*a.xvalue,
#   b*a.yvalue,
#   b*a.zvalue)
# end

# function /(a::Real, b::CellVector)
# CellVector(b.domain,
#   a./b.xvalue,
#   a./b.yvalue,
#   a./b.zvalue)
# end

# function /(a::CellVector, b::Real)
# CellVector(a.domain,
#   a.xvalue/b,
#   a.yvalue/b,
#   a.zvalue/b)
# end
