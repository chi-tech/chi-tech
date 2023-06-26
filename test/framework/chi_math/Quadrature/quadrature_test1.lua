function PrintTable(t, indent)
  if not indent then indent = 0 end
  strform = "%"..tostring(indent).."s"

  for k,v in pairs(t) do
    if (type(v) == "table") then
      print(string.rep(" ", indent)..k.." ".."table")
      PrintTable(v, indent+2)
    else
      print(string.rep(" ", indent)..k.." "..tostring(v))
    end
  end
end

print("GOLD_BEGIN")
q = chi_math.QuadratureGaussLegendre.Create({N = 4, verbose = true})

qdata = chi_math.Get1DQuadratureData(q)

print("qpoints:")
PrintTable(qdata.qpoints, 2)
print("weights:")
PrintTable(qdata.weights, 2)
print()

--################################################
q = chi_math.QuadratureGaussChebyshev.Create({N = 4, verbose = true})

qdata = chi_math.Get1DQuadratureData(q)

print("qpoints:")
PrintTable(qdata.qpoints, 2)
print("weights:")
PrintTable(qdata.weights, 2)

print("GOLD_END")