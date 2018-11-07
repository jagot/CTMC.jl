isdefined(Main, :ElectricField) || @derived_dimension ElectricField Unitful.𝐈^-1*Unitful.𝐋*Unitful.𝐌*Unitful.𝐓^-3
isdefined(Main, :Intensity) || @derived_dimension Intensity Unitful.𝐌*Unitful.𝐓^-3
isdefined(Main, :InvIntensity) || @derived_dimension InvIntensity Unitful.𝐌^-1*Unitful.𝐓^3
isdefined(Main, :Velocity) || @derived_dimension Velocity Unitful.𝐋*Unitful.𝐓^-1

bohr = 5.2917721067e-11u"m"
ec = 1.60217662e-19u"C"
kₑ = 8987551787.3681764u"N"*u"m"^2/(u"C"^2)

atomic_units(I::Intensity) = I/(3.5094452e16*u"W"/(u"cm"^2)) |> NoUnits
atomic_units(iI::InvIntensity) = iI*3.5094452e16*u"W"/(u"cm"^2) |> NoUnits
atomic_units(E::ElectricField) = E/(5.14220651e11*u"V"/u"m") |> NoUnits
atomic_units(t::Unitful.Time) = t/24.1888430u"as" |> NoUnits
atomic_units(Wk::Unitful.Energy) = Wk/27.211u"eV" |> NoUnits
atomic_units(l::Unitful.Length) = l/5.2917721067e-11u"m" |> NoUnits
atomic_units(v::Velocity) = v/(2.1876912633e6u"m"/u"s") |> NoUnits

export bohr, ec
