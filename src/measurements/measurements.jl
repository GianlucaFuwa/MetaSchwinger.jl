module Measurements
    using Arpack
    using LinearAlgebra
    using Printf
    using SparseArrays

    import ..DiracOperators: AbstractDiracOperator
    import ..Gaugefields: Gaugefield, plaquette_sum
    import ..Metadynamics: BiasPotential
    import ..Observables: poly_loop_avg, topological_charge, wilson_loop_all

    defaultmeasures = Array{Dict,1}(undef, 2)

    for i in 1:length(defaultmeasures)
        defaultmeasures[i] = Dict()
    end

	defaultmeasures[1]["methodname"] = "Meta_charge"
    defaultmeasures[1]["measure_every"] = 1
    defaultmeasures[2]["methodname"] = "Topological_charge"
    defaultmeasures[2]["measure_every"] = 1

    struct MeasurementSet
        nummeas::Int64
        D::Union{Nothing, AbstractDiracOperator}
        meas_calls::Array{Dict,1}
        meas_files::Array{IOStream,1}

        function MeasurementSet(
            measure_dir;
            meas_calls = defaultmeasures,
            D = nothing,
            instance = "",
        )
            nummeas = length(meas_calls)
            meas_files = Array{IOStream,1}(undef, nummeas)

            for i in 1:nummeas
                method = meas_calls[i]

                meas_overwrite = "w"
                if method["methodname"] == "Plaquette"
                    meas_files[i] = open(
                        measure_dir * "/Plaquette" * instance * ".txt",
                        meas_overwrite
                    )
                    println(meas_files[i], "itrj\tRe(plaq)")
                elseif method["methodname"] == "Action"
                    meas_files[i] = open(
                        measure_dir * "/Action" * instance * ".txt",
                        meas_overwrite
                    )
                    println(meas_files[i], "itrj\tAction")
                elseif method["methodname"] == "Meta_charge"
                    meas_files[i] = open(
                        measure_dir * "/Meta_charge" * instance * ".txt",
                        meas_overwrite
                    )
                    println(meas_files[i], "itrj\tQmeta")
                elseif method["methodname"] == "Topological_charge"
                    meas_files[i] = open(
                        measure_dir * "/Topological_charge" * instance * ".txt",
                        meas_overwrite
                    )
                    println(meas_files[i], "itrj\tQtop")
                elseif method["methodname"] == "Polyakov_loop"
                    meas_files[i] = open(
                        measure_dir * "/Polyakov_loop" * instance * ".txt",
                        meas_overwrite
                    )
                    println(
                        meas_files[i],
                        "itrj\t$(rpad("Re(poly)", 21, " "))\t$(rpad("Im(poly)", 21, " "))",
                    )
                elseif method["methodname"] == "Dirac_eigenvalues"
                    @assert D !== nothing "Dirac operator has to be specified for EVs"
                    meas_files[i] = open(
                        measure_dir * "/Dirac_eigenvalues" * instance * ".txt",
                        meas_overwrite
                    )
                    println(meas_files[i], "itrj\teigs")
                elseif method["methodname"] == "Dirac_determinant"
                    @assert D !== nothing "Dirac operator has to be specified for det"
                    meas_files[i] = open(
                        measure_dir * "/Dirac_determinant" * instance * ".txt",
                        meas_overwrite
                    )
                    println(
                        meas_files[i],
                        "itrj\t$(rpad("Re(logdetD)", 21, " "))\t$(rpad("Im(logdetD)", 21, " "))",
                    )
                elseif method["methodname"] == "Chiral_condensate"
                    @assert D !== nothing "Dirac operator has to be specified for CC"
                    meas_files[i] = open(
                        measure_dir * "/Chiral_condensate" * instance * ".txt",
                        meas_overwrite
                    )
                    println(
                        meas_files[i],
                        "itrj\t$(rpad("Re(CC)", 21, " "))\t$(rpad("Im(CC)", 21, " "))",
                    )
                elseif occursin("Wilson_loop", method["methodname"])
                    meas_files[i] = open(
                        measure_dir * "/" * method["methodname"] * instance * ".txt",
                        meas_overwrite
                    )
                    println(meas_files[i], "itrj\tWLoops(1 -> NX)")
                else
                    error("$(method["methodname"]) is not supported")
                end
            end

            return new(
                nummeas,
                D,
                meas_calls,
                meas_files,
            )
        end
    end

    function measurements(itr, field::Gaugefield, measset::MeasurementSet)
        has_been_calculated = false
        D = measset.D

        for i in 1:measset.nummeas
            method = measset.meas_calls[i]
            measfile = measset.meas_files[i]
            if itr % method["measure_every"] == 0
                if method["methodname"] == "Plaquette"
                    plaq = plaquette_sum(field)
                    plaq /= field.NV
                    plaqstr = @sprintf("%.15E", plaq)
                    println(measfile, "$(rpad("$itr", 5, " "))\t$plaqstr")
                elseif method["methodname"] == "Action"
                    s = field.Sg / field.NV
                    sstr = @sprintf("%.15E", s)
                    println(measfile, "$(rpad("$itr", 5, " "))\t$sstr")
                elseif method["methodname"] == "Meta_charge"
                    q = field.CV
                    qstr = @sprintf("%.15E", q)
                    println(measfile, "$(rpad("$itr", 5, " "))\t$qstr")
                elseif method["methodname"] == "Topological_charge"
                    qt = topological_charge(field)
                    println(measfile, "$(rpad("$itr", 5, " "))\t$qt")
                elseif method["methodname"] == "Polyakov_loop"
                    poly_re, poly_im = poly_loop_avg(field)
                    poly_restr = @sprintf("%.15E", poly_re)
                    poly_imstr = @sprintf("%.15E", poly_im)
                    println(measfile, "$(rpad("$itr", 5, " "))\t$poly_restr\t$poly_imstr")
                elseif method["methodname"] == "Dirac_eigenvalues"
                    if has_been_calculated == false
                        D(field)
                        Ds = sparse(D.Dop)
                    end

                    vals, _ = eigs(
                        Ds,
                        nev = method["nev"],
                        which = method["which"],
                        maxiter = method["maxiter"],
                    )
                    print(measfile, "$(rpad("$itr", 5, " "))")

                    for λ in vals
                        λrestr = @sprintf("%.15E", real(λ))
                        λimstr = @sprintf("%.15E", imag(λ))
                        print(measfile, "\t$(λrestr) $(λimstr)")
                    end

                    println(measfile, "")
                elseif method["methodname"] == "Dirac_determinant"
                    if has_been_calculated == false
                        D(field)
                        Ds = sparse(D.Dop)
                    end

                    logdetD = logdet(Ds)
                    Drestr = @sprintf("%.15E", real(logdetD))
                    Dimstr = @sprintf("%.15E", imag(logdetD))
                    println(measfile, "$(rpad("$itr", 5, " "))\t$Drestr\t$Dimstr")
                elseif method["methodname"] == "Chiral_condensate"
                    if has_been_calculated == false
                        D(field)
                    end

                    cc = tr(inv(D.Dop)) / field.NV
                    ccrestr = @sprintf("%.15E", real(cc))
                    ccimstr = @sprintf("%.15E", imag(cc))
                    println(measfile, "$(rpad("$itr", 5, " "))\t$ccrestr\t$ccimstr")
                elseif occursin( "Wilson_loop", method["methodname"])
                    num = filter.(isdigit, method["methodname"])
                    LT = parse(Int64, num)
                    wils_re, wils_im = wilson_loop_all(field, LT)
                    wils_re = round.(wils_re, sigdigits = 4)
                    wils_im = round.(wils_im, sigdigits = 4)
                    println(measfile, rpad("$itr", 5, " "), " ", wils_re, " ", wils_im)
                else
                    error("$(method["methodname"]) is not supported")
                end
                flush(measfile)
            end
        end
        return nothing
    end

    function calc_weights(q_vals::Vector{Float64}, b::BiasPotential)
        weights = zeros(length(q_vals))

        for (i, q) in enumerate(q_vals)
            V = b(q)
            weights[i] = exp(V)
        end

        return weights
    end

end
