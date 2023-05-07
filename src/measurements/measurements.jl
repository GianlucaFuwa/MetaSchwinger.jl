module Measurements
    using Arpack
    using LinearAlgebra
    using SparseArrays

    import ..DiracOperatorModule: AbstractDiracOperator, NaiveDiracOperator
    import ..DiracOperatorModule: WilsonDiracOperator
    import ..Gaugefields: Gaugefield, plaquette, plaquette_sum
    import ..Metadynamics: BiasPotential
    import ..Observables: meta_charge, poly_loop_avg, topological_charge, wilson_loop_all

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
                elseif method["methodname"] == "Action"
                    meas_files[i] = open(
                        measure_dir * "/Action" * instance * ".txt",
                        meas_overwrite
                    )
                elseif method["methodname"] == "Meta_charge"
                    meas_files[i] = open(
                        measure_dir * "/Meta_charge" * instance * ".txt",
                        meas_overwrite
                    )
                elseif method["methodname"] == "Topological_charge"
                    meas_files[i] = open(
                        measure_dir * "/Topological_charge" * instance * ".txt",
                        meas_overwrite
                    )
                elseif method["methodname"] == "Polyakov_loop"
                    meas_files[i] = open(
                        measure_dir * "/Polyakov_loop" * instance * ".txt",
                        meas_overwrite
                    )
                elseif method["methodname"] == "Dirac_eigenvalues"
                    meas_files[i] = open(
                        measure_dir * "/Dirac_eigenvalues" * instance * ".txt",
                        meas_overwrite
                    )
                elseif method["methodname"] == "Dirac_determinant"
                    meas_files[i] = open(
                        measure_dir * "/Dirac_determinant" * instance * ".txt",
                        meas_overwrite
                    )
                elseif method["methodname"] == "Chiral_condensate"
                    meas_files[i] = open(
                        measure_dir * "/Chiral_condensate" * instance * ".txt",
                        meas_overwrite
                    )
                elseif occursin("Wilson_loop", method["methodname"])
                    meas_files[i] = open(
                        measure_dir * "/" * method["methodname"] * instance * ".txt",
                        meas_overwrite
                    )
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
        Du = nothing
        Ds = nothing

        for i in 1:measset.nummeas
            method = measset.meas_calls[i]
            measfile = measset.meas_files[i]
            if itr % method["measure_every"] == 0
                if method["methodname"] == "Plaquette"
                    plaq_re, plaq_im = plaquette_sum(field)
                    plaq_re /= field.NV
                    plaq_im /= field.NV
                    println(measfile, "$itr $plaq_re $plaq_im # plaq_re plaq_im")
                elseif method["methodname"] == "Action"
                    s = field.Sg / field.NV
                    println(measfile, "$itr $s # action")
                elseif method["methodname"] == "Meta_charge" 
                    q = field.CV
                    println(measfile, "$itr $q # metacharge")
                elseif method["methodname"] == "Topological_charge"
                    qt = topological_charge(field)
                    println(measfile, "$itr $qt $(qt^2) # Qtop Qtop^2")
                elseif method["methodname"] == "Polyakov_loop"
                    poly_re, poly_im = poly_loop_avg(field)
                    println(measfile, "$itr $poly_re $poly_im # poly_re poly_im")
                elseif method["methodname"] == "Dirac_eigenvalues"
                    if Du === nothing
                        Du = measset.D(field)
                        Ds = sparse(Du)
                    end

                    vals, _ = eigs(
                        Ds,
                        nev = method["nev"],
                        which = method["which"],
                        maxiter = method["maxiter"],
                    )
                    print(measfile, itr, " ")

                    for λ in vals
                        print(measfile, real(λ), " ", imag(λ), " ")
                    end

                    println(measfile, " # eigs")
                elseif method["methodname"] == "Dirac_determinant"
                    if Du === nothing
                        Du = measset.D(field)
                        Ds = sparse(Du)
                    end

                    logdetD = logdet(Ds)
                    println(
                        measfile,
                        itr, " ", real(logdetD), " ", imag(logdetD), " # log(det(D))"
                    )
                elseif method["methodname"] == "Chiral_condensate"
                    if Du === nothing
                        Du = measset.D(field)
                    end

                    cc = tr(inv(Du)) / field.NV
                    println(
                        measfile,
                        itr, " ", real(cc), " ", imag(cc), " # chiral condensate")
                elseif occursin( "Wilson_loop", method["methodname"])
                    num = filter.(isdigit, method["methodname"])
                    LT = parse(Int64, num)
                    wils_re, wils_im = wilson_loop_all(field, LT)
                    wils_re = round.(wils_re, sigdigits = 4)
                    wils_im = round.(wils_im, sigdigits = 4)
                    println(measfile, "$itr $wils_re $wils_im # wilson_re wilson_im")
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