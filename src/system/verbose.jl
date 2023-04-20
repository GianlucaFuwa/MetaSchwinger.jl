module Verbose_print
    import InteractiveUtils

    abstract type Verbose end

    struct Verbose_ <: Verbose
        fp::Union{Nothing,IOStream}
        Verbose_() = new(nothing)
        Verbose_(filename::String) = new(open(filename,"w"))
        Verbose_(fp::IOStream) = new(fp)
    end

    function Base.flush(v::Verbose)
        if v.fp !== nothing
            flush(v.fp)
        end
    end

    function InteractiveUtils.versioninfo(v::Verbose)
        InteractiveUtils.versioninfo()
        if v.fp !== nothing
            InteractiveUtils.versioninfo(v.fp)
        end
    end

    function println_verbose(v::Verbose,val...) 
        println(val...)
        if v.fp !== nothing
            println(v.fp,val...)
        end
    end

    function print2file(v::Verbose,val...)
        if v.fp !== nothing
            println(v.fp,val...)
        end
    end
end