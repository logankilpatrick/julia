# This file is a part of Julia. License is MIT: https://julialang.org/license

let cmd = `$(Base.julia_cmd()) --depwarn=error --startup-file=no threads_exec.jl`
    for test_nthreads in (1, 2, 8, 8) # run once to try single-threaded mode, then try a couple times to trigger bad races
        run(pipeline(setenv(cmd, "JULIA_NUM_THREADS" => test_nthreads), stdout = stdout, stderr = stderr))
    end
end
