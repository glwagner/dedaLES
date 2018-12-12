import sys; sys.path.append("..")

from mpi4py import MPI

import dedaLES

iters = 100
resolution = {'nx': 64, 'ny': 64, 'nz': 16}

model_dns  = dedaLES.init_rayleigh_benard_benchmark(**resolution, closure=None)
model_smag = dedaLES.init_rayleigh_benard_benchmark(**resolution, closure=dedaLES.ConstantSmagorinsky())

buildtime_dns = dedaLES.benchmark_build(model_dns)
buildtime_smag = dedaLES.benchmark_build(model_smag)

for model in (model_dns, model_smag):
    dedaLES.set_ic_rayleigh_benard_benchmark(model)

dedaLES.mpiprint("Running benchmarks:")
runtime_dns  = dedaLES.benchmark_run(model_dns, iterations=iters)
runtime_smag = dedaLES.benchmark_run(model_smag, iterations=iters)

msg = """

Benchmark results for {iters} iterations:

    ** Build time **
    ================

                   Model     Time
                   -----     ----

                     DNS  |  {build_dns:>8.2f} s
    Constant Smagorinsky  |  {build_sma:>8.2f} s


    ** Run time **
    ==============  

                   Model     Time           Slowdown vs DNS
                   -----     ----           ---------------

                     DNS  |  {run_dns:>8.2f} s  |  1
    Constant Smagorinsky  |  {run_sma:>8.2f} s  |  {sma_slowdown:>2f}


"""
formattedmsg = msg.format(
    iters = iters, 
    build_dns = buildtime_dns,
    build_sma = buildtime_smag,
    run_dns = runtime_dns, 
    run_sma = runtime_smag, 
    sma_slowdown = runtime_smag/runtime_dns,
)

dedaLES.mpiprint(formattedmsg)
