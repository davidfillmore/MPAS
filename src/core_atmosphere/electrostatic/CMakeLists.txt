set(ELECTROSTATIC_SOURCES
    mpas_electrostatic_benchmarks.F
    mpas_electrostatic_bcs.F
    mpas_electrostatic_driver.F
    mpas_electrostatic_source.F
    mpas_electrostatic_vertical_remap.F
)

list(TRANSFORM ELECTROSTATIC_SOURCES PREPEND electrostatic/)
