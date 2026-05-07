set(ELECTROSTATIC_SOURCES
    mpas_electrostatic_benchmarks.F
    mpas_electrostatic_bcs.F
    mpas_electrostatic_driver.F
    mpas_electrostatic_source.F
)

list(TRANSFORM ELECTROSTATIC_SOURCES PREPEND electrostatic/)
