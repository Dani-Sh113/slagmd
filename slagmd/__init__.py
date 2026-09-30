"""slagmd: build, run, watch and analyse molecular dynamics of molten oxide slags.

The simulation protocol and the force-field parameters follow the published work of
Verbeeck, Bellemans and co-workers on CaO-Al2O3-SiO2 melts. Every parameter in
forcefields.json carries the printed table it was transcribed from, and every setting a
generated LAMMPS deck writes is tagged [THEIRS], [PAPER] or [MINE] so that a default of
this tool can never be mistaken for a published value.

See README.md for citations.
"""
__version__ = "1.0.0"
