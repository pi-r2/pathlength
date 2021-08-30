# pathlength

A first approach to pathlength computation.

See the resultats.pdf for an explanation of the reasonning behing the algorithm, aswell as some experimental results.

Compile with $cythonize -a -i pathlength.pyx

Test using sage with 
>>from pathlength import run
>>g = Some_sage_graph
>>pl = run(g)
