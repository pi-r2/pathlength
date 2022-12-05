# pathlength

A first approach to pathlength computation, and the result of my internship at COATI (https://team.inria.fr/coati/) under the supervision of David Coudert and Nicolas Nisse.

See the resultats.pdf for an explanation of the reasoning behing the algorithm, aswell as some experimental results, or the overleaf link here : https://fr.overleaf.com/read/ydhpwwjrqvys (in french for now)

Compile with $cythonize -a -i pathlength.pyx

Test using sage with 
>from pathlength import run

>g = Some_sage_graph

>pl = run(g)
