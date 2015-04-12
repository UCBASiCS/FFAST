[FFAST](https://www.eecs.berkeley.edu/~kannanr/project_ffft.html)
====

Fast Fourier Aliasing-based Sparse Transform (FFAST) for calculating the Discrete Time Fourier (DFT) transform signals having sparse spectrum. 

The distribution contains C++ code for the FFAST engine and example runs.

Idea behind FFAST
====

The algorithm cleverly induces sparse graph alias codes in the DFT domain, by sub-sampling of the time-domain samples. The resulting sparse graph alias codes are then exploited in a peeling decoder. 

Besides the exactly sparse spectrum, the algorithm can handle both the noisy setting and signals with approximately sparse spectra.

Example usage
====

- Display help

        ./ffast --help

- Run an FFAST experiment `(-a)` on randomly generated signals of length 124950 `(-n 124950)` having 10 sparse spectrum `(-k 10)` 30 times `(-i 30)` 

        ./ffast -a -n 124950 -k 10 -i 30

- Run FFAST on input data saved in `inFile.txt` (-f inFile.txt)

        ./ffast -f inFile.txt

References
====
1. [Computing a k-sparse n-length Discrete Fourier Transform using at most 4k samples and O(k logk) complexity -- S. Pawar and K. Ramchandran][ffast1]
2. [A robust FFAST framework for computing a k-sparse n-length DFT in O(k log n) sample complexity using sparse-graph codes -- S. Pawar and K. Ramchandran][ffast2]
3. [Robustifying the Sparse Walsh-Hadamard Transform without Increasing the Sample Complexity of O(K log N) -- X. Li, J.K. Bradley, S. Pawar, K. Ramchandran][ffast3]

[ffast1]: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6620269 
[ffast2]: https://www.eecs.berkeley.edu/~kannanr/assets/project_ffft/ISITExtendedNoisyDFT-v5.pdf
[ffast3]: https://www.eecs.berkeley.edu/~kannanr/assets/project_ffft/WHT_noisy.pdf

Authors
====

Quentin Byron, Thibault Derousseaux, Rohan Varma, Xioa (Simon) Li, Orhan Ocal

