
<!-- README.md is generated from README.Rmd. Please edit that file -->

## **DICEPRO (Deconvolution with Iterative Completion for Estimating cellular Proportion from RNA-seq data)**

Is a software package that has been meticulously designed to address the
persistent challenges that cellular deconvolution methods have faced for
years. One of its key features is its ability to overcome the
limitations posed by incomplete reference matrix, which have hampered
the accuracy and efficiency of many existing methodologies. In doing so,
`DICEPRO` hopes to contribute to the advancement of RNA-seq analysis,
offering a more complete and reliable solution for better elucidating
the complexities of cell composition and gene expression.

`DICEPRO`’s core functionality lies in its iterative deconvolution
process, which continues until convergence is reached. This approach
involves joint updates of the reference matrix, which may contain
missing cell populations, and the cell population abundance matrix,
which contains newly discovered populations. Like other deconvolution
methods, `DICEPRO` starts by fundamentally modelling the deconvolution
problem. It relies on established techniques such as `CiberSortX`, `DCQ`
or others for the supervised aspect, while using
`non-negative matrix factorization (NMF)` decomposition for the
unsupervised component. This hybrid approach enables `DICEPRO` to
provide a robust and versatile solution for unravelling complex RNA-seq
data.

## Installation

``` r
BiocManager::install("kalidouBA/DICEPRO")
```

#### Running cibersortx to deduce the fractions of the different cell populations requires convigurations following the steps below:

`Step 1:`

Install `Docker Desktop` :
(<https://www.docker.com/products/docker-desktop>)

Open Docker Desktop on your local computer and log in. Then, you open a
terminal, and you type the following command:
`\> docker pull cibersortx/fractions`

`Step 2:`

The next thing you need is a token that you will provide every time you
run the CIBERSORTx executables.

You can obtain the token from the CIBERSORTx website:
(<https://cibersortx.stanford.edu/getoken.php>).

Please note that each token is uniquely tied to your user account, and
tokens are good for a specific time interval from date of request, so
you will need to request a new token when an existing one has expired.

`Step 3:`

Once you have pulled the CIBERSORTx executable from Docker, and you have
obtained a token from the CIBERSORTx website, you now have access to the
CIBERSORTx Fractions executable and can run it following the
instructions below.

##### Once the configuration is complete and docker is running, you can run `DICEPRO` with the token from the CIBERSORTx website `cibersortx_email` and `cibersortx_token` using the main function `DICEPRO()`.

##### Other cell deconvolution methods are available to see with the help of the function `running_method()`.
