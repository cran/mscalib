\name{A starting point}
\alias{A starting point}
\title{Introductory comments}
\description{
  A brief introduction to the R/qtl package, with a walk-through of an
  analysis.  
}

\section{New to R and/or R/qtl?}{
  \itemize{
    \item In order to use the R/qtl package, you must type (within R)
    \code{library(mscalib)}.  You may wish to include this in a
    \code{\link[base:Startup]{.First}} function or
    \code{\link[base:Startup]{.Rprofile}} file. 

    \item Documention and several tutorials are available
    at the R archive (\url{http://cran.r-project.org}).

    \item Use the \code{\link[base]{help.start}} function to start the
    html version of the R help.  In Windows or MacOS, you may wish to
    use \code{options(htmlhelp=TRUE)} to get easy access to the html
    version of the help files; this could be included in a
    \code{\link[base:Startup]{.First}} function or
    \code{\link[base:Startup]{.Rprofile}} file.

    \item Type \code{library(help=mscalib)} to get a list of the functions
    in R/qtl.

    \item Use the \code{\link[base]{example}} function to run examples
    of the various functions in R/mscalib.

    \item A slide show on the use of R/mscalib is available on 
    \url{http://www.molgen.mpg.de/~wolski/mscalib/}.

    \item Download the latest version of R/mscalib from
    \url{http://www.molgen.mpg.de/~wolski/mscalib/}.
  }
}

\section{Walk-through of an analysis}{ 
Here we briefly describe the use of R/mscalib to analyze,calibrate and filter a set of massvectors.  

A difficult first step in the use of most data analysis software is the
import of data. For loading the Bruker Daltonics File format the package provides the functions
\code{\link[mscalib]{readBruker.massvector},\link[mscalib]{readBruker.massvectorlist}}.
If you have trouble importing data, send an email to Witold Wolski,
\email{wolski@molgen.mpg.de}, perhaps attaching examples of your data
files.  (Such data will be kept confidential.)

We consider the example data \code{\link[mscalib]{mvl}}, a dataset of 333 massvectors
obtained by on run of a sample support on a Bruker Reflex III Instrument.  
Use the \code{\link[base]{data}} function to load the data.

\code{data(mvl)}

The \code{\link[mscalib]{mvl}} data set has class \code{massvectorlist}.  
By visual analysis information about the dataset can be gathered. The \code{\link{hist.massvectorlist}}
functions shows the frequencies of masses in the dataset. The \code{\link{plot.massvectorlist} the masses of all samples
on a stripchart.

\code{hist(mvl)}
\code{plot(mvl)}

The first step of the data

The function \code{\link[qtl]{plot.massvectorlist}} gives a graphical summary of
the data. 


\code{plot(mvl)} \cr
\code{hist(mvl)}



\code{hyper <- drop.nullmarkers(hyper)} \cr
\code{totmar(hyper)}

The function \code{\link[qtl]{est.rf}} estimates the recombination
fraction between each pair of markers, and calculates a LOD score for
the test of \eqn{r} = 1/2.  This is useful for identifying markers that
are placed on the wrong chromosome.  Note that since, for these data,
many markers were typed only on recombinant individuals, the pairwise
recombination fractions show rather odd patterns.

\code{hyper <- est.rf(hyper)} \cr
\code{plot.rf(hyper)} \cr
\code{plot.rf(hyper, chr=c(1,4))}

To re-estimate the genetic map for an experimental cross, use the
function \code{\link[qtl]{est.map}}.  The function
\code{\link[qtl]{plot.map}}, in addition to plotting a single map, can
plot the comparison of two genetic maps (as long as they are composed of
the same numbers of chromosomes and markers per chromosome).  The
function \code{\link[qtl]{replace.map}} map be used to replace the
genetic map in a cross with a new one.

\code{newmap <- est.map(hyper, error.prob=0.01, trace=TRUE)} \cr
\code{plot.map(hyper, newmap)} \cr
\code{hyper <- replace.map(hyper, newmap)}

Before doing QTL analyses, a number of intermediate calculations may
need to be performed.  The function \code{\link[qtl]{calc.genoprob}}
calculates conditional genotype probabilities given the multipoint
marker data.  \code{\link[qtl]{sim.geno}} simulates sequences of
genotypes from their joint distribution, given the observed marker data.
\code{\link[qtl]{argmax.geno}} calculates the most likely sequence of
underlying genotypes, given the observed marker data.  

These three functions return a modified version of the input cross, with
the intermediate calculations included.  

\code{hyper <- calc.genoprob(hyper, step=2.5, error.prob=0.01)} \cr
\code{hyper <- sim.geno(hyper, step=2.5, n.draws=64, error.prob=0.01)} \cr
\code{hyper <- argmax.geno(hyper, error.prob=0.01)}

The function \code{\link[qtl]{calc.errorlod}} may be used to assist in
identifying possible genotyping errors; it calculates the error LOD
scores described by Lincoln and Lander (1992).  It requires the results
of \code{\link[qtl]{calc.genoprob}}, run with a matching value for the
assumed genotyping error probability, \code{error.prob}.

\code{hyper <- calc.errorlod(hyper, error.prob=0.01)} \cr
\code{plot.errorlod(hyper)} \cr
\code{top.errorlod(hyper)} \cr
\code{plot.errorlod(hyper, chr=c(4,11,16))}

The function \code{\link[qtl]{plot.geno}} may be used to inspect the
observed genotypes for a chromosome, with likely genotyping errors
flagged.

\code{plot.geno(hyper, chr=16, ind=71:90, min.sep=4)}

The function \code{\link[qtl]{scanone}} performs a genome scan with a
single QTL model.  By default, it performs standard interval mapping
(Lander and Botstein 1989): use of a normal model and the EM algorithm.
If one specifies \code{method="hk"}, Haley-Knott regression is performed
(Haley and Knott 1992).  These two methods require the results from
\code{\link[qtl]{calc.genoprob}}.

\code{out.em <- scanone(hyper)} \cr
\code{out.hk <- scanone(hyper, method="hk")}

If one specifies \code{method="imp"}, a genome scan is performed by the
multiple imputation method of Sen and Churchill (2001).  This method
requires the results from \code{\link[qtl]{sim.geno}}.

\code{out.imp <- scanone(hyper, method="imp")}

The output of \code{\link[qtl]{scanone}} is a data.frame with class
\code{scanone}.  The function \code{\link[qtl]{plot.scanone}} may be
used to plot the results, and may plot up to three sets of results
against each other, as long as they conform appropriately.

\code{plot(out.em)} \cr
\code{plot(out.hk, col="blue", add=TRUE)} \cr
\code{plot(out.imp, col="red", add=TRUE)} \cr
\code{plot(out.hk, out.imp, out.em, chr=c(1,4), lty=1, col=c("blue","red","black"))}

The function \code{\link[qtl]{summary.scanone}} may be used to list
information on the peak LOD for each chromosome for which the LOD
exceeds a specified threshold.

\code{summary(out.em)} \cr
\code{summary(out.em, 3)} \cr
\code{summary(out.hk, 3)} \cr
\code{summary(out.imp, 3)}

The function \code{\link[qtl]{max.scanone}} returns the maximum LOD
score, genome-wide.

\code{max(out.em)} \cr
\code{max(out.hk)} \cr
\code{max(out.imp)}

One may also use \code{\link[qtl]{scanone}} to perform a permutation
test to get a genome-wide LOD significance threshold.  This will take
some time.

\code{operm.hk <- scanone(hyper, method="hk", n.perm=100)} \cr
\code{quantile(operm.hk, 0.95)}

We should say at this point that the function
\code{\link[base]{save.image}} will save your workspace to disk.  You'll
wish you had used this if R crashes.

\code{save.image()}

The function \code{\link[qtl]{scantwo}} performs a two-dimensional
genome scan with a two-QTL model.  Methods \code{"em"}, \code{"hk"} and
\code{"imp"} are all available.  \code{\link[qtl]{scantwo}} is
considerably slower than \code{\link[qtl]{scanone}}, and can require a
great deal of memory.  Thus, you may wish to create a version of
\code{hyper} for a more coarse grid.

\code{hyper.coarse <- calc.genoprob(hyper, step=10, err=0.01)} \cr
\code{hyper.coarse <- sim.geno(hyper, step=10, n.draws=64, err=0.01)} \cr
\cr
\code{out2.hk <- scantwo(hyper.coarse, method="hk")} \cr
\code{out2.em <- scantwo(hyper.coarse)} \cr
\code{out2.imp <- scantwo(hyper.coarse, method="imp")}

The output is an object with class \code{scantwo}.  The function
\code{\link[qtl]{plot.scantwo}} may be used to plot the results.  The
upper triangle contains LOD scores for tests of epistasis, while the
lower triangle contains joint LOD scores.

\code{plot(out2.hk)} \cr
\code{plot(out2.em)} \cr
\code{plot(out2.imp)}

The function \code{\link[qtl]{summary.scantwo}} lists the interesting
aspects of the output.  You must provide three LOD thresholds: for the joint
LOD, epistasis LOD, and conditional, single-QTL LOD scores.  The locus
pairs giving the highest joint LOD for each pair of chromosomes are
inspected, and those whose LOD score exceed the joint LOD threshold and
for which either the interaction LOD exceeds its threshold or both the
conditional single-QTL LODs exceed their threshold, are printed.

\code{summary(out2.em, c(8, 3, 3))} \cr
\code{summary(out2.em, c(0, 1000, 4))} \cr
\code{summary(out2.em, c(0, 4, 1000))}

The function \code{\link[qtl]{max.scantwo}} returns the maximum joint
and interaction LODs for a two-dimensional genome scan.

\code{max(out2.em)}

Permutation tests may also performed with \code{\link[qtl]{scantwo}};
it may take a few days of CPU time.  The output is a matrix with two
columns: the maximum joint and epistasis LODs, across the
two-dimensional genome scan, for each simulation replicate.

\code{operm2 <- scantwo(hyper.coarse, method="hk", n.perm=100)} \cr
\code{apply(operm2, 2, quantile, 0.95)} \cr
\code{hist(operm.hk,breaks=20)}

}