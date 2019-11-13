\name{NEWS}
\title{News for Package \pkg{relations}}
\newcommand{\cpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 0.6-7}{
  \itemize{
    \item \code{relation_dissimilarity()}, \code{relation_consensus()}
    and \code{relation_choice()} gain methods for generalized paired
    comparison dissimilarities (which include symdiff, CKS, and
    Emond-Mason dissimilarities as special cases).

    \item \code{relation_is()} meta-predicate function added, to replace
    the current single predicate functions in future releases.
  
    \item \code{is.na} method for relations added.
  }
}

\section{Changes in version 0.6-6}{
  \itemize{
    \item No user-level changes.
  }
}

\section{Changes in version 0.6-5}{
  \itemize{
    \item Safeguard code against \code{nchar(NA)} being \code{NA}.
  }
}

\section{Changes in version 0.6-4}{
  \itemize{
    \item Rd improvements.
  }
}

\section{Changes in version 0.6-3}{
  \itemize{
    \item Improve namespace and avoid search path semantics.
    \item Improve message strings.
  }
}

\section{Changes in version 0.6-2}{
  \itemize{
    \item Package metadata improvements.
  }
}

\section{Changes in version 0.6-1}{
  \itemize{
    \item Rd improvements.
  }
}

\section{Changes in version 0.6}{
  \itemize{
    \item Code added for computing (weakly/strongly) connected components,
          condensations and component representations.
    \item \code{transitive_reduction()} and the plot method for
    relations now also work for cyclic relations.
    \item The plot method for relations now also allows the drawing of
    \dQuote{raw} graphs without any transformation.
  }
}

\section{Changes in version 0.5-11}{
  \itemize{
    \item Bug in \code{transitive_reduction()} fixed.
  }
}

\section{Changes in version 0.5-10}{
  \itemize{
    \item Vignette fixes.
    \item Package news converted to Rd format.
  }
}

\section{Changes in version 0.5-9}{
  \itemize{
    \item Add predicates for (a)cyclic relations.
    \item Fix incorrect \code{transitive_reduction()} code.
  }
}

\section{Changes in version 0.5-8}{
  \itemize{
    \item Ensure that examples do not fail when package \pkg{Rgraphviz}
    is unavailable.
  }
}

\section{Changes in version 0.5-7}{
  \itemize{
    \item support for finding up to \eqn{n} consensus relations added
    for most methods.
  }
}

\section{Changes in version 0.5-5}{
  \itemize{
    \item \code{mean()} method added for relation ensembles.
    \item Jaccard dissimilarity added.
    \item \code{relation_precursors()} added.
    \item support for missings added.
    \item \code{relation_choice()} gains a method for computing
      Euclidean choices and Schulze Winners.
    \item \code{relation_impute()} added with methods for linear, weak,
      and partial orders.
  }
}
\section{Changes in version 0.5-4}{
  \itemize{
    \item Rd fixes.
  }
}

\section{Changes in version 0.5-3}{
  \itemize{
    \item Adjust for recent changes in package \cpkg{clue}.
  }
}

\section{Changes in version 0.5-2}{
  \itemize{
    \item Sparse lightweight array and matrix code moved to package
      \cpkg{slam}. 
  }
}

\section{Changes in version 0.5-1}{
  \itemize{
    \item Rd fixes.
  }
}

\section{Changes in version 0.5}{
  \itemize{
    \item New \code{relation_consensus()} method for Cook-Kress-Seiford (CKS)
      dissimilarity.
    \item New symdiff fitters for consensus preorders, and family \samp{G}
      for general crisp endorelations.
    \item Performance enhancements for finding all symdiff consensus solutions, and
      symdiff consensus for equivalences or weak orders with a given number of
      classes.
    \item \code{relation_choice()} gains a method for CKS dissimilarity, and a
      \code{weights} argument.
    \item Performance enhancements for symdiff choice via employing a
    more compact 
      linearization of the underlying binary quadratic program (BQP).
    \item GLPK is now used as the default MILP solver.
    \item Consensus and choice fitters now use sparse constraint matrices by
      default.
    \item New class for representing rankings.
    \item Improvements for plotting relations and relation ensembles.
    \item New function \code{relation_cover()} for computing the covering relation of
      an endorelation.
    \item New function \code{relation_trace()} for computing the left or right trace
      of an endorelation.
    \item New simple choice function \code{relation_elements()} for determining the
      minimal, first, last or maximal elements of an endorelation.
    \item New predicate \code{relation_is_quasitransitive()}.
  }
}

\section{Changes in version 0.4}{
  \itemize{
    \item symdiff fitters can now use the GLPK, SYMPHONY and CPLEX solvers in
      addition to lp_solve, and optionally use sparse constraint matrices
    \item more methods added to \code{relation_violations()}.
      optionally, return all tuples violating a property.
    \item \code{relation_scores()} now also works for relation ensembles
    \item Transitive/reflexive closure and reduction functions are now
      available as methods for closure and reduction generics (provided
      by package \cpkg{sets})
    \item During creation of relations, domains are converted to sets and
      incidences are canonicalized to match the internal order of the domain
    \item \code{relation_consensus()} can now compute all solutions for the Condorcet
      method, and gains a \sQuote{majority} method
    \item subassignment method for ensembles added
    \item \code{set_intersection()} accepts any number of arguments
    \item \code{as.relation()} generic now has a \code{...} argument
    \item the \code{attrs} argument in \code{plot.relation_ensemble()} is now a list, so
      that the Rgraphviz parameters can be changed in each plot
      individually
    \item \code{all.equal()} methods added for relations and ensembles
  }
}

\section{Changes in version 0.3}{
  \itemize{
    \item rank method is default for \code{relation_scores()}
    \item the relations in the SVM Benchmark data are now all reflexive
    \item support for sets moved to separate \cpkg{sets} package
    \item support for fuzzy relations added
    \item fitters can now use a sparse constraint matrix representation
    \item new symdiff fitters for the S (\code{"symmetric"}) and M (\code{"matches"}) families
    \item relations are now subsettable
    \item fitters for Cook-Seiford method and Euclidean consensus added
    \item new \code{summary()} method for relations which computes all implemented
      predicates
    \item \code{relation_choice()} for choosing \dQuote{winner} objects based on an ensemble
      of relations between these
  }
}

\section{Changes in version 0.2}{
  \itemize{
    \item Copeland's method added to \code{relation_consensus()}
    \item \code{relation_classes()} added to extract and pretty-print
      (ordered) classes from preferences and equivalences.
    \item \code{relation_violations()} added to compute a measure of remoteness from a
      specified property (e.g., symmetry, transitivity, etc.)
    \item \code{set_symdiff()} fixed for more than two elements (used to compute 
      the complement of the union and the intersection of a set, which is
      wrong in general).  Now, the argument list is treated recursively.
    \item Vignette added
    \item New symdiff fitter added for the C (\code{"complete"}) and A (\code{"antisymmetric"}) 
      families of relations
    \item Summary methods for tuples and sets added
    \item Ops and \code{t()} methods for the relation_ensemble class added
    \item \code{plot.relation()} fixed 
      (used to plot the dual, now changes the graphviz layout)
    \item \code{Ops.tuple()} now correctly handles unary operators, and preserves class
    \item \code{transform_incidences()} added for transforming incidence matrices of
      endorelations from/to other codings found in literature
  }
}

\section{Changes in version 0.1}{
  \itemize{
    \item Initial release
  }
}
